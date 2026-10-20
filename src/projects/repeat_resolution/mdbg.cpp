//
// Created by Andrey Bzikadze on 11/25/21.
//

#include "mdbg.hpp"

using namespace repeat_resolution;

std::vector<SuccinctEdgeInfo> MultiplexDBG::SparseDBG2SuccinctEdgeInfo(
    dbg::SparseDBG &dbg, const UniqueClassificator &classificator) {
    const std::unordered_map<std::string, uint64_t> vert2ind = [&dbg]() {
      std::unordered_map<std::string, uint64_t> vert2ind;
      uint64_t cnt;
      for (const dbg::Vertex &vertex : dbg.vertices()) {
          const std::string &id = vertex.getId();
          vert2ind.emplace(id, cnt);
          ++cnt;
      }
      return vert2ind;
    }();

    std::vector<SuccinctEdgeInfo> edge_info;
    for (auto it = dbg.edges().begin(); it!=dbg.edges().end(); ++it) {
        const dbg::Edge &edge = *it;
        const RRVertexType start_ind = vert2ind.at(edge.start()->getId());
        const RRVertexType end_ind = vert2ind.at(edge.end()->getId());
        edge_info.push_back(
            {start_ind, end_ind, &edge, classificator.isUnique(edge)});
    }
    return edge_info;
}

void MultiplexDBG::AssertValidity() const {
    int64_t est_max_vert_index = [this]() {
      int64_t est_max_vert_index{-1};
      for (const auto &vertex : *this) {
          est_max_vert_index = std::max(est_max_vert_index, (int64_t) vertex);
      }
      return est_max_vert_index;
    }();
    VERIFY(next_vert_index >= 1 + est_max_vert_index);

    int64_t est_max_edge_index = [this]() {
      int64_t est_max_edge_index{-1};
      for (const auto &vertex : *this) {
          auto[out_nbr_begin, out_nbr_end] = out_neighbors(vertex);
          for (auto it2 = out_nbr_begin; it2!=out_nbr_end; ++it2) {
              est_max_edge_index =
                  std::max(est_max_edge_index,
                           (int64_t) it2->second.prop().Index());
          }
      }
      return est_max_edge_index;
    }();
    VERIFY(next_edge_index >= 1 + est_max_edge_index);

    for (const auto &vertex : *this) {
        if (count_in_neighbors(vertex)==1 and count_out_neighbors(vertex)==1) {
            auto[in_nbr_begin, in_nbr_end] = in_neighbors(vertex);
            VERIFY_MSG(in_nbr_begin->first==vertex,
                       "No 1in-1out vertices are allowed except loops")
            VERIFY_MSG(node_prop(vertex).IsFrozen(),
                       "An isolated loop must be frozen");
        }
    }

    for (const auto &vertex : *this) {
        auto[out_nbr_begin, out_nbr_end] = out_neighbors(vertex);
        const RRVertexProperty &vertex_prop = node_prop(vertex);
        for (auto it = out_nbr_begin; it!=out_nbr_end; ++it) {
            const RREdgeProperty &edge_prop = it->second.prop();
            const int64_t inner_edge_size = edge_prop.Size();
            if (inner_edge_size < 0) {
                const RRVertexProperty &neighbor_prop = node_prop(it->first);
                VERIFY(vertex_prop.Seq().ToSequence().Suffix(-inner_edge_size)==
                    neighbor_prop.Seq().ToSequence().Prefix(-inner_edge_size));
            }
        }
    }

    for (const RRVertexType &vertex : *this) {
        const RRVertexProperty &vertex_prop = node_prop(vertex);
        if (not vertex_prop.IsFrozen()) {
            VERIFY(vertex_prop.size()==n_iter + start_k);
        }
    }

    if (contains_rc) {
        std::map<Sequence, RREdgeIndexType> seq_edge;
        std::unordered_map<RREdgeIndexType, bool> is_unique_edge;
        for (const RRVertexType &vertex : *this) {
            auto[begin, end] = out_neighbors(vertex);
            for (auto it = begin; it!=end; ++it) {
                const RREdgeProperty &edge_prop = it->second.prop();
                seq_edge.emplace(
                    GetEdgeSequence(find(vertex), it, false, false)
                        .ToSequence(),
                        edge_prop.Index());
                is_unique_edge.emplace(edge_prop.Index(), edge_prop.IsUnique());
            }
        }

        for (const RRVertexType &vertex : *this) {
            auto[begin, end] = out_neighbors(vertex);
            for (auto it = begin; it!=end; ++it) {
                const RREdgeProperty &edge_prop = it->second.prop();
                const Sequence seq =
                    GetEdgeSequence(find(vertex), it, false, false)
                        .ToSequence();
                VERIFY_MSG(seq_edge.find(!seq)!=seq_edge.end(),
                           "no rev comp for edge " + itos(edge_prop.Index()));
                VERIFY_MSG(is_unique_edge.at(edge_prop.Index())
                               ==is_unique_edge.at(seq_edge.at(!seq)),
                           "edge_prop " + itos(edge_prop.Index()) + ", unique: "
                               + itos(is_unique_edge.at(edge_prop.Index()))
                               + " . rev compl uniqueness not equal");
            }
        }
    }
}

void MultiplexDBG::SpreadFrost() {
    std::unordered_set<RRVertexType> prev_frozen, new_frozen;
    for (const RRVertexType &vertex : *this) {
        const RRVertexProperty &vertex_prop = node_prop(vertex);
        if (vertex_prop.IsFrozen()) {
            prev_frozen.insert(vertex);
        }
    }

    auto upd_new_frozen = [this, &new_frozen](const RRVertexType &vertex,
                                              const RRVertexProperty &vertex_prop,
                                              NeighborsIterator begin,
                                              NeighborsIterator end) {
      for (auto it = begin; it!=end; ++it) {
          const RREdgeProperty &edge_prop = it->second.prop();
          const RRVertexType &neighbor = it->first;
          const RRVertexProperty &neighbor_prop = node_prop(neighbor);
          if (not neighbor_prop.IsFrozen() and
              FullEdgeSize(find(vertex), it)==1 + neighbor_prop.size()) {
              FreezeVertex(neighbor);
              new_frozen.insert(neighbor);
          }
      }
    };

    while (not prev_frozen.empty()) {
        for (const RRVertexType &vertex : prev_frozen) {
            auto[in_nbr_begin, in_nbr_end] = in_neighbors(vertex);
            auto[out_nbr_begin, out_nbr_end] = out_neighbors(vertex);
            upd_new_frozen(vertex, node_prop(vertex), in_nbr_begin, in_nbr_end);
            upd_new_frozen(vertex,
                           node_prop(vertex),
                           out_nbr_begin,
                           out_nbr_end);
        }
        prev_frozen = std::move(new_frozen);
    }
}

void MultiplexDBG::FreezeUnpairedVertices() {
    for (const RRVertexType &vertex : *this) {
        RRVertexProperty &vertex_prop = node_prop(vertex);
        if (vertex_prop.IsFrozen()) {
            continue;
        }

        auto[in_edges, out_edges] = GetNeighborEdgesIndexes(vertex);
        if (in_edges.size()==1 and out_edges.size()==1) {
            // must be a self-loop
            VERIFY(in_edges==out_edges);
            FreezeVertex(vertex);
        } else if (in_edges.size() >= 2 and out_edges.size() >= 2) {
            auto[ac_s2e, ac_e2s] = GetEdgepairsVertex(vertex);
            for (const RREdgeIndexType &edge : in_edges) {
                if (ac_s2e.find(edge)==ac_s2e.end()) {
                    FreezeVertex(vertex);
                    break;
                }
            }
            for (const RREdgeIndexType &edge : out_edges) {
                if (ac_e2s.find(edge)==ac_e2s.end()) {
                    FreezeVertex(vertex);
                    break;
                }
            }
        }
    }
}

std::unordered_map<RREdgeIndexType, Sequence>
MultiplexDBG::GetEdgeSeqs(size_t threads) const {
    ParallelRecordCollector<std::pair<RREdgeIndexType, Sequence>> res(threads);
    std::vector<ConstIterator> vits;
    for (auto v_it = begin(); v_it!=end(); ++v_it) {
        vits.emplace_back(v_it);
    }
    omp_set_num_threads(threads);
#pragma omp parallel for default(none) shared(vits, res)
    for (size_t i = 0; i < vits.size(); i++) {
        ConstIterator v_it = vits[i];
        auto[e_begin, e_end] = out_neighbors(v_it);
        for (auto e_it = e_begin; e_it!=e_end; ++e_it) {
            const RREdgeProperty &prop = e_it->second.prop();
            res.emplace_back(prop.Index(),
                         GetEdgeSequence(v_it, e_it, false, false)
                         .ToSequence());
        }
    }
    return {res.begin(), res.end()};
}

std::unordered_map<RRVertexType, Sequence> MultiplexDBG::GetVertexSeqs(
    const std::unordered_map<RREdgeIndexType, Sequence> &edge_seq) const {
    std::unordered_map<RRVertexType, Sequence> seqs;
    for (auto it = begin(); it!=end(); ++it) {
        const RRVertexProperty &vertex_prop = node_prop(it);
        const uint64_t vertex_size = vertex_prop.size();
        auto[ibegin, iend] = in_neighbors(it);
        auto[obegin, oend] = out_neighbors(it);
        if (ibegin!=iend) {
            const RREdgeIndexType edge_ind = ibegin->second.prop().Index();
            seqs.emplace(*it, edge_seq.at(edge_ind).Suffix(vertex_size));
        } else if (obegin!=oend) {
            const RREdgeIndexType edge_ind = obegin->second.prop().Index();
            seqs.emplace(*it, edge_seq.at(edge_ind).Prefix(vertex_size));
        } else {
            seqs.emplace(*it, vertex_prop.Seq().ToSequence());
        }
    }
    return seqs;
}

void MultiplexDBG::MoveEdge(const RRVertexType &s1, NeighborsIterator e1_it,
                            const RRVertexType &s2, const RRVertexType &e2) {
    // this method by itself does not update read paths
    add_edge_with_prop(s2, e2, std::move(e1_it->second.prop()));
    ConstIterator s1_it = find(s1);
    remove_edge(s1_it, e1_it);
}

void MultiplexDBG::MergeEdges(const RRVertexType &s1, NeighborsIterator e1_it,
                              NeighborsIterator e2_it) {
    const RRVertexType &s2 = e1_it->first;
    VERIFY_MSG(not node_prop(s2).IsFrozen(),
               "Cannot merge edges via a frozen vertex");
    RREdgeProperty &e1_prop = e1_it->second.prop();
    RREdgeProperty &e2_prop = e2_it->second.prop();
    const RREdgeIndexType e2_index = e2_prop.Index();
    rr_paths->Merge(e1_prop.Index(), e2_prop.Index());
    const RRVertexProperty &v1 = node_prop(s1);
    const RRVertexProperty &v3 = node_prop(e2_it->first);
    e1_prop.Merge(std::move(node_prop(s2)), std::move(e2_prop));
    MoveEdge(s1, e1_it, s1, e2_it->first);
    remove_edge(find(s2), FindOutEdgeIterator(s2, e2_index));
    remove_nodes(s2);
}

RREdgeIndexType MultiplexDBG::AddConnectingEdge(NeighborsIterator eleft_it,
                                                const RRVertexType &vright,
                                                NeighborsIterator eright_it) {
    const RRVertexType &vleft = eleft_it->first;
    VERIFY_MSG(vleft!=vright, "Can only add edge b/w disconnected edges");
    const RRVertexProperty &vleft_prop = node_prop(vleft);
    const RRVertexProperty &vright_prop = node_prop(vright);

    const RREdgeProperty &eleft_prop = eleft_it->second.prop();
    const RREdgeProperty &eright_prop = eright_it->second.prop();
    const RREdgeIndexType new_index = next_edge_index;
    ++next_edge_index;

    RREdgeProperty e_new_prop = Add(vleft_prop, vright_prop, new_index);
    rr_paths->Add(eleft_prop.Index(), eright_prop.Index(), e_new_prop.Index());
    add_edge_with_prop(vleft, vright, std::move(e_new_prop));
    return new_index;
}

RRVertexType MultiplexDBG::GetNewVertex(MDBGSeq seq) {
    RRVertexType new_vertex{next_vert_index};
    ++next_vert_index;
    RRVertexProperty property(std::move(seq), false);
    add_node_with_prop(new_vertex, std::move(property));
    return new_vertex;
}

MultiplexDBG::MultiplexDBG(const std::vector<SuccinctEdgeInfo> &edges,
                           const uint64_t start_k, RRPaths *const rr_paths,
                           bool contains_rc)
    : rr_paths{rr_paths}, start_k{start_k}, contains_rc{contains_rc} {
    for (const SuccinctEdgeInfo &edge_info : edges) {
        const dbg::Edge *edge = edge_info.edge;
        next_vert_index = std::max(next_vert_index, 1 + edge_info.start_ind);
        next_vert_index = std::max(next_vert_index, 1 + edge_info.end_ind);
        add_node_with_prop(edge_info.start_ind,
                           RRVertexProperty(MDBGSeq(edge, 0, start_k), false));
        add_node_with_prop(
            edge_info.end_ind,
            // Anton's edge does not contain prefix
            RRVertexProperty(MDBGSeq(edge,
                                     edge->size(),
                                     edge->size() + start_k),
                             false));

        int64_t infix_size = ((int64_t) edge->size()) - start_k;
        VERIFY(infix_size > 0 or -infix_size < start_k);
        MDBGSeq edge_seq;
        if (infix_size > 0) {
            edge_seq = MDBGSeq(edge, start_k, start_k + infix_size);
        }
        RREdgeProperty edge_property(next_edge_index,
                                     std::move(edge_seq),
                                     infix_size,
                                     edge_info.unique);
        add_edge_with_prop(edge_info.start_ind, edge_info.end_ind,
                           std::move(edge_property));
        ++next_edge_index;
    }

    FreezeUnpairedVertices();
    SpreadFrost();
}

MultiplexDBG::MultiplexDBG(dbg::SparseDBG &dbg, RRPaths *const rr_paths,
                           const uint64_t start_k,
                           UniqueClassificator &classificator)
    : MultiplexDBG(SparseDBG2SuccinctEdgeInfo(dbg, classificator), start_k,
                   rr_paths, true) {}

void MultiplexDBG::ExportToDot(
    const std::experimental::filesystem::path &path) const {
    graph_lite::Serializer serializer(*this);
    serializer.set_max_num_nodes_per_line(1);
    serializer.set_max_num_edges_per_line(1);
    std::ofstream dot_os(path);
    serializer.serialize_to_dot(dot_os);
}

void MultiplexDBG::ExportToGFA(
    const std::experimental::filesystem::path &path, size_t threads) const {
    const std::unordered_map<RREdgeIndexType, Sequence>
        edge_seqs = GetEdgeSeqs(threads);
    const std::unordered_map<RRVertexType, Sequence> vertex_seqs =
        GetVertexSeqs(edge_seqs);

    const std::unordered_map<RRVertexType, RRVertexType> vertex2rc =
        MapSeqs2RC<RRVertexType>(vertex_seqs);
    const std::unordered_map<RREdgeIndexType, RREdgeIndexType> edge2rc =
        MapSeqs2RC<RREdgeIndexType>(edge_seqs);

    const std::unordered_map<RRVertexType, bool> vertex_can =
        AreSeqsCanonical<RRVertexType>(vertex_seqs);
    const std::unordered_map<RREdgeIndexType, bool> edge_can =
        AreSeqsCanonical<RREdgeIndexType>(edge_seqs);

    ExportToGFA(path, vertex_seqs, edge_seqs, vertex2rc, edge2rc, vertex_can,
                edge_can);
}

void MultiplexDBG::ExportToGFA(
    const std::experimental::filesystem::path &path,
    const std::unordered_map<RRVertexType, Sequence> &vertex_seqs,
    const std::unordered_map<RREdgeIndexType, Sequence> &edge_seqs,
    const std::unordered_map<RRVertexType, RRVertexType> &vertex2rc,
    const std::unordered_map<RREdgeIndexType, RREdgeIndexType> &edge2rc,
    const std::unordered_map<RRVertexType, bool> &vertex_can,
    const std::unordered_map<RREdgeIndexType, bool> &edge_can) const {

    std::ofstream os;
    os.open(path);
    os << "H\tVN:Z:1.0" << std::endl;
    // edge2can_id: 任意边索引 -> 规范化边索引（canonical edge id）
    std::unordered_map<RREdgeIndexType, RREdgeIndexType> edge2can_id;
    // 导出一致性统计：用于衡量因映射缺失被跳过的链接规模
    size_t canonical_edge_cnt = 0;
    size_t missing_edge_can = 0;
    size_t missing_edge2rc = 0;
    size_t missing_edge_seq = 0;
    size_t missing_vertex_can = 0;
    size_t missing_vertex2rc = 0;
    size_t missing_out_edge_can = 0;
    size_t missing_in_edge_can = 0;
    size_t missing_out_can = 0;
    size_t missing_in_can = 0;
    size_t candidate_link_cnt = 0;
    size_t emitted_link_cnt = 0;
    for (auto v_it = begin(); v_it!=end(); ++v_it) {
        auto[begin, end] = out_neighbors(v_it);
        for (auto e_it = begin; e_it!=end; ++e_it) {
            const RREdgeProperty &prop = e_it->second.prop();
            const RREdgeIndexType e_ind = prop.Index();
            // 防御性访问：edge_can/edge2rc/edge_seq 任一缺失都统计并跳过当前边
            auto ecan_it = edge_can.find(e_ind);
            if (ecan_it==edge_can.end()) {
                missing_edge_can += 1;
                continue;
            }
            if (ecan_it->second) {
                auto erc_it = edge2rc.find(e_ind);
                if (erc_it==edge2rc.end()) {
                    missing_edge2rc += 1;
                    continue;
                }
                auto eseq_it = edge_seqs.find(e_ind);
                if (eseq_it==edge_seqs.end()) {
                    missing_edge_seq += 1;
                    continue;
                }
                canonical_edge_cnt += 1;
                edge2can_id.emplace(e_ind, e_ind);
                edge2can_id.emplace(erc_it->second, e_ind);
                os << "S\t" << e_ind << "\t" << eseq_it->second << "\n";
            }
        }
    }

    for (auto v_it = begin(); v_it!=end(); ++v_it) {
        // 防御性访问：若 vertex_can 缺失该点，按非 canonical 处理并跳过
        auto vcan_it = vertex_can.find(*v_it);
        if (vcan_it==vertex_can.end()) {
            missing_vertex_can += 1;
            continue;
        }
        if (not vcan_it->second) {
            continue;
        }
        auto[begin, end] = out_neighbors(v_it);
        for (auto out_it = begin; out_it!=end; ++out_it) {
            const RREdgeProperty &out_prop = out_it->second.prop();
            const RREdgeIndexType out_ind = out_prop.Index();
            // 防御性访问：若 edge_can 缺失，默认按负链路符号处理
            auto out_sign_it = edge_can.find(out_ind);
            if (out_sign_it==edge_can.end()) {
                missing_out_edge_can += 1;
            }
            bool out_sign = out_sign_it==edge_can.end() ? false : out_sign_it->second;
            auto out_can_it = edge2can_id.find(out_ind);
            if (out_can_it==edge2can_id.end()) {
                missing_out_can += 1;
                continue;
            }
            const RREdgeIndexType out_can_ind = out_can_it->second;
            // 防御性访问：若缺少 RC 顶点映射，则跳过该顶点的链接输出
            auto vrc_it = vertex2rc.find(*v_it);
            if (vrc_it==vertex2rc.end()) {
                missing_vertex2rc += 1;
                continue;
            }
            const RRVertexType v_rc = vrc_it->second;
            auto[begin_rc, end_rc] = out_neighbors(v_rc);

            for (auto in_it = begin_rc; in_it!=end_rc; ++in_it) {
                candidate_link_cnt += 1;
                const RREdgeProperty &in_prop = in_it->second.prop();
                const RREdgeIndexType in_ind = in_prop.Index();
                auto in_sign_it = edge_can.find(in_ind);
                if (in_sign_it==edge_can.end()) {
                    missing_in_edge_can += 1;
                }
                bool in_sign = in_sign_it==edge_can.end() ? true : (not in_sign_it->second);
                auto in_can_it = edge2can_id.find(in_ind);
                if (in_can_it==edge2can_id.end()) {
                    missing_in_can += 1;
                    continue;
                }
                const RREdgeIndexType in_can_ind = in_can_it->second;
                os << "L\t" << in_can_ind << "\t" << (in_sign ? "+" : "-")
                   << "\t"
                   << out_can_ind << "\t" << (out_sign ? "+" : "-") << "\t"
                   << node_prop(v_it).size() << "M\n";
                emitted_link_cnt += 1;
            }
        }
    }
    const size_t skipped_link_cnt = missing_out_can + missing_in_can;
    // 写出一致性报告，便于外部脚本或人工排查导出质量
    std::ofstream report(path.string() + ".consistency_report.txt");
    report << "canonical_edges\t" << canonical_edge_cnt << "\n";
    report << "missing_edge_can\t" << missing_edge_can << "\n";
    report << "missing_edge2rc\t" << missing_edge2rc << "\n";
    report << "missing_edge_seq\t" << missing_edge_seq << "\n";
    report << "missing_vertex_can\t" << missing_vertex_can << "\n";
    report << "missing_vertex2rc\t" << missing_vertex2rc << "\n";
    report << "missing_out_edge_can\t" << missing_out_edge_can << "\n";
    report << "missing_in_edge_can\t" << missing_in_edge_can << "\n";
    report << "candidate_links\t" << candidate_link_cnt << "\n";
    report << "emitted_links\t" << emitted_link_cnt << "\n";
    report << "skipped_links_total\t" << skipped_link_cnt << "\n";
    report << "missing_out_can\t" << missing_out_can << "\n";
    report << "missing_in_can\t" << missing_in_can << "\n";
    report.close();
    std::cout << "GFA consistency: skipped_links=" << skipped_link_cnt
              << " (missing_out_can=" << missing_out_can
              << ", missing_in_can=" << missing_in_can << ")" << std::endl;
}

[[nodiscard]] bool MultiplexDBG::IsFrozen() const {
    return std::all_of(begin(), end(), [this](const RRVertexType &v) {
      return node_prop(v).IsFrozen();
    });
}

std::vector<RREdgeIndexType>
MultiplexDBG::GetInEdgesIndexes(const RRVertexType &vertex) const {
    std::vector<RREdgeIndexType> indexes;
    auto[in_nbr_begin, in_nbr_end] = in_neighbors(vertex);
    for (auto it = in_nbr_begin; it!=in_nbr_end; ++it) {
        indexes.push_back(it->second.prop().Index());
    }
    return indexes;
}

bool MultiplexDBG::IsVertexComplex(const RRVertexType &vertex) const {
    const int indegree = count_in_neighbors(vertex);
    const int outdegree = count_out_neighbors(vertex);
    return indegree >= 2 and outdegree >= 2;
}

bool MultiplexDBG::IsVertexSimple(const RRVertexType &vertex) const {
    return not IsVertexComplex(vertex);
}

bool MultiplexDBG::IsVertexCanonical(const RRVertexType &vertex) const {
    const RRVertexProperty &vertex_prop = node_prop(vertex);
    return vertex_prop.IsCanonical();
}

bool MultiplexDBG::IsEdgeCanonical(ConstIterator vertex,
                                   NeighborsConstIterator e_it) const {
    const RREdgeProperty &edge_prop = e_it->second.prop();
    const MDBGSeq edge_seq = GetEdgeSequence(vertex, e_it, false, false);
    return edge_seq.IsCanonical();
}

size_t MultiplexDBG::FullEdgeSize(ConstIterator st_v_it,
                                  NeighborsConstIterator e_it) const {
    const RRVertexType &st_v = *st_v_it;
    const RRVertexType &en_v = e_it->first;
    const RRVertexProperty &st_v_prop = node_prop(st_v_it);
    const RRVertexProperty &en_v_prop = node_prop(en_v);
    const RREdgeProperty &edge_prop = e_it->second.prop();
    int64_t inner_edge_size = edge_prop.Size();
    if (inner_edge_size < 0) {
        VERIFY(st_v_prop.size() >= -inner_edge_size);
        VERIFY(en_v_prop.size() >= -inner_edge_size);
    }
    return st_v_prop.size() + inner_edge_size + en_v_prop.size();
}

MDBGSeq MultiplexDBG::ExtractEdgePostStartPrefix(ConstIterator st_v_it,
                                                 NeighborsIterator e_it,
                                                 uint64_t len) {
    const RRVertexType &st_v = *st_v_it;
    const RRVertexType &en_v = e_it->first;
    const RRVertexProperty &st_v_prop = node_prop(st_v);
    const RRVertexProperty &en_v_prop = node_prop(en_v);
    RREdgeProperty &edge_prop = e_it->second.prop();
    VERIFY(len + st_v_prop.size() <= FullEdgeSize(st_v_it, e_it));

    uint64_t inner_part_len =
        std::min(len, (uint64_t) std::max(0L, edge_prop.Size()));
    MDBGSeq prefix = edge_prop.ExtractSeqPrefix(inner_part_len);

    uint64_t en_v_part_len = len - inner_part_len;
    VERIFY(en_v_part_len <= en_v_prop.size());
    prefix.Append(en_v_prop.GetSeqPrefix(en_v_part_len, -edge_prop.Size()));
    if (en_v_part_len) {
        edge_prop.ShortenWithEmptySeq(en_v_part_len);
    }
    return prefix;
}

MDBGSeq MultiplexDBG::ExtractEdgePreEndSuffix(ConstIterator en_v_it,
                                              NeighborsIterator e_it,
                                              uint64_t len) {
    // Edge is reversed
    const RRVertexType &st_v = e_it->first;
    const RRVertexType &en_v = *en_v_it;
    const RRVertexProperty &st_v_prop = node_prop(st_v);
    const RRVertexProperty &en_v_prop = node_prop(en_v);
    RREdgeProperty &edge_prop = e_it->second.prop();
    size_t full_edge_size = FullEdgeSize(find(en_v), e_it);
    VERIFY(len + en_v_prop.size() <= full_edge_size);

    uint64_t inner_part_len =
        std::min(len, (uint64_t) std::max(0L, edge_prop.Size()));
    uint64_t st_v_part_len = len - inner_part_len;
    VERIFY(st_v_part_len <= st_v_prop.size());
    MDBGSeq suffix = st_v_prop.GetSeqSuffix(st_v_part_len, -edge_prop.Size());
    suffix.Append(edge_prop.ExtractSeqSuffix(inner_part_len));
    if (st_v_part_len) {
        edge_prop.ShortenWithEmptySeq(st_v_part_len);
    }
    return suffix;
}

void MultiplexDBG::IncreaseVertex(const RRVertexType &vertex, uint64_t len) {
    const int indegree = count_in_neighbors(vertex);
    const int outdegree = count_out_neighbors(vertex);
    VERIFY((indegree==1)!=(outdegree==1));
    if (indegree==1) {
        NeighborsIterator edge_it = in_neighbors(vertex).first;
        MDBGSeq new_seq = ExtractEdgePreEndSuffix(find(vertex), edge_it, len);
        node_prop(vertex).IncLeft(std::move(new_seq));
    } else {
        VERIFY(outdegree==1);
        NeighborsIterator edge_it = out_neighbors(vertex).first;
        MDBGSeq
            new_seq = ExtractEdgePostStartPrefix(find(vertex), edge_it, len);
        node_prop(vertex).IncRight(std::move(new_seq));
    }
}

std::vector<RREdgeIndexType>
MultiplexDBG::GetOutEdgesIndexes(const RRVertexType &vertex) const {
    std::vector<RREdgeIndexType> indexes;
    auto[out_nbr_begin, out_nbr_end] = out_neighbors(vertex);
    for (auto it = out_nbr_begin; it!=out_nbr_end; ++it) {
        indexes.push_back(it->second.prop().Index());
    }
    return indexes;
}

std::pair<std::vector<RREdgeIndexType>, std::vector<RREdgeIndexType>>
MultiplexDBG::GetNeighborEdgesIndexes(const RRVertexType &vertex) const {
    return {GetInEdgesIndexes(vertex), GetOutEdgesIndexes(vertex)};
}

std::pair<MultiplexDBG::EdgeNeighborMap, MultiplexDBG::EdgeNeighborMap>
MultiplexDBG::GetEdgepairsVertex(const RRVertexType &vertex) const {
    auto get_init_transitions =
        [this](const std::vector<RREdgeIndexType> &in_edges,
               const std::vector<RREdgeIndexType> &out_edges) {
          EdgeNeighborMap ac_s2e, ac_e2s;
          for (const RREdgeIndexType &in_ind : in_edges) {
              for (const RREdgeIndexType &out_ind : out_edges) {
                  if (rr_paths->ContainsPair(in_ind, out_ind)) {
                      ac_s2e[in_ind].emplace(out_ind);
                      ac_e2s[out_ind].emplace(in_ind);
                  }
              }
          }
          return std::make_pair(ac_s2e, ac_e2s);
        };

    auto extend_transitions_single_loop =
        [this, &vertex](const std::vector<RREdgeIndexType> &in_edges,
                        const std::vector<RREdgeIndexType> &out_edges,
                        EdgeNeighborMap &ac_s2e, EdgeNeighborMap &ac_e2s) {
          std::vector<RREdgeIndexType> loops;
          for (const RREdgeIndexType &index : in_edges) {
              if (std::find(out_edges.begin(), out_edges.end(), index)!=
                  out_edges.end()) {
                  loops.push_back(index);
              }
          }

          if (loops.size()==1) {
              const RREdgeIndexType loop = loops.front();
              const RREdgeProperty &loop_prop =
                  FindOutEdgeConstiterator(vertex, loop)->second.prop();
              if (loop_prop.IsUnique()) {
                  if (in_edges.size()==2) {
                      const size_t loop_index = in_edges.back()==loop;
                      const size_t nonloop = in_edges[loop_index ^ 1];
                      ac_s2e[nonloop].emplace(loop);
                      ac_e2s[loop].emplace(nonloop);
                  }
                  if (out_edges.size()==2) {
                      const size_t loop_index = out_edges.back()==loop;
                      const size_t nonloop = out_edges[loop_index ^ 1];
                      ac_s2e[loop].emplace(nonloop);
                      ac_e2s[nonloop].emplace(loop);
                  }
              }
          }
        };

    auto extend_transitions_all_unique =
        [this, &vertex](const std::vector<RREdgeIndexType> &in_edges,
                        const std::vector<RREdgeIndexType> &out_edges,
                        EdgeNeighborMap &ac_s2e, EdgeNeighborMap &ac_e2s) {
          std::vector<RREdgeIndexType> unpaired_in, unpaired_out;
          for (const RREdgeIndexType &index : out_edges) {
              if (ac_e2s.find(index)==ac_e2s.end()) {
                  unpaired_out.push_back(index);
              }
          }
          for (const RREdgeIndexType &index : in_edges) {
              if (ac_s2e.find(index)==ac_s2e.end()) {
                  unpaired_in.push_back(index);
              }
          }

          bool all_in_unique =
              std::all_of(in_edges.begin(), in_edges.end(),
                          [this, &vertex](const RREdgeIndexType &edge_index) {
                            return FindInEdgeConstiterator(vertex, edge_index)
                                ->second.prop()
                                .IsUnique();
                          });
          bool all_out_unique =
              std::all_of(out_edges.begin(), out_edges.end(),
                          [this, &vertex](const RREdgeIndexType &edge_index) {
                            return FindOutEdgeConstiterator(vertex, edge_index)
                                ->second.prop()
                                .IsUnique();
                          });
          if (unpaired_in.size()==1 and unpaired_out.size()==1 and
              (all_in_unique or all_out_unique)) {
              const RRVertexType unp_in = unpaired_in.front();
              const RRVertexType unp_out = unpaired_out.front();
              VERIFY(ac_s2e.find(unp_in)==ac_s2e.end());
              VERIFY(ac_e2s.find(unp_out)==ac_e2s.end());
              ac_s2e[unp_in] = {unp_out};
              ac_e2s[unp_out] = {unp_in};
          }
        };

    const auto[in_edges, out_edges] = GetNeighborEdgesIndexes(vertex);
    auto[ac_s2e, ac_e2s] = get_init_transitions(in_edges, out_edges);
    extend_transitions_single_loop(in_edges, out_edges, ac_s2e, ac_e2s);
    extend_transitions_all_unique(in_edges, out_edges, ac_s2e, ac_e2s);

    return std::make_pair(ac_s2e, ac_e2s);
}

MultiplexDBG::NeighborsIterator
MultiplexDBG::FindInEdgeIterator(const RRVertexType &v,
                                 const RREdgeIndexType &edge) {
    auto[it, end] = in_neighbors(v);
    while (it!=end and it->second.prop().Index()!=edge) {
        ++it;
    }
    return it;
};

MultiplexDBG::NeighborsConstIterator
MultiplexDBG::FindInEdgeConstiterator(const RRVertexType &v,
                                      const RREdgeIndexType &edge) const {
    auto[it, end] = in_neighbors(v);
    while (it!=end and it->second.prop().Index()!=edge) {
        ++it;
    }
    return it;
};

MultiplexDBG::NeighborsIterator
MultiplexDBG::FindOutEdgeIterator(const RRVertexType &v,
                                  const RREdgeIndexType &edge) {
    auto[it, end] = out_neighbors(v);
    while (it!=end and it->second.prop().Index()!=edge) {
        ++it;
    }
    return it;
};

MultiplexDBG::NeighborsConstIterator
MultiplexDBG::FindOutEdgeConstiterator(const RRVertexType &v,
                                       const RREdgeIndexType &edge) const {
    auto[it, end] = out_neighbors(v);
    while (it!=end and it->second.prop().Index()!=edge) {
        ++it;
    }
    return it;
}

int64_t MultiplexDBG::GetInnerEdgeSize(ConstIterator vertex,
                                       NeighborsConstIterator e_it) const {
    return e_it->second.prop().Size();
}

MDBGSeq MultiplexDBG::GetEdgeSequence(ConstIterator vertex,
                                      NeighborsConstIterator e_it,
                                      bool trim_left, bool trim_right) const {
    const RRVertexProperty &vertex_prop = node_prop(vertex);
    const RRVertexProperty &neighbor_prop = node_prop(e_it->first);
    const RREdgeProperty &edge_prop = e_it->second.prop();

    if (edge_prop.Size() < 0) {
        if (trim_left and trim_right) {
            return {};
        }

        MDBGSeq neighbor_seq = neighbor_prop.Seq();
        if (trim_left) {
            neighbor_seq.TrimLeft(-edge_prop.Size());
            return neighbor_seq;
        }

        MDBGSeq vertex_seq = vertex_prop.Seq();
        vertex_seq.TrimRight(-edge_prop.Size());
        if (not trim_right) {
            vertex_seq.Append(neighbor_seq);
        }
        return vertex_seq;
    }

    MDBGSeq seq = edge_prop.Seq();
    if (not trim_left) {
        MDBGSeq vertex_seq = vertex_prop.Seq();
        seq.Prepend(std::move(vertex_seq));
    }
    if (not trim_right) {
        MDBGSeq neighbor_seq = neighbor_prop.Seq();
        seq.Append(std::move(neighbor_seq));
    }
    return seq;
}

std::vector<Contig>
MultiplexDBG::GetContigs(size_t threads) const {
    const std::unordered_map<RREdgeIndexType, Sequence>
        edge_seqs = GetEdgeSeqs(threads);
    const std::unordered_map<RRVertexType, Sequence> vertex_seqs =
        GetVertexSeqs(edge_seqs);

    const std::unordered_map<RRVertexType, RRVertexType> vertex2rc =
        MapSeqs2RC<RRVertexType>(vertex_seqs);

    std::unordered_map<RRVertexType, bool> vertex_can =
        AreSeqsCanonical<RRVertexType>(vertex_seqs);
    std::unordered_map<RREdgeIndexType, bool> edge_can =
        AreSeqsCanonical<RREdgeIndexType>(edge_seqs);

    return GetContigs(vertex_seqs, edge_seqs, vertex2rc, vertex_can, edge_can);
}

[[nodiscard]] std::vector<Contig> MultiplexDBG::GetContigs(
    const std::unordered_map<RRVertexType, Sequence> &vertex_seqs,
    const std::unordered_map<RREdgeIndexType, Sequence> &edge_seqs,
    const std::unordered_map<RRVertexType, RRVertexType> &vertex2rc,
    const std::unordered_map<RRVertexType, bool> &vertex_can,
    const std::unordered_map<RREdgeIndexType, bool> &edge_can) const {

    const std::unordered_map<RRVertexType, bool> trim = [this, &vertex_can,
        &vertex2rc]() {
      std::unordered_map<RRVertexType, bool> trim;
      for (const RRVertexType &vertex : *this) {
          const RRVertexProperty &vertex_prop = node_prop(vertex);
          // 防御性访问：顶点 canonical 标记缺失时，默认不参与 trim 计算
          auto vcan_it = vertex_can.find(vertex);
          if (vcan_it!=vertex_can.end() && vcan_it->second) {
              const bool trim_vertex = count_out_neighbors(vertex)!=1;
              trim.emplace(vertex, trim_vertex);
              auto vrc_it = vertex2rc.find(vertex);
              if (vrc_it!=vertex2rc.end()) {
                  trim.emplace(vrc_it->second, not trim_vertex);
              }
          }
      }
      return trim;
    }();

    std::vector<Contig> contigs;
    for (const RRVertexType &vertex : *this) {
        auto vertex_it = find(vertex);
        const RRVertexProperty &v_prop = node_prop(vertex_it);
        auto[out_begin, out_end] = out_neighbors(vertex);
        for (auto it = out_begin; it!=out_end; ++it) {
            const RREdgeProperty &e_prop = it->second.prop();
            const RREdgeIndexType e_ind = e_prop.Index();
            auto ecan_it = edge_can.find(it->second.prop().Index());
            // 防御性访问：仅在 edge 标记为 canonical 时输出 contig
            if (ecan_it!=edge_can.end() && ecan_it->second) {

                auto edge_seq_it = edge_seqs.find(e_ind);
                if (edge_seq_it==edge_seqs.end()) {
                    continue;
                }
                Sequence edge_str = edge_seq_it->second;
                uint64_t left = 0;
                auto trim_vertex_it = trim.find(vertex);
                bool trim_vertex = trim_vertex_it==trim.end() ? true : trim_vertex_it->second;
                if (trim_vertex) {
                    left += v_prop.size();
                }
                uint64_t right = edge_str.size();
                auto trim_next_it = trim.find(it->first);
                bool trim_next = trim_next_it==trim.end() ? true : trim_next_it->second;
                if (not trim_next) {
                    right -= node_prop(it->first).size();
                }
                if (left >= right) {
                    continue;
                }
                edge_str = edge_str.Subseq(left, right);
                contigs.emplace_back(std::move(edge_str),
                                     itos(it->second.prop().Index()));
            }
        }
    }
    return contigs;
}

std::vector<Contig> MultiplexDBG::ExportContigs(
    const std::experimental::filesystem::path &f,
    const std::unordered_map<RRVertexType, Sequence> &vertex_seqs,
    const std::unordered_map<RREdgeIndexType, Sequence> &edge_seqs,
    const std::unordered_map<RRVertexType, RRVertexType> &vertex2rc,
    const std::unordered_map<RRVertexType, bool> &vertex_can,
    const std::unordered_map<RREdgeIndexType, bool> &edge_can) const {
    std::ofstream os;
    os.open(f);
    std::vector<Contig> edges =
        GetContigs(vertex_seqs, edge_seqs, vertex2rc, vertex_can, edge_can);
    for (const Contig &contig : edges) {
        os << ">" << contig.id << "\n" << contig.seq << "\n";
    }
    os.close();
    return edges;
}

std::vector<Contig> MultiplexDBG::ExportContigsAndGFA(
    const std::experimental::filesystem::path &contigs_fn,
    const std::experimental::filesystem::path &gfa_fn, size_t threads) const {

    const std::unordered_map<RREdgeIndexType, Sequence>
        edge_seqs = GetEdgeSeqs(threads);
    const std::unordered_map<RRVertexType, Sequence> vertex_seqs =
        GetVertexSeqs(edge_seqs);

    const std::unordered_map<RRVertexType, RRVertexType> vertex2rc =
        MapSeqs2RC<RRVertexType>(vertex_seqs);
    const std::unordered_map<RREdgeIndexType, RREdgeIndexType> edge2rc =
        MapSeqs2RC<RREdgeIndexType>(edge_seqs);

    std::unordered_map<RRVertexType, bool> vertex_can =
        AreSeqsCanonical<RRVertexType>(vertex_seqs);
    std::unordered_map<RREdgeIndexType, bool> edge_can =
        AreSeqsCanonical<RREdgeIndexType>(edge_seqs);

    // 统计 contig 导出路径中的所有缺失映射，避免只看到 GFA 侧统计
    struct ContigExportStats {
        size_t missing_vertex_can = 0;
        size_t missing_vertex2rc = 0;
        size_t missing_edge_can = 0;
        size_t missing_edge_seq = 0;
        size_t missing_trim_vertex = 0;
        size_t missing_trim_next = 0;
        size_t emitted_contigs = 0;
        size_t skipped_contigs = 0;
    } contig_stats;
    {
        std::unordered_map<RRVertexType, bool> trim;
        for (const RRVertexType &vertex : *this) {
            auto vcan_it = vertex_can.find(vertex);
            if (vcan_it==vertex_can.end()) {
                contig_stats.missing_vertex_can += 1;
                continue;
            }
            if (!vcan_it->second) {
                continue;
            }
            const bool trim_vertex = count_out_neighbors(vertex)!=1;
            trim.emplace(vertex, trim_vertex);
            auto vrc_it = vertex2rc.find(vertex);
            if (vrc_it==vertex2rc.end()) {
                contig_stats.missing_vertex2rc += 1;
            } else {
                trim.emplace(vrc_it->second, !trim_vertex);
            }
        }
        for (const RRVertexType &vertex : *this) {
            auto vertex_it = find(vertex);
            const RRVertexProperty &v_prop = node_prop(vertex_it);
            auto[out_begin, out_end] = out_neighbors(vertex);
            for (auto it = out_begin; it!=out_end; ++it) {
                const RREdgeIndexType e_ind = it->second.prop().Index();
                auto ecan_it = edge_can.find(e_ind);
                if (ecan_it==edge_can.end()) {
                    contig_stats.missing_edge_can += 1;
                    contig_stats.skipped_contigs += 1;
                    continue;
                }
                if (!ecan_it->second) {
                    continue;
                }
                auto eseq_it = edge_seqs.find(e_ind);
                if (eseq_it==edge_seqs.end()) {
                    contig_stats.missing_edge_seq += 1;
                    contig_stats.skipped_contigs += 1;
                    continue;
                }
                Sequence edge_str = eseq_it->second;
                auto trim_vertex_it = trim.find(vertex);
                if (trim_vertex_it==trim.end()) {
                    contig_stats.missing_trim_vertex += 1;
                }
                bool trim_vertex = trim_vertex_it==trim.end() ? true : trim_vertex_it->second;
                uint64_t left = trim_vertex ? v_prop.size() : 0;
                auto trim_next_it = trim.find(it->first);
                if (trim_next_it==trim.end()) {
                    contig_stats.missing_trim_next += 1;
                }
                bool trim_next = trim_next_it==trim.end() ? true : trim_next_it->second;
                uint64_t right = edge_str.size();
                if (!trim_next) {
                    right -= node_prop(it->first).size();
                }
                if (left >= right) {
                    contig_stats.skipped_contigs += 1;
                    continue;
                }
                contig_stats.emitted_contigs += 1;
            }
        }
    }

    ExportToGFA(gfa_fn, vertex_seqs, edge_seqs, vertex2rc, edge2rc, vertex_can,
                edge_can);
    // 将 contig 侧一致性统计追加到与 GFA 相同的报告文件
    std::ofstream report(gfa_fn.string() + ".consistency_report.txt", std::ios::app);
    report << "contigs_emitted\t" << contig_stats.emitted_contigs << "\n";
    report << "contigs_skipped\t" << contig_stats.skipped_contigs << "\n";
    report << "contigs_missing_vertex_can\t" << contig_stats.missing_vertex_can << "\n";
    report << "contigs_missing_vertex2rc\t" << contig_stats.missing_vertex2rc << "\n";
    report << "contigs_missing_edge_can\t" << contig_stats.missing_edge_can << "\n";
    report << "contigs_missing_edge_seq\t" << contig_stats.missing_edge_seq << "\n";
    report << "contigs_missing_trim_vertex\t" << contig_stats.missing_trim_vertex << "\n";
    report << "contigs_missing_trim_next\t" << contig_stats.missing_trim_next << "\n";
    report.close();
    std::cout << "Contig consistency: skipped=" << contig_stats.skipped_contigs
              << " (missing_edge_can=" << contig_stats.missing_edge_can
              << ", missing_edge_seq=" << contig_stats.missing_edge_seq
              << ", missing_trim_vertex=" << contig_stats.missing_trim_vertex
              << ", missing_trim_next=" << contig_stats.missing_trim_next << ")" << std::endl;
    return ExportContigs(contigs_fn, vertex_seqs, edge_seqs, vertex2rc,
                         vertex_can, edge_can);
}

void MultiplexDBG::ExportActiveTransitions(
    const std::experimental::filesystem::path &path) const {
    rr_paths->ExportActiveTransitions(path);
}
