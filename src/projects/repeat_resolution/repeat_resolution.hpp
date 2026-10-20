//
// Created by Andrey Bzikadze on 10/13/21.
//

#pragma once

#include "../dbg/graph_alignment_storage.hpp"
#include "../dbg/sparse_dbg.hpp"
#include "../error_correction/multiplicity_estimation.hpp"
#include "mdbg.hpp"
#include "mdbg_inc.hpp"
#include "paths.hpp"
#include <algorithm>
#include <graphlite/serialize.hpp>
#include <unordered_set>
#include <unordered_map>
#include <vector>

namespace repeat_resolution {

class RepeatResolver {
    dbg::SparseDBG &dbg;
    RecordStorage *reads_storage;
    std::vector<RecordStorage *> extra_storages{};
    std::uint64_t start_k{1};
    std::uint64_t saturating_k{1};
    // ONT 对每条边的支持计数（read-level）
    std::unordered_map<RREdgeIndexType, size_t> ont_edge_support;
    // ONT 对边转移(edge1->edge2)的支持计数（用于 gap bridge 证据）
    std::unordered_map<PairEdgeIndexType, size_t, PairEdgeIndexHash> ont_transition_support;
    // 顶点冻结阈值
    std::uint64_t ont_min_vertex_support{0};
    // 转移注入阈值
    std::uint64_t ont_min_transition_support{0};
    std::experimental::filesystem::path dir;
    uint64_t unique_threshold{0};
    bool diploid{false};
    bool debug{false};
    UniqueClassificator classificator;

    [[nodiscard]] std::vector<RecordStorage *> get_storages() const {
        // 保持原始流程行为：HiFi/read_paths 对应存储始终参与 RRPaths 构建
        std::vector<RecordStorage *> storages = extra_storages;
        storages.push_back(reads_storage);
        return storages;
    }

 public:
    // 在增 k 之前，把 ONT 支持不足的顶点冻结，避免在弱证据区域过度扩展
    void ApplyOntGuidedFreezing(MultiplexDBG &mdbg, logging::Logger &logger) const {
        if (ont_edge_support.empty() || ont_min_vertex_support == 0) {
            return;
        }
        size_t considered_vertices = 0;
        size_t frozen_vertices = 0;
        auto edge_sup = [this](RREdgeIndexType edge) {
            auto it = ont_edge_support.find(edge);
            return it==ont_edge_support.end() ? size_t(0) : it->second;
        };
        for (const RRVertexType &vertex : mdbg) {
            if (mdbg.node_prop(vertex).IsFrozen()) {
                continue;
            }
            considered_vertices += 1;
            const auto [in_edges, out_edges] = mdbg.GetNeighborEdgesIndexes(vertex);
            size_t vertex_support = 0;
            for (RREdgeIndexType edge : in_edges) {
                vertex_support = std::max(vertex_support, edge_sup(edge));
            }
            for (RREdgeIndexType edge : out_edges) {
                vertex_support = std::max(vertex_support, edge_sup(edge));
            }
            if (vertex_support < ont_min_vertex_support) {
                mdbg.FreezeVertex(vertex);
                frozen_vertices += 1;
            }
        }
        logger.info() << "ONT-guided freezing finished: frozen " << frozen_vertices
                      << " of " << considered_vertices
                      << " non-frozen vertices (min support = "
                      << ont_min_vertex_support << ")" << std::endl;
    }

    // 将高支持的 ONT 边转移注入 RRPaths，给 MDBG 提供“可跨 gap 的转移证据”
    void InjectOntTransitions(RRPaths &rr_paths, logging::Logger &logger) const {
        if (ont_transition_support.empty() || ont_min_transition_support == 0) {
            return;
        }
        // 仅在“困难区域”注入：边的任一端顶点出现分叉（in/out>=2）即视为重复/歧义候选区
        std::unordered_set<RREdgeIndexType> difficult_edges;
        {
            RREdgeIndexType edge_i = 0;
            for (auto it = dbg.edges().begin(); it!=dbg.edges().end(); ++it) {
                const dbg::Edge &edge = *it;
                const dbg::Vertex *s = edge.start();
                const dbg::Vertex *e = edge.end();
                const bool difficult =
                    (s->inDeg() >= 2 || s->outDeg() >= 2 || e->inDeg() >= 2 || e->outDeg() >= 2);
                if (difficult) {
                    difficult_edges.insert(edge_i);
                }
                edge_i += 1;
            }
        }
        // 按源边聚合转移支持，后续进行 Top-K 筛选
        std::unordered_map<RREdgeIndexType, size_t> src_total_support;
        std::unordered_map<RREdgeIndexType, std::vector<std::pair<RREdgeIndexType, size_t>>> src2candidates;
        size_t difficult_candidates = 0;
        for (const auto &pair_and_support : ont_transition_support) {
            const PairEdgeIndexType &edge_pair = pair_and_support.first;
            const size_t support = pair_and_support.second;
            if (support < ont_min_transition_support) {
                continue;
            }
            // 只保留困难边相关的转移，避免在易区引入不必要结构扰动
            if (difficult_edges.find(edge_pair.first)==difficult_edges.end() &&
                difficult_edges.find(edge_pair.second)==difficult_edges.end()) {
                continue;
            }
            // 若该转移已被原始路径支持，则不重复注入，降低噪声放大
            if (rr_paths.ContainsPair(edge_pair.first, edge_pair.second)) {
                continue;
            }
            src_total_support[edge_pair.first] += support;
            src2candidates[edge_pair.first].emplace_back(edge_pair.second, support);
            difficult_candidates += 1;
        }

        const size_t max_injected_per_source = 3;
        std::vector<RRPath> augmented_paths = rr_paths.GetPaths();
        size_t injected = 0;
        size_t filtered_by_fraction = 0;
        for (auto &src_and_list : src2candidates) {
            const RREdgeIndexType src = src_and_list.first;
            std::vector<std::pair<RREdgeIndexType, size_t>> &cands = src_and_list.second;
            const size_t total_support = src_total_support[src];
            std::sort(cands.begin(), cands.end(),
                      [](const auto &lhs, const auto &rhs) { return lhs.second > rhs.second; });
            size_t kept_for_src = 0;
            for (const auto &dst_and_support : cands) {
                if (kept_for_src >= max_injected_per_source) {
                    break;
                }
                const RREdgeIndexType dst = dst_and_support.first;
                const size_t support = dst_and_support.second;
                // 占比约束：仅保留该源边上占比 >=20% 的强候选，抑制长尾噪声转移
                if (support * 5 < total_support) {
                    filtered_by_fraction += 1;
                    continue;
                }
                RRPath synthetic_path;
                synthetic_path.id = "ont_bridge_" + std::to_string(src) + "_" +
                                    std::to_string(dst) + "_" + std::to_string(support);
                synthetic_path.edge_list.push_back(src);
                synthetic_path.edge_list.push_back(dst);
                augmented_paths.emplace_back(std::move(synthetic_path));
                injected += 1;
                kept_for_src += 1;
            }
        }
        if (injected == 0) {
            logger.info() << "ONT-guided transitions: no transition passed support threshold (min support = "
                          << ont_min_transition_support << ")" << std::endl;
            return;
        }
        rr_paths = PathsBuilder::FromPathVector(std::move(augmented_paths));
        logger.info() << "ONT-guided transitions injected into RRPaths: " << injected
                      << " synthetic transition paths (difficult_candidates=" << difficult_candidates
                      << ", filtered_by_fraction=" << filtered_by_fraction
                      << ", max_per_source=" << max_injected_per_source << ")" << std::endl;
    }

    RepeatResolver(dbg::SparseDBG &dbg,
                   RecordStorage *reads_storage,
                   std::vector<RecordStorage *> extra_storages,
                   uint64_t start_k,
                   uint64_t saturating_k,
                   std::unordered_map<RREdgeIndexType, size_t> ont_edge_support,
                   std::unordered_map<PairEdgeIndexType, size_t, PairEdgeIndexHash> ont_transition_support,
                   uint64_t ont_min_vertex_support,
                   uint64_t ont_min_transition_support,
                   const std::experimental::filesystem::path &dir,
                   uint64_t unique_threshold,
                   bool diploid,
                   bool debug,
                   logging::Logger &logger)
        : dbg{dbg}, reads_storage{std::move(reads_storage)},
          extra_storages{std::move(extra_storages)}, start_k{start_k},
          saturating_k{saturating_k}, ont_edge_support{std::move(ont_edge_support)},
          ont_transition_support{std::move(ont_transition_support)},
          ont_min_vertex_support{ont_min_vertex_support},
          ont_min_transition_support{ont_min_transition_support}, dir{std::move(dir)},
          unique_threshold{unique_threshold}, diploid{diploid}, debug{debug},
          classificator{dbg, *(this->reads_storage), diploid, debug} {
        std::experimental::filesystem::create_directory(this->dir);
        classificator.classify(logger, unique_threshold, dir/"mult_dir");
        // TODO reactivate filtering
//        for (RecordStorage *const storage : get_storages()) {
//            storage->invalidateSubreads(logger, 1);
//        }
    }

    void ResolveRepeats(logging::Logger &logger, size_t threads) {
        logger.info() << "Resolving repeats" << std::endl;
        logger.info() << "Constructing paths" << std::endl;
        RRPaths rr_paths = PathsBuilder::FromDBGStorages(dbg, get_storages());
        InjectOntTransitions(rr_paths, logger);

        logger.info() << "Building graph" << std::endl;
        MultiplexDBG mdbg(dbg, &rr_paths, start_k, classificator);
        if (debug) {
            logger.trace() << "Checking validity of graph" << std::endl;
            mdbg.AssertValidity();
            logger.trace() << "Graph has passed validity check" << std::endl;
        }
        // logger.info() << "Export to dot" << std::endl;
        // mdbg.ExportToDot(dir/"init_graph.dot");
        // logger.info() << "Export to GFA" << std::endl;
        // mdbg.ExportToGFA(dir/"init_graph.gfa");

        ApplyOntGuidedFreezing(mdbg, logger);
        if (debug) {
            mdbg.AssertValidity();
        }

        logger.info() << "Increasing k" << std::endl;
        MultiplexDBGIncreaser k_increaser{start_k, saturating_k, logger, debug};
        k_increaser.IncreaseUntilSaturation(mdbg, true);
        logger.info() << "Finished increasing k" << std::endl;

        logger.info() << "Exporting remaining active transitions" << std::endl;
        mdbg.ExportActiveTransitions(dir/"mdbg_remaining_trans.txt");

        logger.info() << "Export to Dot" << std::endl;
        mdbg.ExportToDot(dir/"mdbg.hpc.dot");
        logger.info() << "Export to GFA and compressed contigs" << std::endl;
        std::vector<Contig> contigs = mdbg.ExportContigsAndGFA(
            dir/"assembly.hpc.fasta", dir/"mdbg.hpc.gfa", threads);
        logger.info() << "Finished repeat resolution" << std::endl;
    }
};

} // End namespace repeat_resolution