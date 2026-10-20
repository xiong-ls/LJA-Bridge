#include "ont_gap_bridge.hpp"

#include "common/rolling_hash.hpp"

#include <spoa/spoa.hpp>
#include <algorithm>
#include <cstdint>
#include <numeric>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {
struct GapBridgeKey {
    size_t left = 0;
    size_t right = 0;
    bool operator==(const GapBridgeKey &other) const {
        return left == other.left && right == other.right;
    }
};

struct GapBridgeKeyHash {
    size_t operator()(const GapBridgeKey &k) const noexcept {
        return std::hash<size_t>{}(k.left) ^ (std::hash<size_t>{}(k.right) << 1);
    }
};

size_t OriginalContigIndex(size_t oriented_node) {
    return oriented_node / 2;
}

size_t ReverseComplementNode(size_t oriented_node) {
    return oriented_node ^ size_t(1);
}

bool IsReverseComplementNode(size_t oriented_node) {
    return oriented_node % 2 != 0;
}

Sequence OrientedContigSequence(const std::vector<Contig> &assembly, size_t oriented_node) {
    const Sequence &seq = assembly[OriginalContigIndex(oriented_node)].seq;
    return IsReverseComplementNode(oriented_node) ? !seq : seq;
}

std::string OrientedContigId(const std::vector<Contig> &assembly, size_t oriented_node) {
    const std::string &id = assembly[OriginalContigIndex(oriented_node)].id;
    return IsReverseComplementNode(oriented_node) ? id + "_rc" : id;
}

GapBridgeKey ReverseComplementBridgeKey(const GapBridgeKey &key) {
    return {ReverseComplementNode(key.right), ReverseComplementNode(key.left)};
}

bool BridgeKeyLess(const GapBridgeKey &lhs, const GapBridgeKey &rhs) {
    return lhs.left < rhs.left || (lhs.left == rhs.left && lhs.right < rhs.right);
}

// 计算 a 的后缀与 b 的前缀最大重叠长度（用于拼接时去除重复片段）
size_t LongestSuffixPrefixOverlap(const Sequence &a, const Sequence &b, size_t min_overlap, size_t max_overlap) {
    if (a.size() == 0 || b.size() == 0) {
        return 0;
    }
    const size_t upper = std::min(max_overlap, std::min(a.size(), b.size()));
    if (upper < min_overlap) {
        return 0;
    }
    for (size_t ov = upper;; ov--) {
        bool ok = true;
        const size_t a_start = a.size() - ov;
        for (size_t i = 0; i < ov; i++) {
            if (a[a_start + i] != b[i]) {
                ok = false;
                break;
            }
        }
        if (ok) {
            return ov;
        }
        if (ov == min_overlap) {
            break;
        }
    }
    return 0;
}

// 将 left + bridge + right 合并为一条序列，并尽量去除端点重叠
Sequence MergeWithBridge(const Sequence &left, const Sequence &bridge, const Sequence &right) {
    const size_t ov_left_bridge = LongestSuffixPrefixOverlap(left, bridge, 21, 500);
    Sequence merged = left + bridge.Subseq(ov_left_bridge);
    const size_t ov_merged_right = LongestSuffixPrefixOverlap(merged, right, 21, 500);
    merged = merged + right.Subseq(ov_merged_right);
    return merged;
}

struct BridgeAnchorInfo {
    size_t index = 0;
    size_t left_overlap = 0;
    size_t right_overlap = 0;
};

size_t BridgeAnchorScore(const BridgeAnchorInfo &anchor) {
    return std::min(anchor.left_overlap, anchor.right_overlap);
}

BridgeAnchorInfo EvaluateBridgeAnchor(const Sequence &left, const Sequence &bridge, const Sequence &right, size_t max_overlap,
                                      size_t index = 0) {
    BridgeAnchorInfo info;
    info.index = index;
    info.left_overlap = LongestSuffixPrefixOverlap(left, bridge, 21, max_overlap);
    info.right_overlap = LongestSuffixPrefixOverlap(bridge, right, 21, max_overlap);
    return info;
}

Sequence BuildBridgeConsensus(const std::vector<Sequence> &bridges) {
    if (bridges.empty()) {
        return Sequence();
    }
    auto alignment_engine = spoa::AlignmentEngine::Create(
            spoa::AlignmentType::kNW, 10, -8, -8, -1);
    spoa::Graph graph{};
    for (const Sequence &bridge : bridges) {
        std::string seq = bridge.str();
        if (seq.empty()) {
            continue;
        }
        std::int32_t score = 0;
        auto alignment = alignment_engine->Align(seq, graph, &score);
        graph.AddAlignment(alignment, seq);
    }
    std::vector<uint32_t> coverages;
    std::string consensus = graph.GenerateConsensus(&coverages);
    if (consensus.empty()) {
        return Sequence();
    }
    return Sequence(consensus);
}

struct BridgePickResult {
    Sequence bridge;
    size_t selected_for_consensus = 0;
    size_t best_anchor_score = 0;
    size_t final_anchor_score = 0;
    bool used_consensus = false;
};

struct SelectedBridgeEdge {
    size_t left = 0;
    size_t right = 0;
    Sequence bridge;
    size_t support = 0;
    size_t selected_for_consensus = 0;
    size_t final_anchor_score = 0;
    bool used_consensus = false;
};

BridgePickResult PickBridgeByAnchorAndConsensus(const Sequence &left, const Sequence &right,
                                                const std::vector<Sequence> &bridges) {
    const size_t kBridgeAnchorMinOverlap = 21;
    const size_t kBridgeAnchorMaxOverlap = 500;
    const size_t kBridgeConsensusMaxReads = 12;
    const size_t kBridgeConsensusMinReads = 3;

    BridgePickResult result;
    if (bridges.empty()) {
        return result;
    }
    std::vector<BridgeAnchorInfo> scored;
    scored.reserve(bridges.size());
    for (size_t i = 0; i < bridges.size(); i++) {
        scored.push_back(EvaluateBridgeAnchor(left, bridges[i], right, kBridgeAnchorMaxOverlap, i));
    }
    std::sort(scored.begin(), scored.end(),
              [&bridges](const BridgeAnchorInfo &lhs, const BridgeAnchorInfo &rhs) {
                  const size_t lhs_score = BridgeAnchorScore(lhs);
                  const size_t rhs_score = BridgeAnchorScore(rhs);
                  if (lhs_score != rhs_score) {
                      return lhs_score > rhs_score;
                  }
                  const size_t lhs_sum = lhs.left_overlap + lhs.right_overlap;
                  const size_t rhs_sum = rhs.left_overlap + rhs.right_overlap;
                  if (lhs_sum != rhs_sum) {
                      return lhs_sum > rhs_sum;
                  }
                  return bridges[lhs.index].size() > bridges[rhs.index].size();
              });
    result.bridge = bridges[scored.front().index];
    result.best_anchor_score = BridgeAnchorScore(scored.front());
    result.final_anchor_score = result.best_anchor_score;

    std::vector<Sequence> consensus_inputs;
    consensus_inputs.reserve(std::min(kBridgeConsensusMaxReads, scored.size()));
    for (const BridgeAnchorInfo &anchor : scored) {
        if (anchor.left_overlap < kBridgeAnchorMinOverlap || anchor.right_overlap < kBridgeAnchorMinOverlap) {
            continue;
        }
        consensus_inputs.push_back(bridges[anchor.index]);
        if (consensus_inputs.size() >= kBridgeConsensusMaxReads) {
            break;
        }
    }
    result.selected_for_consensus = consensus_inputs.size();
    if (consensus_inputs.size() < kBridgeConsensusMinReads) {
        return result;
    }

    Sequence consensus_bridge = BuildBridgeConsensus(consensus_inputs);
    if (consensus_bridge.empty()) {
        return result;
    }
    BridgeAnchorInfo consensus_anchor = EvaluateBridgeAnchor(left, consensus_bridge, right, kBridgeAnchorMaxOverlap);
    const size_t consensus_score = BridgeAnchorScore(consensus_anchor);
    if (consensus_anchor.left_overlap < kBridgeAnchorMinOverlap || consensus_anchor.right_overlap < kBridgeAnchorMinOverlap ||
        consensus_score < result.best_anchor_score) {
        return result;
    }
    result.bridge = std::move(consensus_bridge);
    result.final_anchor_score = consensus_score;
    result.used_consensus = true;
    return result;
}
} // namespace

// 在LJA代码中调用OntGapBridgeAssembly(logger, assembly, ont_reads_for_gap_bridge);
// 利用 ONT 跨端锚定信息对 assembly 进行保守 gap bridge 合并
std::vector<Contig> OntGapBridgeAssembly(logging::Logger &logger, const std::vector<Contig> &assembly,
                                         const io::Library &ont_lib) {
    if (assembly.empty() || ont_lib.empty()) {
        return assembly;
    }
    const size_t seed_k = 31;
    const size_t end_window = 4000;
    const size_t edge_seed_step = 7;
    const size_t read_seed_step = 5;
    // 收紧阈值：提高最小seed票数，降低随机短锚命中导致的伪桥接
    const size_t min_seed_votes = 5;
    // 收紧阈值：提高最小桥接支持读数，仅保留更高置信bridge
    const size_t min_bridge_support = 6;
    // 收紧阈值：限制跨距，避免超长跨度导致错误连接
    const size_t max_bridge_span = 20000;
    struct EndHit {
        size_t oriented_node = 0;
        bool is_suffix = false;
    };
    std::unordered_map<hashing::htype, std::vector<EndHit>, hashing::alt_hasher<hashing::htype>> seed2hits;
    hashing::RollingHash seed_hasher(seed_k, 239);
    const size_t oriented_node_count = assembly.size() * 2;
    for (size_t node = 0; node < oriented_node_count; node++) {
        Sequence seq = OrientedContigSequence(assembly, node);
        if (seq.size() < seed_k) {
            continue;
        }
        const size_t pref_len = std::min(end_window, seq.size());
        const size_t suf_len = std::min(end_window, seq.size());
        Sequence pref = seq.Subseq(0, pref_len);
        Sequence suf = seq.Subseq(seq.size() - suf_len);
        hashing::htype hp = seed_hasher.hash(pref, 0);
        for (size_t pos = 0;; pos++) {
            if (pos % edge_seed_step == 0) {
                seed2hits[hp].push_back({node, false});
            }
            if (!seed_hasher.hasNext(pref, pos)) {
                break;
            }
            hp = seed_hasher.next(pref, pos, hp);
        }
        hashing::htype hs = seed_hasher.hash(suf, 0);
        for (size_t pos = 0;; pos++) {
            if (pos % edge_seed_step == 0) {
                seed2hits[hs].push_back({node, true});
            }
            if (!seed_hasher.hasNext(suf, pos)) {
                break;
            }
            hs = seed_hasher.next(suf, pos, hs);
        }
    }
    struct BridgeEvidence {
        size_t support = 0;
        std::vector<Sequence> bridges;
    };
    std::unordered_map<GapBridgeKey, BridgeEvidence, GapBridgeKeyHash> evidence_map;
    size_t reads_with_bridge_signal = 0;
    size_t total_reads = 0;
    for (StringContig sc : io::SeqReader(ont_lib)) {
        // 这里使用原始（未同聚体压缩）ONT 序列，避免与 uncompressed assembly 域不一致导致零命中
        Contig read(Sequence(sc.seq), sc.id);
        if (read.size() < seed_k) {
            continue;
        }
        total_reads += 1;
        struct ReadBridgeCall {
            bool valid = false;
            GapBridgeKey key;
            Sequence bridge;
            size_t score = 0;
            size_t vote_sum = 0;
        };

        auto find_bridge_call = [&](const Sequence &seq) {
            ReadBridgeCall call;
            if (seq.size() < seed_k) {
                return call;
            }
            std::unordered_map<size_t, std::pair<size_t, size_t>> suffix_votes;
            std::unordered_map<size_t, std::pair<size_t, size_t>> prefix_votes;
            hashing::htype h = seed_hasher.hash(seq, 0);
            for (size_t pos = 0;; pos++) {
                if (pos % read_seed_step == 0) {
                    auto it = seed2hits.find(h);
                    if (it != seed2hits.end()) {
                        for (const EndHit &hit : it->second) {
                            auto &bucket = hit.is_suffix ? suffix_votes : prefix_votes;
                            auto &v = bucket[hit.oriented_node];
                            v.first += 1;
                            if (v.first == 1) {
                                v.second = pos;
                            }
                        }
                    }
                }
                if (!seed_hasher.hasNext(seq, pos)) {
                    break;
                }
                h = seed_hasher.next(seq, pos, h);
            }

            size_t best_suffix_idx = size_t(-1);
            size_t best_suffix_votes = 0;
            size_t best_suffix_pos = 0;
            for (const auto &it : suffix_votes) {
                if (it.second.first >= min_seed_votes && it.second.first > best_suffix_votes) {
                    best_suffix_idx = it.first;
                    best_suffix_votes = it.second.first;
                    best_suffix_pos = it.second.second;
                }
            }
            size_t best_prefix_idx = size_t(-1);
            size_t best_prefix_votes = 0;
            size_t best_prefix_pos = 0;
            for (const auto &it : prefix_votes) {
                if (it.second.first >= min_seed_votes && it.second.first > best_prefix_votes) {
                    best_prefix_idx = it.first;
                    best_prefix_votes = it.second.first;
                    best_prefix_pos = it.second.second;
                }
            }
            if (best_suffix_idx == size_t(-1) || best_prefix_idx == size_t(-1) ||
                OriginalContigIndex(best_suffix_idx) == OriginalContigIndex(best_prefix_idx)) {
                return call;
            }
            if (best_prefix_pos <= best_suffix_pos + seed_k) {
                return call;
            }
            const size_t bridge_start = best_suffix_pos;
            const size_t bridge_end = best_prefix_pos + seed_k;
            if (bridge_end - bridge_start > max_bridge_span) {
                return call;
            }
            Sequence bridge_seq = seq.Subseq(bridge_start, bridge_end);
            if (bridge_seq.size() < seed_k) {
                return call;
            }
            GapBridgeKey key{best_suffix_idx, best_prefix_idx};
            const GapBridgeKey rc_key = ReverseComplementBridgeKey(key);
            call.valid = true;
            if (BridgeKeyLess(rc_key, key)) {
                call.key = rc_key;
                call.bridge = !bridge_seq;
            } else {
                call.key = key;
                call.bridge = std::move(bridge_seq);
            }
            call.score = std::min(best_suffix_votes, best_prefix_votes);
            call.vote_sum = best_suffix_votes + best_prefix_votes;
            return call;
        };

        // 分别扫描正向和反向互补；bridge 序列必须从产生命中的同一方向截取。
        ReadBridgeCall forward_call = find_bridge_call(read.seq);
        Sequence rc_read = !read.seq;
        ReadBridgeCall rc_call = find_bridge_call(rc_read);
        const ReadBridgeCall *best_call = nullptr;
        if (forward_call.valid) {
            best_call = &forward_call;
        }
        if (rc_call.valid &&
            (best_call == nullptr || rc_call.score > best_call->score ||
             (rc_call.score == best_call->score && rc_call.vote_sum > best_call->vote_sum) ||
             (rc_call.score == best_call->score && rc_call.vote_sum == best_call->vote_sum &&
              rc_call.bridge.size() < best_call->bridge.size()))) {
            best_call = &rc_call;
        }
        if (best_call == nullptr) {
            continue;
        }
        evidence_map[best_call->key].support += 1;
        evidence_map[best_call->key].bridges.push_back(best_call->bridge);
        reads_with_bridge_signal += 1;
    }
    std::vector<std::pair<GapBridgeKey, BridgeEvidence>> accepted;
    for (auto &it : evidence_map) {
        if (it.second.support >= min_bridge_support) {
            accepted.emplace_back(it.first, std::move(it.second));
        }
    }
    std::vector<SelectedBridgeEdge> candidate_edges;
    candidate_edges.reserve(accepted.size());
    size_t candidate_consensus_pairs = 0;
    size_t candidate_fallback_pairs = 0;
    size_t low_anchor_pairs = 0;
    size_t consensus_input_reads = 0;
    for (const auto &cand : accepted) {
        const size_t left = cand.first.left;
        const size_t right = cand.first.right;
        const auto &bridges = cand.second.bridges;
        Sequence left_seq = OrientedContigSequence(assembly, left);
        Sequence right_seq = OrientedContigSequence(assembly, right);
        BridgePickResult bridge_pick = PickBridgeByAnchorAndConsensus(left_seq, right_seq, bridges);
        if (bridge_pick.bridge.empty()) {
            continue;
        }
        if (bridge_pick.best_anchor_score < 21) {
            low_anchor_pairs += 1;
        }
        if (bridge_pick.used_consensus) {
            candidate_consensus_pairs += 1;
            consensus_input_reads += bridge_pick.selected_for_consensus;
        } else {
            candidate_fallback_pairs += 1;
        }
        candidate_edges.push_back({left, right, std::move(bridge_pick.bridge), cand.second.support,
                                   bridge_pick.selected_for_consensus, bridge_pick.final_anchor_score,
                                   bridge_pick.used_consensus});
    }
    auto bridge_confidence_score = [](const SelectedBridgeEdge &edge) {
        const long long support_score = static_cast<long long>(edge.support) * 500;
        const long long anchor_score = static_cast<long long>(std::min<size_t>(edge.final_anchor_score, 100)) * 5;
        const long long consensus_score = edge.used_consensus ? 150 : 0;
        const long long consensus_read_score = static_cast<long long>(edge.selected_for_consensus) * 10;
        const long long span_penalty = static_cast<long long>(std::min<size_t>(edge.bridge.size() / 100, 200));
        return support_score + anchor_score + consensus_score + consensus_read_score - span_penalty;
    };
    std::sort(candidate_edges.begin(), candidate_edges.end(),
              [&bridge_confidence_score](const SelectedBridgeEdge &lhs, const SelectedBridgeEdge &rhs) {
                  const long long lhs_score = bridge_confidence_score(lhs);
                  const long long rhs_score = bridge_confidence_score(rhs);
                  if (lhs_score != rhs_score) {
                      return lhs_score > rhs_score;
                  }
                  if (lhs.support != rhs.support) {
                      return lhs.support > rhs.support;
                  }
                  if (lhs.final_anchor_score != rhs.final_anchor_score) {
                      return lhs.final_anchor_score > rhs.final_anchor_score;
                  }
                  if (lhs.used_consensus != rhs.used_consensus) {
                      return lhs.used_consensus;
                  }
                  return lhs.bridge.size() < rhs.bridge.size();
              });

    const size_t npos = size_t(-1);
    std::vector<size_t> outgoing(oriented_node_count, npos);
    std::vector<size_t> incoming(oriented_node_count, npos);
    std::vector<size_t> edge_by_left(oriented_node_count, npos);
    std::vector<SelectedBridgeEdge> selected_edges;
    selected_edges.reserve(candidate_edges.size() * 2);
    size_t selected_bridge_pairs = 0;
    size_t multi_edge_conflicts_skipped = 0;
    size_t cycle_edges_skipped = 0;
    size_t orientation_conflicts_skipped = 0;
    auto would_create_cycle = [&](size_t left, size_t right) {
        for (size_t cur = right; cur != npos; cur = outgoing[cur]) {
            if (cur == left) {
                return true;
            }
        }
        return false;
    };
    auto would_reuse_original_contig = [&](size_t left, size_t right) {
        std::unordered_set<size_t> upstream_contigs;
        for (size_t cur = left; cur != npos; cur = incoming[cur]) {
            upstream_contigs.insert(OriginalContigIndex(cur));
        }
        for (size_t cur = right; cur != npos; cur = outgoing[cur]) {
            if (upstream_contigs.count(OriginalContigIndex(cur)) != 0) {
                return true;
            }
        }
        return false;
    };
    for (SelectedBridgeEdge &edge : candidate_edges) {
        const size_t left = edge.left;
        const size_t right = edge.right;
        const size_t rc_left = ReverseComplementNode(right);
        const size_t rc_right = ReverseComplementNode(left);
        if (outgoing[left] != npos || incoming[right] != npos ||
            outgoing[rc_left] != npos || incoming[rc_right] != npos) {
            multi_edge_conflicts_skipped += 1;
            continue;
        }
        if (would_create_cycle(left, right) || would_create_cycle(rc_left, rc_right)) {
            cycle_edges_skipped += 1;
            continue;
        }
        if (would_reuse_original_contig(left, right) ||
            would_reuse_original_contig(rc_left, rc_right)) {
            orientation_conflicts_skipped += 1;
            continue;
        }
        SelectedBridgeEdge rc_edge{rc_left, rc_right, !edge.bridge, edge.support,
                                   edge.selected_for_consensus, edge.final_anchor_score,
                                   edge.used_consensus};
        edge_by_left[left] = selected_edges.size();
        outgoing[left] = right;
        incoming[right] = left;
        selected_edges.push_back(std::move(edge));
        edge_by_left[rc_left] = selected_edges.size();
        outgoing[rc_left] = rc_right;
        incoming[rc_right] = rc_left;
        selected_edges.push_back(std::move(rc_edge));
        selected_bridge_pairs += 1;
    }

    std::vector<Contig> merged;
    std::vector<bool> consumed(assembly.size(), false);
    size_t bridge_paths = 0;
    size_t path_consumed_original_contigs = 0;
    for (size_t start = 0; start < oriented_node_count; start++) {
        if (incoming[start] != npos || outgoing[start] == npos) {
            continue;
        }
        if (consumed[OriginalContigIndex(start)]) {
            continue;
        }
        size_t cur = start;
        Sequence merged_seq = OrientedContigSequence(assembly, cur);
        std::string merged_id = OrientedContigId(assembly, cur);
        consumed[OriginalContigIndex(cur)] = true;
        size_t path_nodes = 1;
        while (outgoing[cur] != npos) {
            const size_t edge_idx = edge_by_left[cur];
            const SelectedBridgeEdge &edge = selected_edges[edge_idx];
            const size_t next = edge.right;
            Sequence next_seq = OrientedContigSequence(assembly, next);
            merged_seq = MergeWithBridge(merged_seq, edge.bridge, next_seq);
            merged_id += "_ONTgap_" + OrientedContigId(assembly, next);
            cur = next;
            consumed[OriginalContigIndex(cur)] = true;
            path_nodes += 1;
        }
        merged.emplace_back(merged_seq, merged_id);
        bridge_paths += 1;
        path_consumed_original_contigs += path_nodes;
    }
    std::vector<Contig> result = merged;
    size_t unused_original_contigs = 0;
    for (size_t i = 0; i < assembly.size(); i++) {
        if (!consumed[i]) {
            result.push_back(assembly[i]);
            unused_original_contigs += 1;
        }
    }
    logger.info() << "ONT gap-bridge summary: total_ont_reads=" << total_reads
                  << ", reads_with_bridge_signal=" << reads_with_bridge_signal
                  << ", oriented_contig_nodes=" << oriented_node_count
                  << ", accepted_bridge_candidates=" << accepted.size()
                  << ", candidate_bridge_edges=" << candidate_edges.size()
                  << ", selected_bridge_pairs=" << selected_bridge_pairs
                  << ", selected_oriented_edges=" << selected_edges.size()
                  << ", bridge_paths=" << bridge_paths
                  << ", path_merged_contigs=" << merged.size()
                  << ", path_consumed_original_contigs=" << path_consumed_original_contigs
                  << ", unused_original_contigs=" << unused_original_contigs
                  << ", multi_edge_conflicts_skipped=" << multi_edge_conflicts_skipped
                  << ", cycle_edges_skipped=" << cycle_edges_skipped
                  << ", orientation_conflicts_skipped=" << orientation_conflicts_skipped
                  << ", consensus_pairs=" << candidate_consensus_pairs
                  << ", fallback_pairs=" << candidate_fallback_pairs
                  << ", low_anchor_pairs=" << low_anchor_pairs
                  << ", consensus_input_reads=" << consensus_input_reads
                  << ", contigs_before=" << assembly.size()
                  << ", contigs_after=" << result.size() << std::endl;
    return result;
}
