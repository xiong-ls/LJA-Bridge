#include "multi_graph.hpp"
#include "subdataset_processing.hpp"
#include "uncompressed_output.hpp"
#include "gap_closing.hpp"
#include "ont_gap_bridge.hpp"
#include "polishing/perfect_alignment.hpp"
#include "error_correction/mult_correction.hpp"
#include "error_correction/mitochondria_rescue.hpp"
#include "error_correction/initial_correction.hpp"
#include "error_correction/manyk_correction.hpp"
#include "repeat_resolution/repeat_resolution.hpp"
#include "error_correction/precorrection.hpp"
#include "sequences/seqio.hpp"
#include "dbg/dbg_construction.hpp"
#include "common/rolling_hash.hpp"
#include "common/dir_utils.hpp"
#include "common/cl_parser.hpp"
#include "common/logging.hpp"
#include <wait.h>
#include <error_correction/dimer_correction.hpp>
#include <polishing/homopolish.hpp>
#include <algorithm>
#include <unordered_map>

using namespace dbg;

static size_t stage_num = 0;
std::vector<Contig> ref;
// ONT 引导冻结阈值：顶点相邻边的最大 ONT 支持小于该值则冻结
static const size_t kOntGuidedMinVertexSupport = 2;
// ONT 引导转移阈值：edge1->edge2 被 ONT 支持次数达到该值才注入 RRPaths
static const size_t kOntGuidedMinTransitionSupport = 3;

//  在每个组装阶段输出图结构.dot、Contig序列.fa及比对信息.ref，用于调试和结果追溯
void PrintPaths(logging::Logger &logger, const std::experimental::filesystem::path &dir, const std::string &stage,
                SparseDBG &dbg, RecordStorage &readStorage, const io::Library &paths_lib, bool small) {
    stage_num += 1;
    std::string stage_name = itos(stage_num) + "_" + stage;
    logger.info() << "Dumping current state. Stage id: " << stage_name << std::endl;
    ensure_dir_existance(dir);
    ensure_dir_existance(dir / "paths");
    printDot(dir / (stage_name + ".dot"), Component(dbg), readStorage.labeler());
    dbg.printFastaOld(dir / (stage_name + ".fasta"));
    if(!small)
        readStorage.printFullAlignments(logger, dir / (stage_name + ".als"));
    std::vector<Contig> paths;
    for(StringContig sc : io::SeqReader(paths_lib)) {
        Contig contig = sc.makeContig();
        if(contig.size() > 100000) {
            paths.emplace_back(contig.seq.Subseq(0, 50000), contig.id + "_start");
            paths.emplace_back(contig.seq.Subseq(contig.size() - 50000), contig.id + "_end");
        } else {
            paths.emplace_back(std::move(contig));
        }
    }
    GraphAlignmentStorage storage(dbg);
    for(Contig &contig : paths) {
        storage.fill(contig);
    }
    for(Contig &contig : paths) {
        ensure_dir_existance(dir / "paths" / contig.getId());
        Component comp = small ? Component::neighbourhood(dbg, contig, dbg.hasher().getK() + 500) :
                Component::longEdgeNeighbourhood(dbg, contig, 20000);
        std::function<std::string(Edge &)> labeler = readStorage.labeler() + storage.labeler();
        printDot(dir / "paths" / contig.getId() / (stage_name + ".dot"), comp, labeler);
    }
    std::ofstream ref_os;
    ref_os.open(dir / (stage_name + ".ref"));
    for(Contig &contig : ref){
        ref_os << contig.getId() << std::endl;
        for(const PerfectAlignment<Contig, Edge> &al : GraphAligner(dbg).carefulAlign(contig)) {
            ref_os << al.seg_from << " " << al.seg_to << std::endl;
        }
    }
    ref_os.close();
}

// 第一阶段纠错：对输入的HiFi reads 进行初始纠错
std::pair<std::experimental::filesystem::path, std::experimental::filesystem::path>
AlternativeCorrection(logging::Logger &logger, const std::experimental::filesystem::path &dir,
            const io::Library &reads_lib, const io::Library &pseudo_reads_lib, const io::Library &paths_lib,
        size_t threads, size_t k, size_t w, double threshold, double reliable_coverage,
        bool close_gaps, bool remove_bad, bool skip, bool debug, bool load) {
    logger.info() << "Performing initial correction with k = " << k << std::endl;
    if (k % 2 == 0) {
        logger.info() << "Adjusted k from " << k << " to " << (k + 1) << " to make it odd" << std::endl;
        k += 1;
    }
    ensure_dir_existance(dir);
    hashing::RollingHash hasher(k, 239);
    std::function<void()> ic_task = [&dir, &logger, &hasher, close_gaps, load, remove_bad, k, w, &reads_lib,
            &pseudo_reads_lib, &paths_lib, threads, threshold, reliable_coverage, debug] {
        io::Library construction_lib = reads_lib + pseudo_reads_lib;
        SparseDBG dbg = load ? DBGPipeline(logger, hasher, w, reads_lib, dir, threads, (dir/"disjointigs.fasta").string(), (dir/"vertices.save").string()) :
                        DBGPipeline(logger, hasher, w, reads_lib, dir, threads);
        dbg.fillAnchors(w, logger, threads);
        size_t extension_size = std::max<size_t>(k * 2, 1000);
        ReadLogger readLogger(threads, dir/"read_log.txt");
        RecordStorage readStorage(dbg, 0, extension_size, threads, readLogger, true, true, false);
        RecordStorage refStorage(dbg, 0, extension_size, threads, readLogger, false, false);
        io::SeqReader reader(reads_lib);
        readStorage.fill(reader.begin(), reader.end(), dbg, w + k - 1, logger, threads);
        coverageStats(logger, dbg);
        if(debug) {
            PrintPaths(logger, dir / "state_dump", "initial", dbg, readStorage, paths_lib, true);
        }
        Precorrect(logger, threads, dbg, readStorage, reliable_coverage);
        RemoveUncovered(logger, threads, dbg, {&readStorage, &refStorage}, extension_size);
        readStorage.trackSuffixes(logger, threads);
//        CorrectDimers(logger, readStorage, k, threads, reliable_coverage);
        correctAT(logger, readStorage, k, threads);
        ManyKCorrect(logger, dbg, readStorage, threshold, reliable_coverage, 800, 4, threads);
        if(debug)
            PrintPaths(logger, dir/ "state_dump", "mk800", dbg, readStorage, paths_lib, true);
        RemoveUncovered(logger, threads, dbg, {&readStorage, &refStorage}, std::max<size_t>(k * 5 / 2, 3000));
        ManyKCorrect(logger, dbg, readStorage, threshold, reliable_coverage, 2000, 4, threads);
        if(debug)
            PrintPaths(logger, dir/ "state_dump", "mk2000", dbg, readStorage, paths_lib, true);
        RemoveUncovered(logger, threads, dbg, {&readStorage, &refStorage}, std::max<size_t>(k * 7 / 2, 5000));
//        CorrectDimers(logger, readStorage, k, threads, reliable_coverage);
        correctAT(logger, readStorage, k, threads);
        correctLowCoveredRegions(logger, dbg, readStorage, refStorage, "/dev/null", threshold, reliable_coverage, k, threads, false);
        ManyKCorrect(logger, dbg, readStorage, threshold, reliable_coverage, 3500, 4, threads);
        RemoveUncovered(logger, threads, dbg, {&readStorage, &refStorage});
        coverageStats(logger, dbg);
        if(debug)
            PrintPaths(logger, dir/ "state_dump", "mk3500", dbg, readStorage, paths_lib, false);
        readStorage.printReadFasta(logger, dir / "corrected.fasta");
        if(debug)
            DrawSplit(Component(dbg), dir / "split");
        dbg.printFastaOld(dir / "graph.fasta");
    };
    if(!skip)
        runInFork(ic_task);
    std::experimental::filesystem::path res;
    res = dir / "corrected.fasta";
    logger.info() << "Initial correction results with k = " << k << " printed to " << res << std::endl;
    return {res, dir / "graph.fasta"};
}

// 无纠错模式: 跳过纠错步骤，直接对原始 reads 构建 de Bruijn 图
std::vector<std::experimental::filesystem::path> NoCorrection(logging::Logger &logger, const std::experimental::filesystem::path &dir,
                const io::Library &reads_lib, const io::Library &pseudo_reads_lib, const io::Library &paths_lib,
                size_t threads, size_t k, size_t w, bool skip, bool debug, bool load) {
    logger.info() << "Performing initial correction with k = " << k << std::endl;
    if (k % 2 == 0) {
        logger.info() << "Adjusted k from " << k << " to " << (k + 1) << " to make it odd" << std::endl;
        k += 1;
    }
    ensure_dir_existance(dir);
    hashing::RollingHash hasher(k, 239);
    std::function<void()> ic_task = [&dir, &logger, &hasher, load, k, w, &reads_lib,
            &pseudo_reads_lib, &paths_lib, threads, debug] {
        io::Library construction_lib = reads_lib + pseudo_reads_lib;
        SparseDBG dbg = load ? DBGPipeline(logger, hasher, w, reads_lib, dir, threads, (dir/"disjointigs.fasta").string(), (dir/"vertices.save").string()) :
                        DBGPipeline(logger, hasher, w, reads_lib, dir, threads);
        dbg.fillAnchors(w, logger, threads);
        size_t extension_size = std::max<size_t>(k * 2, 1000);
        ReadLogger readLogger(threads, dir/"read_log.txt");
        RecordStorage readStorage(dbg, 0, extension_size, threads, readLogger, true, true, false);
        RecordStorage extra_reads(dbg, 0, extension_size, threads, readLogger, false, true, false);
        io::SeqReader reader(reads_lib);
        readStorage.fill(reader.begin(), reader.end(), dbg, w + k - 1, logger, threads);
        coverageStats(logger, dbg);
        if(debug) {
            PrintPaths(logger, dir / "state_dump", "initial", dbg, readStorage, paths_lib, true);
        }
        dbg.printFastaOld(dir / "final_dbg.fasta");
        printDot(dir / "final_dbg.dot", Component(dbg), readStorage.labeler());
        printGFA(dir / "final_dbg.gfa", Component(dbg), true);
        SaveAllReads(dir/"final_dbg.aln", {&readStorage, &extra_reads});
        readStorage.printReadFasta(logger, dir / "corrected_reads.fasta");
    };
    if(!skip)
        runInFork(ic_task);

    return {dir/"corrected_reads.fasta", dir / "final_dbg.fasta", dir / "final_dbg.aln"};
}

// 第二阶段纠错与 Gap 填补
std::vector<std::experimental::filesystem::path> SecondPhase(
    logging::Logger &logger, const std::experimental::filesystem::path &dir,
    const io::Library &reads_lib, const io::Library &pseudo_reads_lib,
    const io::Library &paths_lib, size_t threads, size_t k, size_t w, double threshold, double reliable_coverage,
    size_t unique_threshold, bool diploid, bool skip, bool debug, bool load) {
    logger.info() << "Performing second phase of error correction using k = " << k << std::endl;
    if (k%2==0) {
        logger.info() << "Adjusted k from " << k << " to " << (k + 1)
                      << " to make it odd" << std::endl;
        k += 1;
    }
    ensure_dir_existance(dir);
    hashing::RollingHash hasher(k, 239);
    std::function<void()> ic_task = [&dir, &logger, &hasher, load, k, w,
                                     &reads_lib, &pseudo_reads_lib, &paths_lib,
                                     threads, threshold, reliable_coverage,
                                     debug, unique_threshold, diploid]
                                     {
        io::Library construction_lib = reads_lib + pseudo_reads_lib;
        SparseDBG dbg =
            load ? DBGPipeline(logger, hasher, w, reads_lib, dir, threads,
                               (dir/"disjointigs.fasta").string(),
                               (dir/"vertices.save").string())
                 : DBGPipeline(logger, hasher, w, reads_lib, dir, threads);
        dbg.fillAnchors(w, logger, threads);
        size_t extension_size = 10000000;
        ReadLogger readLogger(threads, dir/"read_log.txt");
        RecordStorage readStorage(dbg, 0, extension_size, threads, readLogger, true, debug);
        RecordStorage refStorage(dbg, 0, extension_size, threads, readLogger, false, false);
        io::SeqReader reader(reads_lib);
        readStorage.fill(reader.begin(), reader.end(), dbg, w + k - 1, logger, threads);
        if(debug) {
            DrawSplit(Component(dbg), dir / "before_figs", readStorage.labeler(), 25000);
            PrintPaths(logger, dir / "state_dump", "initial", dbg, readStorage, paths_lib, false);
        }
        initialCorrect(dbg, logger, dir / "correction.txt", readStorage, refStorage,
                       threshold, 2 * threshold, reliable_coverage, threads, false);
        if(debug)
            PrintPaths(logger, dir/ "state_dump", "low", dbg, readStorage, paths_lib, false);
        GapColserPipeline(logger, threads, dbg, {&readStorage, &refStorage});
        if(debug)
            PrintPaths(logger, dir/ "state_dump", "gap1", dbg, readStorage, paths_lib, false);
        readStorage.invalidateBad(logger, threads, threshold, "after_gap1");
        if(debug)
            PrintPaths(logger, dir/ "state_dump", "bad", dbg, readStorage, paths_lib, false);
        RemoveUncovered(logger, threads, dbg, {&readStorage, &refStorage});
        if(debug)
            PrintPaths(logger, dir/ "state_dump", "uncovered1", dbg, readStorage, paths_lib, false);
        RecordStorage extra_reads = MultCorrect(dbg, logger, dir, readStorage, unique_threshold, threads, diploid, debug);
        MRescue(logger, threads, dbg, readStorage, unique_threshold, 0.05);
        if(debug)
            PrintPaths(logger, dir/ "state_dump", "mult", dbg, readStorage, paths_lib, false);
        RemoveUncovered(logger, threads, dbg, {&readStorage, &extra_reads, &refStorage});
        if(debug)
            PrintPaths(logger, dir/ "state_dump", "uncovered2", dbg, readStorage, paths_lib, false);
        GapColserPipeline(logger, threads, dbg, {&readStorage, &extra_reads, &refStorage});
        if(debug) {
            PrintPaths(logger, dir / "state_dump", "gap2", dbg, readStorage, paths_lib, false);
            DrawSplit(Component(dbg), dir / "split_figs", readStorage.labeler());
        }
        dbg.printFastaOld(dir / "final_dbg.fasta");
        printDot(dir / "final_dbg.dot", Component(dbg), readStorage.labeler());
        printGFA(dir / "final_dbg.gfa", Component(dbg), true);
        SaveAllReads(dir/"final_dbg.aln", {&readStorage, &extra_reads});
        readStorage.printReadFasta(logger, dir / "corrected_reads.fasta");
    };
    if(!skip)
        runInFork(ic_task);
    std::experimental::filesystem::path res;
    res = dir / "corrected_reads.fasta";
    logger.info() << "Second phase results with k = " << k << " printed to "
                  << res << std::endl;
    return {res, dir / "final_dbg.fasta", dir / "final_dbg.aln"};
}

// 重复序列解析（多重de Bruijn 图）
std::vector<std::experimental::filesystem::path> MDBGPhase(
        logging::Logger &logger, size_t threads, size_t k, size_t kmdbg, size_t w, size_t unique_threshold, bool diploid,
        const std::experimental::filesystem::path &dir,
        const std::experimental::filesystem::path &graph_fasta,
        const std::experimental::filesystem::path &read_paths, const io::Library &ont_lib,
        size_t ont_min_len_for_mdbg, bool skip, bool debug) {
    logger.info() << "Performing repeat resolution by transforming de Bruijn graph into Multiplex de Bruijn graph" << std::endl;
    std::function<void()> ic_task = [&logger, threads, debug, k, kmdbg, &graph_fasta, unique_threshold, diploid, &read_paths, &dir,
                                     &ont_lib, ont_min_len_for_mdbg] {
        hashing::RollingHash hasher(k, 239);
        SparseDBG dbg = dbg::LoadDBGFromFasta({graph_fasta}, hasher, logger, threads);
        size_t extension_size = 10000000;
        ReadLogger readLogger(threads, dir/"read_log.txt");
        RecordStorage readStorage(dbg, 0, extension_size, threads, readLogger, true, debug);
        RecordStorage extra_reads(dbg, 0, extension_size, threads, readLogger, false, debug);
        LoadAllReads(read_paths, {&readStorage, &extra_reads}, dbg);

        std::unordered_map<repeat_resolution::RREdgeIndexType, size_t> ont_edge_support;
        std::unordered_map<repeat_resolution::PairEdgeIndexType, size_t, repeat_resolution::PairEdgeIndexHash> ont_transition_support;
        /* ONT-guided MDBG 暂时禁用（结果碎片化，待排查）
        if (!ont_lib.empty()) {
            // 只要提供 ONT reads，就自动启用 ONT 引导统计（容错，不因单条 read 失败中止）
            logger.info() << "ONT-guided MDBG: collecting ONT-supported edges/transitions on current DBG"
                          << " (min ONT read length = " << ont_min_len_for_mdbg << ")" << std::endl;
            std::unordered_map<std::string, repeat_resolution::RREdgeIndexType> edgeid2ind;
            std::unordered_map<repeat_resolution::RREdgeIndexType, repeat_resolution::RREdgeIndexType> edge2rc;
            repeat_resolution::RREdgeIndexType edge_i = 0;
            for (auto it = dbg.edges().begin(); it!=dbg.edges().end(); ++it) {
                edgeid2ind[(*it).getId()] = edge_i;
                ++edge_i;
            }
            // 建立边索引到其反向互补边索引的映射，用于保持 ONT 证据的 RC 对称性
            for (auto it = dbg.edges().begin(); it!=dbg.edges().end(); ++it) {
                const dbg::Edge &edge = *it;
                const repeat_resolution::RREdgeIndexType idx = edgeid2ind.at(edge.getId());
                const auto rc_it = edgeid2ind.find(edge.rc().getId());
                if (rc_it != edgeid2ind.end()) {
                    edge2rc[idx] = rc_it->second;
                }
            }
            // 使用较短锚点进行模糊支持统计，避免 k=5001 这类超长精确锚导致 ONT 全部失配
            const size_t seed_k = 31;
            const size_t edge_seed_step = 5;
            const size_t read_seed_step = 3;
            const size_t max_edges_per_seed = 64;
            const size_t min_seed_hits_per_edge = 2;
            std::unordered_map<hashing::htype, std::vector<repeat_resolution::RREdgeIndexType>, hashing::alt_hasher<hashing::htype>>
                seed2edges;
            hashing::RollingHash seed_hasher(seed_k, 239);
            for (auto it = dbg.edges().begin(); it!=dbg.edges().end(); ++it) {
                const dbg::Edge &edge = *it;
                if (edge.size() < seed_k) {
                    continue;
                }
                const repeat_resolution::RREdgeIndexType idx = edgeid2ind.at(edge.getId());
                hashing::htype h = seed_hasher.hash(edge.seq, 0);
                for (size_t pos = 0;; pos++) {
                    if (pos % edge_seed_step == 0) {
                        seed2edges[h].push_back(idx);
                    }
                    if (!seed_hasher.hasNext(edge.seq, pos)) {
                        break;
                    }
                    h = seed_hasher.next(edge.seq, pos, h);
                }
            }
            io::SeqReader ont_reader(ont_lib);
            size_t ont_total = 0;
            size_t ont_with_hits = 0;
            for (StringContig sc : ont_reader) {
                Contig contig = sc.makeContig();
                if (contig.size() < ont_min_len_for_mdbg || contig.size() < seed_k) {
                    continue;
                }
                ont_total += 1;
                std::unordered_map<repeat_resolution::RREdgeIndexType, std::pair<size_t, size_t>> edge_votes;
                hashing::RollingHash seed_hasher(seed_k, 239);
                hashing::htype h = seed_hasher.hash(contig.seq, 0);
                for (size_t pos = 0;; pos++) {
                    if (pos % read_seed_step == 0) {
                        auto sit = seed2edges.find(h);
                        if (sit != seed2edges.end() && sit->second.size() <= max_edges_per_seed) {
                            for (repeat_resolution::RREdgeIndexType idx : sit->second) {
                                auto &vote = edge_votes[idx];
                                vote.first += 1;
                                if (vote.first == 1) {
                                    vote.second = pos;
                                }
                            }
                        }
                    }
                    if (!seed_hasher.hasNext(contig.seq, pos)) {
                        break;
                    }
                    h = seed_hasher.next(contig.seq, pos, h);
                }
                bool has_hit = false;
                std::vector<std::pair<size_t, repeat_resolution::RREdgeIndexType>> ordered_hits;
                for (const auto &edge_vote : edge_votes) {
                    if (edge_vote.second.first < min_seed_hits_per_edge) {
                        continue;
                    }
                    ont_edge_support[edge_vote.first] += 1;
                    auto rc_it = edge2rc.find(edge_vote.first);
                    if (rc_it != edge2rc.end()) {
                        ont_edge_support[rc_it->second] += 1;
                    }
                    ordered_hits.emplace_back(edge_vote.second.second, edge_vote.first);
                    has_hit = true;
                }
                if (!ordered_hits.empty()) {
                    std::sort(ordered_hits.begin(), ordered_hits.end(),
                              [](const auto &lhs, const auto &rhs) { return lhs.first < rhs.first; });
                    repeat_resolution::RREdgeIndexType prev_edge = ordered_hits.front().second;
                    for (size_t i = 1; i < ordered_hits.size(); i++) {
                        const repeat_resolution::RREdgeIndexType cur_edge = ordered_hits[i].second;
                        if (cur_edge == prev_edge) {
                            continue;
                        }
                        ont_transition_support[{prev_edge, cur_edge}] += 1;
                        auto rc_prev_it = edge2rc.find(prev_edge);
                        auto rc_cur_it = edge2rc.find(cur_edge);
                        if (rc_prev_it != edge2rc.end() && rc_cur_it != edge2rc.end()) {
                            // 反向互补路径方向相反：prev->cur 的 RC 对应为 rc(cur)->rc(prev)
                            ont_transition_support[{rc_cur_it->second, rc_prev_it->second}] += 1;
                        }
                        prev_edge = cur_edge;
                    }
                }
                if (has_hit) {
                    ont_with_hits += 1;
                }
            }
            logger.info() << "ONT-guided MDBG: ONT reads with graph hits = " << ont_with_hits
                          << "/" << ont_total << ", supported edges = " << ont_edge_support.size()
                          << ", supported transitions = " << ont_transition_support.size() << std::endl;
        }
        */
        repeat_resolution::RepeatResolver rr(dbg, &readStorage, {&extra_reads},
                                             k, kmdbg, std::move(ont_edge_support), std::move(ont_transition_support),
                                             kOntGuidedMinVertexSupport, kOntGuidedMinTransitionSupport, dir, unique_threshold,
                                             diploid, debug, logger);
        rr.ResolveRepeats(logger, threads);
    };
    if(!skip)
        runInFork(ic_task);
    return {dir / "assembly.hpc.fasta", dir / "mdbg.hpc.gfa"};
}

// 抛光 + 均聚物解压缩
std::vector<std::experimental::filesystem::path> PolishingPhase( // 抛光并解压同聚物，返回关键输出文件路径
        logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &dir, // 日志器、并行线程数，以及基础工作目录
        const std::experimental::filesystem::path &output_dir, // 输出目录（最终 assembly 等结果写到这里）
        const std::experimental::filesystem::path &gfa_file, // 经过前置步骤生成的 GFA 文件路径
        const std::experimental::filesystem::path &corrected_reads, // 已校正后的读段（用于对齐corrected_reads.fasta）
        const io::Library &reads, const io::Library &ont_reads_for_gap_bridge, size_t dicompress, size_t min_alignment,
        bool skip, bool debug) { // 读段库、解压参数、最小比对长度、是否跳过、是否调试
    logger.info() << "Performing polishing and homopolymer uncompression" << std::endl; // 记录开始进行抛光与同聚物解压的日志
    std::function<void()> ic_task = [&logger, threads, &output_dir, debug, &gfa_file, &corrected_reads, &reads,
                                     &ont_reads_for_gap_bridge, dicompress, min_alignment, &dir] { // 将核心流程封装成一个可运行任务（可交给 runInFork）
        io::SeqReader reader(corrected_reads); // 用校正后的 reads 构建序列读取器
        multigraph::MultiGraph vertex_graph; // 顶点层图结构
        vertex_graph.LoadGFA(gfa_file, true); // 从 GFA 文件加载图（true 表示启用某种构建/解析模式）
        multigraph::MultiGraph edge_graph = vertex_graph.DBG(); // 将顶点图转换为 DBG（边图）形式
        std::vector<Contig> contigs = edge_graph.getEdges(false); // 获取边（contig）集合；false 表示使用默认/非特殊过滤
        auto res = PrintAlignments(logger, threads, contigs, reader.begin(), reader.end(), min_alignment, dir); // 为每个 contig 打印/统计比对结果，并返回用于后续解压的中间信息

        std::vector<Contig> uncompressed = Polish(logger, threads, contigs, res.first, reads, dicompress); // 根据比对结果进行抛光 + 同聚物解压，得到未压缩 contigs
        std::vector<Contig> assembly = printUncompressedResults(logger, threads, edge_graph, uncompressed, output_dir, debug); // 将解压后的结果输出为最终 assembly（可选调试输出）
        
        if (!ont_reads_for_gap_bridge.empty()) {
            // 在 polishing 末端用 ONT 跨端证据执行保守 bridge，提升连续性
            assembly = OntGapBridgeAssembly(logger, assembly, ont_reads_for_gap_bridge);
        }
        
        logger.info() << "Printing final assembly to " << (output_dir / "assembly.fasta") << std::endl; // 记录最终 assembly 写入目标文件的日志
        
        std::ofstream os_cut; // 用于写出 assembly.fasta 的文件流
        os_cut.open(output_dir / "assembly.fasta"); // 打开（或创建）输出的 fasta 文件
        for(Contig &contig : assembly) { // 遍历所有 assembly contig
            if(contig.size() > 1500) // 仅保留长度超过阈值的 contig
                os_cut << ">" << contig.id << "\n" << contig.seq << "\n"; // 以 FASTA 格式写出 contig 记录
        }
        os_cut.close(); // 关闭输出文件流

    }; // ic_task 任务结束
    if(!skip) // 如果没有要求跳过该阶段
        runInFork(ic_task); // 在 fork/子进程环境中执行该任务（隔离计算与资源占用）
    return {output_dir / "assembly.fasta", output_dir / "mdbg.gfa"}; // 返回关键输出文件路径：final assembly fasta 与 mdbg gfa
} 



std::string constructMessage() {
    std::stringstream ss;
    ss << "LJA-Bridge: HiFi genome assembler based on LJA, with optional ONT-guided gap bridging.\n";
    ss << "Usage: LJA-Bridge [options] -o <output-dir> --reads <reads_file> [--reads <reads_file2> ...] [--ont-reads <ont_file> ...]\n\n";
    ss << "Basic options:\n";
    ss << "  -o <file_name> (or --output-dir <file_name>)  Name of output folder. Resulting graph will be stored there.\n";
    ss << "  --reads <file_name>                           Name of file that contains HiFi reads in fasta or fastq format. This option can be used any number of times in the same command line. In this case reads from all specified files will be used as an input.\n";
    ss << "  --ont-reads <file_name>                       Optional ONT reads in fasta/fastq format. If provided together with --use-ont-in-polish, ONT reads are used in polishing and in conservative gap bridging after polishing.\n";
    ss << "  --use-ont-in-polish                           Optional switch. If set, ONT reads are used in polishing and enable ONT gap bridging.\n";
    ss << "  --ont-min-len-for-mdbg <int>                  Minimum ONT read length used in ONT-guided MDBG (currently disabled; default: max(3000, 2*K)).\n";
    ss << "  -h (or --help)                                Print this help message.\n";
    ss << "\nAdvanced options:\n";
    ss << "  -t <int> (or --threads <int>)                 Number of threads. The default value is 16.\n";
    ss << "  -k <int>                                      Value of k used for initial error correction.\n";
    ss << "  -K <int>                                      Value of k used for final error correction and initialization of multiDBG.\n";
    ss << "  --diploid                                     Use this option for diploid genomes. By default LJA assumes that the genome is haploid or inbred.\n";
    return ss.str();
}

int main(int argc, char **argv) {
    CLParser parser({"output-dir=",
                     "threads=16",
                     "k-mer-size=501",
                     "window=2000",
                     "K-mer-size=5001",
                     "KmDBG=40001",
                     "Window=500",
                     "cov-threshold=3",
                     "rel-threshold=10",
                     "Cov-threshold=3",
                     "Rel-threshold=7",
                     "unique-threshold=40000",
                     "ont-min-len-for-mdbg=0",
                     "dump",
                     "dimer-compress=32,32,1",
                     "restart-from=none",
                     "load",
                     "noec",
                     "alternative",
                     "diploid",
                     "use-ont-in-polish",
                     "debug",
                     "help"},  //带默认值的单值选项
                    {"reads", "ont-reads", "paths", "ref"}, //多值选项
                    {"o=output-dir", "t=threads", "k=k-mer-size","w=window", "K=K-mer-size","W=Window", "h=help"},
                    constructMessage()); //生成工具帮助信息的函数
    parser.parseCL(argc, argv); // 传入参数解析用户命令
    if (parser.getCheck("help")) {
        std::cout << parser.message() << std::endl;
        return 0;
    }
    if (!parser.check().empty()) {
        std::cout << "Failed to parse command line parameters." << std::endl;
        std::cout << parser.check() << "\n" << std::endl;
        std::cout << parser.message() << std::endl;
        return 1;
    }

    if(parser.getListValue("reads").size() == 0) {
        std::cout << "Please provide at least one file with reads." << std::endl;
        std::cout << parser.message() << std::endl;
        return 1;
    }

    bool debug = parser.getCheck("debug");
    StringContig::homopolymer_compressing = true; // 开启均聚物压缩
    StringContig::SetDimerParameters(parser.getValue("dimer-compress"));  // 初始化压缩规则
    
    // 初始化输出目录与日志系统
    const std::experimental::filesystem::path dir(parser.getValue("output-dir"));
    ensure_dir_existance(dir);
    logging::LoggerStorage ls(dir, "dbg");
    logging::Logger logger;
    logger.addLogFile(ls.newLoggerFile(), debug ? logging::debug : logging::trace);
    for(size_t i = 0; i < argc; i++) {
        logger << argv[i] << " ";
    }
    logger << std::endl;
    logger.info() << "Hello! You are running LJA-Bridge, a HiFi genome assembler based on LJA with optional ONT gap bridging\n";
    logging::logGit(logger, dir / "version.txt"); //将当前代码的 Git 版本信息写入日志和 dir/version.txt 文件
    // 可选模式系列
    bool diploid = parser.getCheck("diploid"); //二倍体模式
    std::string first_stage = parser.getValue("restart-from"); //中间阶段重启流水线 可选 alternative/phase2/rr/polishing）
    bool skip = first_stage != "none";
    bool load = parser.getCheck("load"); // 是否加载之前保存的 de Bruijn 图数据
    bool noec = parser.getCheck("noec"); // 是否跳过纠错步骤
    bool use_ont_in_polish = parser.getCheck("use-ont-in-polish");
    logger.info() << "LJA-Bridge pipeline started" << std::endl;

    size_t threads = std::stoi(parser.getValue("threads"));

    // 自定义的 “数据库” 结构
    io::Library lib = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("reads"));
    io::Library ont_lib = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("ont-reads"));
    io::Library paths = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("paths"));
    io::Library ref_lib = oneline::initialize<std::experimental::filesystem::path>(parser.getListValue("ref")); //存储参考基因组文件路径
    io::Library mandatory_libs = lib + paths + ref_lib;
    if(!io::CheckLibrary(mandatory_libs)) {
        exit(1);
    }
    if(!ont_lib.empty() && !io::CheckLibrary(ont_lib)) {
        exit(1);
    }
    ref = io::SeqReader(ref_lib).readAllContigs(); // 读取参考基因组序列
    

    size_t k = std::stoi(parser.getValue("k-mer-size")); // 1、小k图
    size_t w = std::stoi(parser.getValue("window")); //窗口筛选kmer
    size_t K = std::stoi(parser.getValue("K-mer-size")); // 2、大K图
    size_t W = std::stoi(parser.getValue("Window"));
    size_t KmDBG = std::stoi(parser.getValue("KmDBG")); // 3、重复序列解析（MDBG 阶段）的 k-mer 长度（默认 40001，专门用于区分重复序列）
    size_t unique_threshold = std::stoi(parser.getValue("unique-threshold")); //重复序列解析的 “唯一阈值”（默认 40000
    size_t ont_min_len_for_mdbg = std::stoull(parser.getValue("ont-min-len-for-mdbg"));
    if (ont_min_len_for_mdbg == 0) {
        ont_min_len_for_mdbg = std::max<size_t>(3000, 2 * K);
    }
    size_t effective_kmdbg = KmDBG;

    // (1)
    std::vector<std::experimental::filesystem::path> corrected_final; // 存储最终纠错后的文件路径
    if(noec) { // 分支 1：无纠错模式（--noec）
        corrected_final = NoCorrection(logger, dir / ("k" + itos(K)), lib, {}, paths, threads, K, W,
                                       skip, debug, load);
    } else { // 分支 2：启用纠错模式（默认）
        // 第一段纠错
        double threshold = std::stod(parser.getValue("cov-threshold"));
        double reliable_coverage = std::stod(parser.getValue("rel-threshold"));
        std::pair<std::experimental::filesystem::path, std::experimental::filesystem::path> corrected1;
        if (first_stage == "alternative")
            skip = false;
        corrected1 = AlternativeCorrection(logger, dir / ("k" + itos(k)), lib, {}, paths, threads, k, w,
                                           threshold, reliable_coverage, false, false, skip, debug, load);
        if (first_stage == "alternative" || first_stage == "none")
            load = false; //（后续阶段不加载旧图）

        // 第二段纠错
        double Threshold = std::stod(parser.getValue("Cov-threshold"));
        double Reliable_coverage = std::stod(parser.getValue("Rel-threshold"));

        if (first_stage == "phase2")
            skip = false;
        corrected_final = SecondPhase(logger, dir / ("k" + itos(K)), {corrected1.first}, {corrected1.second}, paths,
                            threads, K, W, Threshold, Reliable_coverage, unique_threshold, diploid, skip, debug, load);
        if (first_stage == "phase2")
            load = false;
    }

    // (2) MDBG
    if(first_stage == "rr")
        skip = false;
    /* ONT-guided MDBG 暂时禁用
    if (!ont_lib.empty()) {
        logger.info() << "ONT reads are provided, ONT-guided MDBG is enabled by default. ONT files: " << ont_lib.size()
                      << ", min ONT read length for MDBG: " << ont_min_len_for_mdbg
                      << ", min vertex support: " << kOntGuidedMinVertexSupport
                      << ", min transition support: " << kOntGuidedMinTransitionSupport << std::endl;
    }
    */
    std::vector<std::experimental::filesystem::path> resolved =
            MDBGPhase(logger, threads, K, effective_kmdbg, W, unique_threshold, diploid, dir / "mdbg", corrected_final[1],
                      corrected_final[2], ont_lib, ont_min_len_for_mdbg, skip, debug);
    if(first_stage == "rr")
        load = false;



    // (3) Polish
    io::Library polishing_reads = lib;
    if(use_ont_in_polish && !ont_lib.empty()) {
        polishing_reads = lib + ont_lib;
        logger.info() << "Using ONT reads in polishing stage: " << ont_lib.size() << " file(s)" << std::endl;
    } else if(use_ont_in_polish && ont_lib.empty()) {
        logger.info() << "--use-ont-in-polish is set but no --ont-reads are provided. Polishing with HiFi reads only." << std::endl;
    } else if(!use_ont_in_polish && !ont_lib.empty()) {
        logger.info() << "ONT reads are provided but --use-ont-in-polish is not set. Polishing with HiFi reads only." << std::endl;
    }
    if(first_stage == "polishing")
        skip = false;
    std::vector<std::experimental::filesystem::path> uncompressed_results =
            PolishingPhase(logger, threads, dir/ "uncompressing", dir, resolved[1],
                           corrected_final[0],
                           polishing_reads, (use_ont_in_polish ? ont_lib : io::Library{}),
                           StringContig::max_dimer_size / 2, K, skip, debug);
    if(first_stage == "polishing")
        load = false;
    logger.info() << "Final homopolymer compressed and corrected reads can be found here: " << corrected_final[0] << std::endl;
    logger.info() << "Final graph with homopolymer compressed edges can be found here: " << resolved[1] << std::endl;
    logger.info() << "Final graph can be found here: " << uncompressed_results[1] << std::endl;
    logger.info() << "Final assembly can be found here: " << uncompressed_results[0] << std::endl;
    logger.info() << "LJA-Bridge pipeline finished" << std::endl;
    return 0;
}
  