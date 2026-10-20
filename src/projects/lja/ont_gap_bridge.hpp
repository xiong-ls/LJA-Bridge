#pragma once

#include "common/logging.hpp"
#include "sequences/seqio.hpp"
#include <vector>

std::vector<Contig> OntGapBridgeAssembly(logging::Logger &logger, const std::vector<Contig> &assembly,
                                         const io::Library &ont_lib);
