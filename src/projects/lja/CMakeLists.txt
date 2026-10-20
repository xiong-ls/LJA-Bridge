project(debruijn)
set(CMAKE_CXX_STANDARD 14)

add_executable(LJA-Bridge lja.cpp subdataset_processing.cpp gap_closing.cpp uncompressed_output.cpp ont_gap_bridge.cpp)
target_link_libraries(LJA-Bridge lja_ec lja_dbg lja_homopolish lja_common lja_sequence m repeat_resolution)
target_compile_features(LJA-Bridge PRIVATE cxx_std_17) # https://cmake.org/cmake/help/latest/manual/cmake-compile-features.7.html#requiring-language-standards

add_executable(jumboDBG dbg.cpp subdataset_processing.cpp)
target_link_libraries(jumboDBG lja_ec lja_dbg lja_common lja_sequence m)

install(TARGETS LJA-Bridge jumboDBG
        RUNTIME DESTINATION bin
        COMPONENT runtime)
install(PROGRAMS "${PYTHON_DIR}/run_rr.py"
        DESTINATION bin
        COMPONENT runtime)
