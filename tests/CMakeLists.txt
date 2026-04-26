# Copy non-cpp files from source directory to the build directory's "tests" directory so that
# it mirrors the source directory's "tests" directory and same relative paths can be used for accessing
# files in tests.
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/resources DESTINATION ${CMAKE_CURRENT_BINARY_DIR}) # copy resource files
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/python DESTINATION ${CMAKE_CURRENT_BINARY_DIR}) # copy python scripts

# cpp test utility classes
add_subdirectory(testutil)

# tests written in c++, using catch
add_executable(cpp-test
        cpp/cpp-test.cpp
        cpp/configuration-tests.cpp
        cpp/math-tests.cpp
        cpp/meshrefinement-test.cpp
        cpp/meshrefinement-config-tests.cpp
        cpp/refinement/midpoint-node-lookup-tests.cpp
        cpp/refinement/tet-splitter-tests.cpp
        cpp/refinement/tet-classifier-tests.cpp
        cpp/refinement/periodic-edge-expander-tests.cpp
        cpp/refinement/q-tensor-interpolator-tests.cpp
        cpp/refinement/refinement-spec-tests.cpp
        cpp/refinement/refinement-event-tests.cpp
        cpp/refinement/element-selector-tests.cpp
        cpp/initial-lc-orientation-tests.cpp
        cpp/init-events-tests.cpp
        cpp/io/meshio-tests.cpp
        cpp/io/reader-test.cpp
        cpp/io/vtk-unstructured-grid-writer-tests.cpp
        cpp/io/gmsh-read-test.cpp
        cpp/io/meshreader-test.cpp
        cpp/io/result-output-test.cpp
        cpp/geom/geom-tests.cpp
        cpp/geom/tet-mesh-search-tests.cpp
        cpp/potential/potential-solver-test.cpp
        cpp/fe/gaussian-integration-tests.cpp
        cpp/lc/lc-solver-tests.cpp
        cpp/lc/lc-solver-tests.cpp
        cpp/adaptive-time-step-tests.cpp
        cpp/alignment-tests.cpp
        cpp/expression-tests.cpp
        cpp/reader-tests.cpp
        cpp/geom/periodicity-tests.cpp
        cpp/dofmap-tests.cpp
        cpp/init-geometry-tests.cpp
        cpp/mesh/mesh-tests.cpp
        cpp/sparsematrix-tests.cpp
        cpp/mesh/element-split-convert-tests.cpp)

target_compile_features(cpp-test PUBLIC cxx_std_17)

target_include_directories(cpp-test
        PRIVATE ../extern/catch)
target_link_libraries(cpp-test lib-qlc3d lib-test-util)