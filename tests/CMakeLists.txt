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
        cpp/io/meshio-tests.cpp
        cpp/math-tests.cpp
        cpp/io/reader-test.cpp
        cpp/io/vtk-unstructured-grid-writer-tests.cpp
        cpp/io/gmsh-read-test.cpp
        cpp/io/meshreader-test.cpp cpp/init-events-tests.cpp cpp/meshrefinement-test.cpp)

target_compile_features(cpp-test PUBLIC cxx_std_17)

target_include_directories(cpp-test
        PRIVATE ../extern/catch)
target_link_libraries(cpp-test lib-qlc3d lib-test-util)