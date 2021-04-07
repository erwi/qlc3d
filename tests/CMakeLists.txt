file(COPY ${CMAKE_CURRENT_SOURCE_DIR}
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_subdirectory(testutil)

# tests written in c++, using catch
add_executable(cpp-test
        cpp-test.cpp
        configuration-tests.cpp meshio-tests.cpp)

target_compile_features(cpp-test PUBLIC cxx_std_17)

target_include_directories(cpp-test
        PRIVATE ../extern/catch)
target_link_libraries(cpp-test lib-qlc3d lib-test-util)