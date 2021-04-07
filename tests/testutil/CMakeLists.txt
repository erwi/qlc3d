add_library(lib-test-util STATIC test-util.cpp)
target_include_directories(lib-test-util
    PUBLIC
        ./)

target_compile_features(lib-test-util PUBLIC cxx_std_17)