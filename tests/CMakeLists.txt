file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/resources/settings.txt
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# add_executable(temp-test temp-test.cpp)
#target_link_libraries(temp-test lib-qlc3d)

# tests written in c++, using catch
add_executable(cpp-test cpp-test.cpp cpp-test.cpp)
target_include_directories(cpp-test
        PRIVATE ../extern/catch)
target_link_libraries(cpp-test lib-qlc3d)