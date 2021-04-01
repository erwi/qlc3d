file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/resources/settings.txt
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_executable(temp_test temp_test.cpp)
target_link_libraries(temp_test lib-qlc3d)