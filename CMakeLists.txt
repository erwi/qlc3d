cmake_minimum_required(VERSION 3.10.0)
project(project-qlc3d)

add_definitions(-DLinux) ## Must choose either -DLinux OR -DWindows

add_subdirectory(extern/SpaMtrix)
add_subdirectory(qlc3d)

#include(CTest)
include(FindPythonInterp)
enable_testing()
add_subdirectory(tests)
set(QLC3D_EXE $<TARGET_FILE:qlc3d>)
add_test(NAME temp_test
        COMMAND ${PYTHON_EXECUTABLE} python/e2e.py ${QLC3D_EXE}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests
       )
# COMMAND temp_test tests/resources/settings.txt