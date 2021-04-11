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

# python end-to-end test
# NOTE: The working directory is source dir so that the python script run the most up to date version without
# having to run CMake to copy the script to the build directory. Is this a good idea?
add_test(NAME temp-test
        COMMAND ${PYTHON_EXECUTABLE} python/e2e_single-newton-iteration.py ${QLC3D_EXE}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests)

# native c++ tests using catch library
add_test(NAME cpp-test
        COMMAND cpp-test
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/tests)
