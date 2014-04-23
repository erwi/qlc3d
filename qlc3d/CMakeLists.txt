CMAKE_MINIMUM_REQUIRED (VERSION 2.6)
PROJECT (qlc3d)
# Set path to SpaMtrix header files
set(SPAMTRIX_INCLUDES_DIR /home/eero/Projects/SpaMtrix/include)
# Set path to SpaMtrix static library file directory
set(SPAMTRIX_BINARY_DIR  /home/eero/Projects/SpaMtrix/build/lib)
add_definitions(-DLinux) ## Must choose either -DLinux OR -DWindows
add_definitions(-DDEBUG) # DISABLES OPENMP PRAGMAS
set(CMAKE_BUILD_TYPE debug)
#
# Should be no need to change settings below
#
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++0x")
file(GLOB SOURCE_FILES ./includes/*.h ./src/*.cpp)
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMALE_EXE_LINKER_FLAGS "${OpenMP_CXX_FLAGS}")
    add_definitions(-DUSES_OPENMP)
endif(OPENMP_FOUND)
link_directories(${SPAMTRIX_BINARY_DIR})
include_directories(./includes ${SPAMTRIX_INCLUDES_DIR})
add_executable(qlc3d ${SOURCE_FILES})
target_link_libraries(qlc3d spamtrix gomp pthread)
