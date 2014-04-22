CMAKE_MINIMUM_REQUIRED (VERSION 2.6)
PROJECT (qlc3d)

set(SPAMTRIX_INCLUDES_DIR /home/eero/Projects/SpaMtrix/include)
set(SPAMTRIX_BINARY_DIR  /home/eero/Projects/SpaMtrix/build/lib)

add_definitions(-DLinux) ## MUST DEFINE EITHER Linux OR Windows

set(CMAKE_BUILD_TYPE release)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++0x")

file(GLOB SOURCE_FILES ./includes/*.h ./src/*.cpp)

link_directories(${SPAMTRIX_BINARY_DIR})
include_directories(./includes ${SPAMTRIX_INCLUDES_DIR})
add_executable(qlc3d ${SOURCE_FILES})
target_link_libraries(qlc3d spamtrix gomp pthread)
