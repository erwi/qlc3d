CMAKE_MINIMUM_REQUIRED (VERSION 3.0)
PROJECT (qlc3d)

add_definitions(-DLinux) ## Must choose either -DLinux OR -DWindows
add_definitions(-DDEBUG) # DISABLES OPENMP PRAGMAS
#
# Should be no need to change_ settings below
#
#set(CMAKE_BUILD_TYPE )
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++0x")
file(GLOB SOURCE_FILES ./includes/*.h ./src/*.cpp)
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${OpenMP_CXX_FLAGS}")
    add_definitions(-DUSES_OPENMP)
endif(OPENMP_FOUND)
#==========================
# qlc3d library
add_library(lib-qlc3d STATIC
        ${SOURCE_FILES}
        src/simulation-container.cpp
        src/eventhandler.cpp
        src/resultoutput.cpp includes/lc-representation.h src/lc-representation.cpp)
target_compile_features(lib-qlc3d PUBLIC cxx_std_17)
target_include_directories(lib-qlc3d
        PUBLIC includes)
target_link_libraries(lib-qlc3d SpaMtrix gomp pthread)

#========================
# qlc3d application
add_executable(qlc3d src/main-app-qlc3d.cpp)
target_compile_features(qlc3d PUBLIC cxx_std_17)
target_include_directories(qlc3d
        PRIVATE includes)
target_link_libraries(qlc3d lib-qlc3d)
