CMAKE_MINIMUM_REQUIRED (VERSION 3.0)
PROJECT (qlc3d)

# set compiler preprocessor macro QLC3D_SHA to current git branch SHA. This will be logged by the compiled executable
execute_process(COMMAND git rev-parse --short HEAD OUTPUT_VARIABLE QLC3D_COMMIT OUTPUT_STRIP_TRAILING_WHITESPACE)
add_compile_definitions(QLC3D_SHA=\"${QLC3D_COMMIT}\")

#if (LOG_DEBUG_HASH) # enables logging of some hash values. Helpful for detecting unintended changes until better test coverage
#    add_definitions(-DLOG_DEBUG_HASH=1)
#endif(LOG_DEBUG_HASH)

set (CMAKE_CXX_FLAGS "-Wall -std=c++0x")
set (CMAKE_CXX_FLAGS_DEBUG "-g")
set (CMAKE_CXX_FLAGS_RELEASE "-Ofast -DNDEBUG")
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
        src/resultio.cpp
        includes/lc-representation.h
        src/lc-representation.cpp
        includes/io/gmsh-read.h
        src/io/gmsh-read.cpp
        includes/io/meshreader.h
        src/io/meshreader.cpp
        src/io/ReadGiDMesh3D.cpp
        includes/io/matlabiofun.h
        src/io/matlabiofun.cpp
        includes/io/vtkiofun.h
        src/io/vtkiofun.cpp
        includes/inits.h
        includes/util/logging.h
        includes/util/exception.h
        src/util/logging.cpp
        includes/util/stringutil.h
        src/util/stringutil.cpp
        includes/io/result-output.h
        src/io/result-output.cpp
        src/io/writeLCD.cpp
        includes/io/lcview-result-output.h
        src/io/lcview-result-output.cpp
        includes/geom/coordinates.h
        src/geom/coordinates.cpp
        includes/util/hash.h
        includes/potential/potential-solver.h
        src/potential/calcpot3d.cpp
        src/potential/potential-solver.cpp)
target_compile_features(lib-qlc3d PUBLIC cxx_std_17)
set_target_properties(lib-qlc3d PROPERTIES CXX_STANDARD 17)
target_include_directories(lib-qlc3d
        PUBLIC includes)
target_link_libraries(lib-qlc3d SpaMtrix gomp pthread fmt)

#========================
# qlc3d application
add_executable(qlc3d src/main-app-qlc3d.cpp)
target_compile_features(qlc3d PUBLIC cxx_std_17)
set_target_properties(qlc3d PROPERTIES CXX_STANDARD 17)
target_include_directories(qlc3d
        PRIVATE includes)
target_link_libraries(qlc3d lib-qlc3d)
