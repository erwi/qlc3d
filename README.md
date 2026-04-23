# Building qlc3d

To compile qlc3d from source, you need a C++17-capable compiler, CMake, and Git. To run the full test suite, you also need Python 3.x.

The general steps are:

1. Get the code: clone the repository and initialise its submodules.
2. Configure and build: use CMake for your platform, compiler, and build type.

## Building on Linux (Command Line)

**Getting the Code**

First clone the repository, including its submodules.
```
git clone --recursive https://github.com/erwi/qlc3d.git
cd qlc3d
```

Then configure and build the project with CMake.
```
cmake -S . -B build
cmake --build build
```
The executable should now be available at `build/qlc3d/qlc3d` on Linux.

**Running Tests From the Command Line (Optional)**

Run the tests from the repository root with:
```
ctest --test-dir build --output-on-failure
```
The test suite includes native C++ tests and Python integration tests, so it can take a few minutes to finish.

## Creating a Distributable Release Build in Windows
Following the steps below creates a single executable that should run on any Windows computer. It statically links the dependency libraries with MinGW, so there is no need to worry about missing DLL files. The executable is compiled in release mode, so it includes optimisations and multi-threading.

This assumes you have installed CMake, MinGW-w64, Git, and a version of Bash such as Git Bash. It also assumes the qlc3d source is in `c:\qlc3d`. A Release build on Linux should follow the same CMake workflow.


First, check your Mingw version by running `g++ --version` in Bash, and you should see something like below.
```
$ g++ --version
g++.exe (GCC) 11.2.0
Copyright (C) 2021 Free Software Foundation, Inc.
This is free software; see the source for copying conditions.  There is NO
warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
```

### 1. Create a build directory `c:/qlc3d/build-release`
```
$ cd /c/qlc3d && mkdir build-release && cd build-release
```

### 2. Build the executable `qlc3d.exe`
First run CMake to create the Mingw makefiles for a Release build type:
```
$ cmake -DCMAKE_BUILD_TYPE=Release \
	-DCMAKE_EXE_LINKER_FLAGS="-static-libgcc -static-libstdc++ -static" \
	.. -G "MinGW Makefiles" 
```
This should set the code optimisation switch `-O3`. It may be possible to improve performance further with other optimisation flags. 

Then compile and link the code to produce the executable:
```
$ cmake --build .
```

This should result in `qlc3d.exe` being created in the `qlc3d` subdirectory. For a quick sanity check, run it with:
```
$ ./qlc3d/qlc3d.exe
```
You should see output similar to the following:
```
[INFO] qlc3d. Build date=Aug 19 2023, build time=15:47:51, git commit SHA=6c2cfce, build type=RELEASE.
[INFO] Current directory: C:\qlc3d\build-release
Error in file :./meshes/test.txt
Settings file does not exist: ./meshes/test.txt
[ERROR] An error has occurred
```
It is safe to ignore the errors; they are printed because no settings file was provided.

## Building on a Mac
This is not currently guaranteed to work because some code is still specific to Windows/Linux file systems, but it should be straightforward to replace those parts with standard C++17 code.

# Running Qlc3d
For instructions on running and configuring qlc3d, see [this document](qlc3d/doc/README.md)

# Meshing and Geometry
For instructions on creating meshes and geometry, see [this document](qlc3d/doc/mesh.md)

# Examples
Some example projects are included in the [examples](examples/README.md) subdirectory of this project.

