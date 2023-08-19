# Building qlc3d


To compile qlc3d from source code, you will required at least a c++ compiler (with support for c++17) and cmake and Git installed on your computer. To optionally also run all tests, you will need python 3.x installed.

The general steps are

1. Get the code: Clone the repository and initialise its submodules.
2. Compile the code: Use CMake to configure the project for your environemnt which may depend on operating system, compiler used, and any optional build options (RELEASE/DEBUG build etc.)

## Building on Linux (Command Line)

**Getting the Code**

First clone the main code and then dependencies.
```
git clone https://github.com/erwi/qlc3d.git
cd qlc3d
git submodule update --init --remote
```

Then, compile the code using CMake. This will apply the default configuration for the project.
```
mkdir build && cd build
cmake ..
make
```
The executable `qlc3d` file should now exist in the subdirectory `build/qlc3d`. You can then manually copy the file to where you want it to be installed (and possibly add it to the PATH environment variable?).

**Running Tests From the Command Line (Optional)**

While in the build directory created as above, all tests can be executed:
```
ctest ..
```
This starts the tests, and you can expect them to take some minutes to complete.

## Creating a Distributable Release Build in Windows
Following the steps below creates a single executable that should be runnable on any Windows computer. It includes all dependency libraries
statically compiled using Mingw, so there is no need to worry about missing dll files. The executable is compiled in release mode, so has various optimisations and multi-threading
turned on.

This assumes you have installed CMake, Mingw-w64, Git and a version of Bash (e.g. Git Bash) on your build machine. 
It is assumed you have the qlc3d code in directory `c:\qlc3d`. It's likely that creating a Release build in Linux is a very
similar process (not currently tested).


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
$ mingw32-make.exe -j10
```

This should result in the `qlc3d.exe` being created in a subdirectory `qlc3d`. For a quick sanity test, you can run it by 
typing `./qlc3d/qlc3d.exe` in bash and should see something like below:
```
[INFO] qlc3d. Build date=Aug 19 2023, build time=15:47:51, git commit SHA=6c2cfce, build type=RELEASE.
[INFO] Current directory: C:\qlc3d\build-release
Error in file :./meshes/test.txt
Settings file does not exist: ./meshes/test.txt
[ERROR] An error has occurred
```
It's safe to ignore the errors, these are printed because no settings file was provided.  

## Building on a Mac
This will probably not currently work as there are some Windows/Linux file system specific code, but it ***should*** be simple to replace these with standard c++17 code (TODO).

# Running Qlc3d
For instructions on running and configuring qlc3d, see [this document](qlc3d/doc/README.md)

# Meshing and Geometry
For instructions on creating meshes and geometry, see [this document](qlc3d/doc/mesh.md)

# Examples
Some example projects are included int the [examples](examples/README.md) subdirectory of this project. 

