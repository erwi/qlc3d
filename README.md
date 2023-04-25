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

## Building in Windows
In the past, qlc3d has successfully been compiled using a 64 bit MinGW compiler, but this hasn't been attempted recently. The approach should be similar to building on Linux.


## Building on a Mac
This will probably not currently work as there are some Windows/Linux file system specific code, but it ***should*** be simple to replace these with standard c++17 code (TODO).

# Running Qlc3d
For instructions on running and configuring qlc3d, see [this document](qlc3d/doc/README.md)

# Examples
Some example projects are included int the [examples](examples/README.md) subdirectory of this project. 

