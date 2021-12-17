# Building

To compile qlc3d from source code, you will required at least a c++ compiler (with support for c++17) and cmake installed on your computer. 

Edit the makefile.def file to set paths and compiler flags, then run Qlc3d make all.

## Building on Linux
The default version of GCC in Ubuntu 20.04 can be used. 

TODO:

## Building in Windows
In the past, qlc3d has successfully been compiled using a 64 bit MinGW compiler, but this hasn't been attempted recently.

TODO:

## Building on a Mac
This will probably not currently work as there are some Windows/Linux file system specific code, but it ***should*** be simple to replace these with standard c++17 code (TODO).

# Running Qlc3d
For instructions on running and configuring qlc3d, see [this document](qlc3d/doc/README.md)

# Examples
Some example projects are included int the [examples](examples/README.md) subdirectory of this project. 

