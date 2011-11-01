#ifndef FILESYSFUN_H
#define FILESYSFUN_H


/*
    OS independent functions for manipulating, changing etc. of directories are defined here

    Must have either 'Linux' or 'Windows' defined

*/
#include <stdlib.h>
#include <string>
#include <string.h>

#ifdef Linux
    #include <unistd.h>
    #include <sys/stat.h> // for making directories
#endif



bool setCurrentDirectory(const std::string& destdir);

std::string getCurrentDirectory();


bool dirExists(const std::string& dir);

bool createDirectory(const std::string& newdir);




#endif // FILESYSFUN_H
