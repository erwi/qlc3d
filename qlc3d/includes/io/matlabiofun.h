#ifndef MATLABIOFUN_H
#define MATLABIOFUN_H
//
//  FUNCTIONS FOR MATLAB FILE WRITING (READING TOO?)
//
#include <fstream>
#include <iostream>
#include <globals.h>
namespace MatlabIOFun {
// WRITES A SINGLE DIMENSIONAL ARRAY
bool writeNumberArray(std::ofstream &fid,
                      const char *varName,
                      const double *values,
                      const idx n);

// WRITES MATRIX, WHERE EACH ROW/COLUMN CORRESPONDS TO A COLUMN ON THE REGULAR GRID
bool writeNumberColumns(std::ofstream &fid,
                        const char *varName,
                        const double *values,
                        const idx nx,
                        const idx ny,
                        const idx nz);
}
#endif // MATLABIOFUN_H
