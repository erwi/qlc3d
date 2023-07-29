#ifndef MATLABIOFUN_H
#define MATLABIOFUN_H
//
//  FUNCTIONS FOR MATLAB FILE WRITING (READING TOO?)
//
#include <fstream>
#include <iostream>
#include <vector>
#include <globals.h>
namespace MatlabIOFun {
// WRITES A SINGLE DIMENSIONAL ARRAY
bool writeNumberArray(std::ofstream &fid,
                      const char *varName,
                      const double *values,
                      const idx n);

/**
 * Writes a matrix to a matlab file, where each row/column corresponds to a column on the regular grid.
 * @param fid output file stream
 * @param varName column variable name
 * @param values column data values
 * @param nx number of x grid points
 * @param ny number of y grid points
 * @param nz number of z grid points
 * @return true if successful
 */
bool writeNumberColumns(std::ofstream &fid, const char *varName, const std::vector<double> &values, idx nx, idx ny, idx nz);
}
#endif // MATLABIOFUN_H
