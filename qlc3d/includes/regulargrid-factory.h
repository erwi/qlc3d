#pragma once
#include <memory>
#include <cstddef>

class Geometry;
class RegularGrid;

/**
 * Build a RegularGrid lookup table from a Geometry object.
 *
 * This is the preferred way to create a RegularGrid outside of the Geometry
 * class hierarchy. Ownership of the returned grid belongs to the caller.
 *
 * @param nx   Number of regular grid points in the x direction (0 disables the grid).
 * @param ny   Number of regular grid points in the y direction (0 disables the grid).
 * @param nz   Number of regular grid points in the z direction (0 disables the grid).
 * @param geom Source geometry (must have coordinates and tetrahedra set).
 * @return Populated RegularGrid, or nullptr if nx, ny, or nz is zero.
 */
std::unique_ptr<RegularGrid> buildRegularGrid(size_t nx, size_t ny, size_t nz,
                                               const Geometry& geom);

