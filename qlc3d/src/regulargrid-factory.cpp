#include <regulargrid-factory.h>
#include <regulargrid.h>
#include <geometry.h>
#include <util/logging.h>

std::unique_ptr<RegularGrid> buildRegularGrid(size_t nx, size_t ny, size_t nz,
                                               const Geometry& geom) {
  if (nx == 0 || ny == 0 || nz == 0) {
    Log::info("Regular grid generation disabled (nx={}, ny={}, nz={})", nx, ny, nz);
    return nullptr;
  }

  Log::info("Building regular grid with nx={}, ny={}, nz={}", nx, ny, nz);
  auto grid = std::make_unique<RegularGrid>();
  grid->createFromTetMesh(static_cast<unsigned>(nx),
                          static_cast<unsigned>(ny),
                          static_cast<unsigned>(nz),
                          geom.getTetrahedra(),
                          geom.getCoordinates(),
                          geom.getBoundingBox());
  return grid;
}

