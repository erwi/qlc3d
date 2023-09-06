#include <resultio.h>
#include <fstream>
#include <stdexcept>
#include <lc-representation.h>
#include <io/vtkiofun.h>
#include <geom/coordinates.h>
#include <geom/vec3.h>
#include <assert.h>

void ResultIO::writeCsvUnstructured(const Coordinates &coordinates,
                                    const SolutionVector &v,
                                    const SolutionVector &q,
                                    const std::string &fileName) {
    std::fstream fs(fileName, fs.out);
    if (!fs.is_open()) {
        throw std::runtime_error("could not open output file:" + fileName);
    }

    fs << "x, y, z, v, nx, ny, nz" << std::endl;
    for (unsigned int i = 0; i < v.getnDoF(); i++) {
        // TODO: handle dielectric regions where no LC material exists

        auto n = qlc3d::TTensor {
            q.getValue(i, 0),
            q.getValue(i, 1),
            q.getValue(i, 2),
            q.getValue(i, 3),
            q.getValue(i, 4)
        }.toDirector();
        auto &p = coordinates.getPoint(i);
        fs << p.x() << "," << p.y() << "," << p.z() << "," <<
           v.getValue(i) << "," << n.nx() << "," << n.ny() << "," << n.nz() << std::endl;
    }

    fs.close();
}

void ResultIO::writeVtkUnstructuredAsciiGrid(
        const Coordinates &coordinates,
        size_t numLcPoints,
        const Mesh &tetMesh,
        const SolutionVector &v,
        const SolutionVector &q,
        const std::string &fileName) {

    assert(numLcPoints <= coordinates.size());
    using namespace vtkIOFun;
    UnstructuredGridWriter writer;
    writer.write(fileName, numLcPoints, coordinates, tetMesh, v, q);
}