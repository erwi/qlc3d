#ifndef PROJECT_QLC3D_INITS_H
#define PROJECT_QLC3D_INITS_H
#include <vector>
#include <filesystem>
#include <geom/vec3.h>

class Geometry;
class Coordinates;
class Simu;
class Alignment;
class Electrodes;
class MeshRefinement;
class EventList;
class SolutionVector;
class LC;
class InitialVolumeOrientation;
class Vec3;

void prepareGeometry(Geometry& geom,
                    const std::filesystem::path &meshFileName,
                    Electrodes& electrodes,
                    const Alignment& alignment,
                    const Vec3& stretchVector = {1, 1, 1},
                    unsigned int regularGridCountX = 0,
                    unsigned int regularGridCountY = 0,
                    unsigned int regularGridCountZ = 0);

/**
 * Set up geometry from mesh file, with any electrode potentials set to all 0 and anchoring to to strong with 0 tilt, twist.
 * This is mainly for test purposes.
 */
void prepareGeometryWithDefaultBoundaries(Geometry &geom,
                                          const std::filesystem::path &meshFileName);


FILE* createOutputEnergyFile(Simu& simu);

/**
 * Sets up initial LC solution vector including volume and surface orientations.
 */
void initialiseLcSolutionVector(SolutionVector &q, const Simu &simu, const LC &lc, const InitialVolumeOrientation &boxes, Alignment &alignment, Geometry &geom);

void setVolumeQ(SolutionVector &q, double S0, const InitialVolumeOrientation &boxes, const Coordinates &coordinates);

/** Sets the q-tensor values for the alignment surfaces. Does not mark the nodes as fixed. */
void setSurfacesQ(SolutionVector &q, Alignment &alignment, double S0, const Geometry &geom);

void setStrongSurfacesQ(SolutionVector &q, const Alignment &alignment, double S0, const Geometry &geom);
/**
 * converts RefinementConfig objects, as defined in settings file, to mesh refinement events that
 * are inserted to the provided event list
 * @param refinement refinement inputs
 * @param eventListOut created events are placed here
 */
void createMeshRefinementEvents(const MeshRefinement &refinement,
                                EventList &eventListOut);

void createElectrodeSwitchingEvents(const Electrodes &electrodes,
                                    EventList &eventListOut);
#endif //PROJECT_QLC3D_INITS_H
