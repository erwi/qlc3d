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
                    const Vec3& stretchVector,
                    unsigned int regularGridCountX = 1,
                    unsigned int regularGridCountY = 1,
                    unsigned int regularGridCountZ = 1);

FILE* createOutputEnergyFile(Simu& simu);

/**
 * Sets up initial LC solution vector including volume and surface orientations.
 */
void initialiseLcSolutionVector(SolutionVector &q, const Simu &simu, const LC &lc, const InitialVolumeOrientation &boxes, const Alignment &alignment, const Geometry &geom);

void setVolumeQ(SolutionVector &q, double S0, const InitialVolumeOrientation &boxes, const Coordinates &coordinates);
void setSurfacesQ(SolutionVector &q, const Alignment &alignment, double S0, const Geometry &geom);
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
