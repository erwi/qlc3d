#ifndef PROJECT_QLC3D_INITS_H
#define PROJECT_QLC3D_INITS_H
#include <vector>
#include <filesystem>

class Geometry;
class Simu;
class Alignment;
class Electrodes;
class MeshRefinement;
class EventList;

void prepareGeometry(Geometry& geom,
                    const std::filesystem::path &meshFileName,
                    Simu& simu,
                    Alignment& alignment,
                    Electrodes& electrodes);

FILE* createOutputEnergyFile(Simu& simu);

/**
 * converts RefinementConfig objects, as defined in settings file, to mesh refinement events that
 * are inserted to the provided event list
 * @param refinement refinement inputs
 * @param eventListOut created events are placed here
 */
void createMeshRefinementEvents(const MeshRefinement &refinement,
                                EventList &eventListOut);
#endif //PROJECT_QLC3D_INITS_H
