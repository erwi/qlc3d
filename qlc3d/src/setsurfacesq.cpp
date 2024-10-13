#include "../includes/qlc3d.h"
#include <geom/vec3.h>
#include <lc-representation.h>
#include <util/logging.h>
#include <util/exception.h>
#include <inits.h>

void setGlobalAngles(SolutionVector &q,
                     double S0,
                     double tiltDegrees,
                     double twistDegrees,
                     const std::set<idx>& ind_nodes) {
    auto director = qlc3d::Director::fromDegreeAngles(tiltDegrees, twistDegrees, S0);
    for (unsigned int i : ind_nodes) {
      q.setValue(i, director);
    }
}

void setHomeotropic(SolutionVector& q,
                    double S0,
                    const std::set<idx> &ind_nodes,
                    const Geometry &geom){
    for (auto i : ind_nodes) {
      Vec3 normal = geom.getNodeNormal(i);

      qlc3d::Director n(normal.x(), normal.y(), normal.z(), S0);
      auto t = qlc3d::TTensor::fromDirector(n);

      q.setValue(i, 0, t.t1());
      q.setValue(i, 1, t.t2());
      q.setValue(i, 2, t.t3());
      q.setValue(i, 3, t.t4());
      q.setValue(i, 4, t.t5());
    }
}

void setStrongSurfacesQ(SolutionVector &q,
                        const Alignment &alignment,
                        double S0,
                        const Geometry &geom) {
    // sets only strong anchoring surfaces
    for (int i = 0 ; i < alignment.getnSurfaces() ; i++ ){
        std::set<idx> indSurfaceNodes = geom.getTriangles().listFixLCSurfaceNodes(i + 1);

        if (!indSurfaceNodes.empty()) { // if nodes found
          // get type of current surface
            AnchoringType aType = alignment.surface[i].getAnchoringType();
            // depending on type, do different things...
            if (aType == Strong) {
                double tilt = alignment.surface[i].getEasyTilt();
                double twist= alignment.surface[i].getEasyTwist();
                setGlobalAngles(q, S0, tilt, twist, indSurfaceNodes);
            }
            else if (aType == Homeotropic) {
                setHomeotropic(q, S0, indSurfaceNodes, geom);
            }
            else if (aType == Freeze) {
              // nothing needs to be done
            }
        }// end if alignment nodes found
    }// end for loop over alignment surfaces
}
//end void setStrongSurfaces


void setManualNodesAnchoring(SolutionVector &q, double S0, const Surface& surf){
    double tilt = surf.getEasyTilt();
    double twist = surf.getEasyTwist();

    // CONVERT SURFACE PARAMS TO VALID NODE INDEX VECTOR
    std::set<idx> nodes_idx;
    for (size_t i = 0 ; i < surf.Params.size() ; i++){
        if ( (surf.Params[i] < 0 ) ){
            RUNTIME_ERROR("Negative node index when setting ManualNodesAnchoring.");
        }
        idx nodeIdx = static_cast<idx>(surf.Params[i]);
        nodes_idx.insert(nodeIdx);
    }
    setGlobalAngles(q, S0, tilt, twist, nodes_idx);
}

/*! sets q-tensor values for all surfaces. */
void setSurfacesQ(SolutionVector &q, const Alignment &alignment, double S0,  const Geometry &geom) {
  Log::info("Setting initial LC configuration for {} surfaces.", alignment.getnSurfaces());

  // loop over all surfaces loaded from settings file
  for (int i = 0 ; i < alignment.getnSurfaces() ; i++) {
    const Surface &surf = alignment.getSurface(i);
    Log::info("Setting surface {}", surf.toString());

    if (!surf.getOverrideVolume()) {
      Log::info("not setting initial orientation for FixLC{} because overrideVolume is false", i + 1);
      continue;
    }

    std::set<idx> indSurfaceNodes = geom.getTriangles().listFixLCSurfaceNodes(i + 1);

    if (surf.getAnchoringType() == ManualNodes) {
      Log::info("FIXLC{} is manual nodes anchoring.", i + 1);
      // MANUAL NODES SHOULD NOT DE DEFINED FOR A FIXLC# THAT IS PRESENT IN THE MESH
      if (!indSurfaceNodes.empty()) {
        RUNTIME_ERROR(fmt::format("Can not set ManualNodes for FIXLC{}. Manual nodes should not be defined for a "
                                  "FIXLC number that is present in the mesh. First available surface for this type is FIXLC{}",
                                  i + 1, alignment.getnSurfaces() + 1));
      }
      setManualNodesAnchoring(q, S0, alignment.surface[i]);
      continue;
    }
    //
    // creates index of all nodes of this alignment surface type

    if (!indSurfaceNodes.empty() ) { // if nodes found
      AnchoringType aType = surf.getAnchoringType();
      double strength = surf.getStrength();

      if (aType == Strong ||
          aType == Weak ||
          (aType == Degenerate && strength >= 0)) {
        double tilt = alignment.surface[i].getEasyTilt();
        double twist= alignment.surface[i].getEasyTwist();
        setGlobalAngles(q, S0, tilt, twist, indSurfaceNodes);
      }
        // if homeotropic OR degenerate with negative strength
      else if (aType == Homeotropic ||
               (aType == Degenerate && strength < 0) ||
               aType == WeakHomeotropic) {
        setHomeotropic(q, S0, indSurfaceNodes, geom);
      }
      else if (aType == Freeze) {
        // nothing needs to be done, just use the current q-tensor values
      }
      else if ( aType == Polymerise) { // FREEZES ALL NODES WHOSE ORDER IS BELOW VALUE DEFINED IN STRENGTH
        // nothing needs to be done, just use the current q-tensor values
      }
      else {
        RUNTIME_ERROR("Unhandled anchoring type " + surf.getAnchoringTypeName() + ".")
      }
    } else {
      RUNTIME_ERROR(fmt::format("FIXLC{} has ben defined in settings file, but no such surface found in "
                                "the mesh.", i + 1));
    }
  }
}





