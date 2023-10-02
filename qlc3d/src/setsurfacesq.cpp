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
                     const vector<idx>& ind_nodes) {
    auto director = qlc3d::Director::fromDegreeAngles(tiltDegrees, twistDegrees, S0);
    auto tensor = qlc3d::TTensor::fromDirector(director);
    for (unsigned int ind_node : ind_nodes) {
        q.setValue(ind_node,0, tensor.t1());
        q.setValue(ind_node,1, tensor.t2());
        q.setValue(ind_node,2, tensor.t3());
        q.setValue(ind_node,3, tensor.t4());
        q.setValue(ind_node,4, tensor.t5());
    }
}

void setHomeotropic(SolutionVector& q,
                    double S0,
                    const vector<idx> &ind_nodes,
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
    vector<idx> ind_nodes;
    for (int i = 0 ; i < alignment.getnSurfaces() ; i++ ){
        ind_nodes.clear();
        // creates index of all nodes of this alignment surface type
        // 08/02/12 geom->e->FindIndexToMaterialNodes((i+1)*MAT_FIXLC1, &ind_nodes);
        //geom->e->listNodesOfMaterial( ind_nodes, (i+1)*MAT_FIXLC1 );
        geom.getTriangles().listFixLCSurfaces(ind_nodes, i+1);

        if (!ind_nodes.empty()) { // if nodes found
            // get type of current surface
            AnchoringType aType = alignment.surface[i]->getAnchoringType();
            // depending on type, do different things...
            if (aType == Strong) {
                double tilt = alignment.surface[i]->getEasyTilt();
                double twist= alignment.surface[i]->getEasyTwist();
                setGlobalAngles(q, S0, tilt, twist, ind_nodes);
            }
            else if (aType == Homeotropic) {
                setHomeotropic(q, S0, ind_nodes, geom);
            }
            else if (aType == Freeze) {
              // nothing needs to be done
            }
        }// end if alignment nodes found
    }// end for loop over alignment surfaces
}
//end void setStrongSurfaces


void setManualNodesAnchoring(SolutionVector &q, double S0, Surface& surf){
    double tilt = surf.getEasyTilt();
    double twist = surf.getEasyTwist();

    // CONVERT SURFACE PARAMS TO VALID NODE INDEX VECTOR
    vector<idx> nodes_idx;
    for (size_t i = 0 ; i < surf.Params.size() ; i++){
        if ( (surf.Params[i] < 0 ) ){
            RUNTIME_ERROR("Negative node index when setting ManualNodesAnchoring.");
        }
        idx nodeIdx = static_cast<idx>(surf.Params[i]);
        nodes_idx.push_back( nodeIdx);
    }
    setGlobalAngles(q, S0, tilt, twist, nodes_idx);
}

void setSurfacesQ(SolutionVector &q, const Alignment &alignment, double S0,  const Geometry &geom) {
  /*! sets q-tensor values for all surfaces. */
  vector<idx> ind_nodes;

  // loop over all surfaces loaded from settings file
  for (int i = 0 ; i < alignment.getnSurfaces() ; i++) {
    const Surface &surf = alignment.getSurface(i);
    if (!surf.getOverrideVolume()) {
      Log::info("not setting initial orientation for FixLC{} because overrideVolume is false", i + 1);
      continue;
    }

    ind_nodes.clear();
    geom.getTriangles().listFixLCSurfaces(ind_nodes, i+1);
    // MANUAL NODES HAVE NO SURFACE TRIANGLES IN MESH AND MUST BE HANDLED SEPARATELY
    if (surf.getAnchoringType() == ManualNodes) {
      Log::info("FIXLC{} is manual nodes anchoring.", i + 1);
      // MANUAL NODES SHOULD NOT DE DEFINED FOR A FIXLC# THAT IS PRESENT IN THE MESH
      if (!ind_nodes.empty()) {
        RUNTIME_ERROR(fmt::format("Can not set ManualNodes for FIXLC{}. Manual nodes should not be defined for a "
                                  "FIXLC number that is present in the mesh. First available surface for this type is FIXLC{}",
                                  i + 1, alignment.getnSurfaces() + 1));
      }
      setManualNodesAnchoring(q, S0, *(alignment.surface[i]));
      continue;
    }
    //
    // creates index of all nodes of this alignment surface type
    if ( !ind_nodes.empty() ) { // if nodes found

      AnchoringType aType = surf.getAnchoringType();
      double strength = surf.getStrength();

      if (aType == Strong ||
          aType == Weak ||
          (aType == Degenerate && strength >= 0)) {
        double tilt = alignment.surface[i]->getEasyTilt();
        double twist= alignment.surface[i]->getEasyTwist();
        setGlobalAngles(q, S0, tilt, twist, ind_nodes);
      }
        // if homeotropic OR degenerate with negative strength
      else if ((aType == Homeotropic) ||
               ((aType == Degenerate) && (strength < 0))) {
        setHomeotropic(q, S0, ind_nodes, geom);
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
  }// end for loop over alignment surfaces
}





