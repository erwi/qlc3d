#ifndef GEOMETRY_H
#define GEOMETRY_H

# include <mesh.h>
# include <material_numbers.h>
# include <float.h> // MAXIMUM DOUBLE VALUES NEEDED IN COORDINATE COMPARISON
# include <vector>
# include <list>
# include <iostream>
# include <algorithm>
# include <limits>
# include <alignment.h>
# include <regulargrid.h>
# include <globals.h>
# define EPS 1e-7

class RegularGrid; // DECLARE HERE TO AVOID CIRCULAR #inlcude

using namespace std;
class Geometry
{
private:
    unsigned int np;						// number of nodes
    unsigned int npLC;					// number of LC nodes
    double* p;					// nodal coordinates x,y,z
    double* NodeNormals;

    double Xmin;
    double Xmax;
    double Ymin;
    double Ymax;
    double Zmin;
    double Zmax;
		
    bool left_right_is_periodic;
    bool front_back_is_periodic;
    bool top_bottom_is_periodic;
	
    vector < list <int> > peri_equ_nodes;

    size_t numWeakSurf;
    size_t* indWeakSurf;     // index to all weak surface triangles
    vector<size_t> periNodes_;

    void setEdgePeriNodes(  list <size_t>& edge0,
                             list <size_t>& edge1,
                            // list <size_t>& edge2,
                            // list <size_t>& edge3,
                             const int& dim);   // edge direction 0,1,2 -> x,y,z

    void setFacePeriNodes( list <size_t>& face0,
                           list <size_t>& face1,
                           const int& norm);    // face normal 0,1,2 -> x,y,z

    void updateMaxNodeNumbers(); // Updates MaxNodeNumbers for surface and tet meshes after a node renumbering
public:
    // UNFORTUNATE HACKERY... SPECIAL INDEX VALUE FOR AN UNSIGNED INDEX THAT WAS NOT FOUND
    static const unsigned int NOT_AN_INDEX;// = std::numeric_limits<unsigned int>::max();

    Mesh* t;						// volume mesh
    Mesh* e;						// surface mesh
    RegularGrid* regularGrid;
    Geometry();
    ~Geometry();

    void setCoordinates(double * coords, const size_t& np); // copies coords to p, sizeof(coords) is 3 * np
    void addCoordinates(double * coords, const size_t& np); // adds new coordinates to existing ones by extending p
    void addCoordinates( vector<double>& coords);  // adds new coordinates to end of existing ones
    void setNodeNormals();		// calculates surface node normals
    void setnp(int n);
    void setnpLC(int n);		// set number of LC nodes
    void ReorderDielectricNodes();	// reorder nodes so that dielectric material nodes are last
    void makePeriEquNodes();	// generates periodic equivalent nodes index
    void ClearGeometry();       // clears all data for geometry
    void CreateOctree();
    bool getleft_right_is_periodic();
    bool getfront_back_is_periodic();
    bool gettop_bottom_is_periodic();


    void genIndToTetsByCoords(vector <unsigned int>& ind,   // Return index
                              double* coord,                // searchable coorinate values
                              const unsigned int& nc,      // length of coord / 3
                              const bool& terminateOnError = true, // terminate app. if coord not found
                              const bool& requireLCElement = false); // only LC elements are considered

    bool brute_force_search(unsigned int &ind,  // found index, return value
                            double* coord,      // pointer to x,y,coords to search
                            const bool& terminateOnError = true, // whether application terminates if coord is not found
                            const bool& requireLCElement = false // only LC elements are considered
                            );
    size_t recursive_neighbour_search(double crd[3],
                                      const vector< set < unsigned int> > & p_to_t,
                                      const size_t& currentTet,
                                      std::set<size_t>& tetHistory,
                                      const bool& requireLCElement = false   // only LC elements are considered
                                      );

    bool getContainingTet( vector<set< unsigned int> >& p_to_t, double crd[3], unsigned int& t0);

    void makeRegularGrid(const size_t& nx, // GENERATES REGULAR GRID LOOKUP INDEXES AND WEIGHTS
                         const size_t& ny,
                         const size_t& nz);

    void setTo(Geometry* geom);						// makes this = geom
    void checkForPeriodicGeometry();	// detects type of periodicity of the strucuture

    size_t getPeriodicEquNode(const size_t& i) const // RETURNS INDEX TO NODE PERIODIC TO i
    {
        if (i<periNodes_.size())
            return periNodes_[i];
        else
            return i;
    }

    unsigned int getnp() const {return np;}
    unsigned int getnpLC()const  {return npLC;}
    double* getPtrTop();
    inline double* getPtrTop(const size_t& i){ if (i<np) return &p[3*i]; return NULL; } // pointer to node i
    double getpX(int i)const;	// return node coordinates at node i
    double getpY(int i)const;
    double getpZ(int i)const;
    double getXmin();
    double getXmax();
    double getYmin();
    double getYmax();
    double getZmin();
    double getZmax();
		
    double getAbsXDist(int i , double x);	// gets absolute distance between x-coord of node i and x
    double getAbsYDist(int i , double y);	//
    double getAbsZDist(int i , double z);	//
    double getAbsDistSqr(const unsigned int i , const double *const coord) const;
    double* getPtrToNodeNormals();

    size_t getTotalSize(); // returns memory consumption

    void getTetBaryCentre(double* x, const unsigned int& it ); // calculates barycentre x,y,z components of tetrahedron it, x must be array of length 3
    void isValidNodeIndex(const unsigned int& i) const;
    void genIndWeakSurfaces(Alignment& alignment);  // generates index to weak surface elements



    // NodeNormal methods
     double getNodeNormalsX(int i);
     double getNodeNormalsY(int i);
     double getNodeNormalsZ(int i);
     void PrintNodeNormals();
     void PrintNodes();
     void PrintNode(int i);
     void PrintPeriodicNodes();

     bool checkForOverlapingNodes(); // Debug function that chaecks makes sure not nodes are overlapping. Returns TRUE if some are, false if everyting is OK
     void countNodeReferences(vector <int>& refc, Mesh& mesh); // counts the number of times each node is used in mesh. DEBUG

};


void prepareGeometry();




#endif

