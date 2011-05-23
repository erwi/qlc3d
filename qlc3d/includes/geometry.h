#ifndef GEOMETRY_H
#define GEOMETRY_H

# include <mesh.h>
# include <material_numbers.h>
# include <float.h> // MAXIMUM DOUBLE VALUES NEEDED IN COORDINATE COMPARISON
# include <vector>
# include <list>
# include <iostream>
# include <algorithm>
#define EPS 1e-7

using namespace std;
class Geometry
{
	private:
		int np;						// number of nodes
		int npLC;					// number of LC nodes
		double* p;					// nodal coordinates x,y,z
		double* NodeNormals;
                //Oct_Box* oct;				// octree mesh index for fast search
                //Oct_Box* this_oct;
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
	
	
	public:
            Mesh* t;						// volume mesh
            Mesh* e;						// surface mesh
		
            Geometry();
            ~Geometry();

            void setCoordinates(double * coords, int np); 	// copies coords to p, sizeof(coords) is 3 * np
            void addCoordinates(double * coords, int np); 	// adds new coordinates to existing ones by extending p
            void addCoordinates( vector<double>& coords);   // adds new coordinates to end of existing ones
            void setNodeNormals();							// calculates surface node normals
            void setnp(int n);
            void setnpLC(int n);							// set number of LC nodes
            void ReorderDielectricNodes();					// reorder nodes so that dielectric material nodes are last
            void MakePeriEquNodes();						// generates periodic equivalent nodes data structure
            void ClearGeometry();							// clears all data for geometry
            void CreateOctree();
            bool getleft_right_is_periodic();
            bool getfront_back_is_periodic();
            bool gettop_bottom_is_periodic();

            //bool IsCoordinateInElement(double x, double y, double z, int i); // is index to element
            void genIndToTetsByCoords(vector <unsigned int>& ind, double* coord, const unsigned int& nc);
            void brute_force_search(unsigned int &ind, double* coord);
            bool getContainingTet( vector<set< unsigned int> >& p_to_t, double* crd, unsigned int& t0);

            void setTo(Geometry* geom);						// makes this = geom
            void checkForPeriodicGeometry();	// detects type of periodicity of the strucuture
		
            int getnp();
            int getnpLC();
            double* getPtrTop();
            double getpX(int i);	// return node coordinates at node i
            double getpY(int i);
            double getpZ(int i);
            double getXmin();
            double getXmax();
            double getYmin();
            double getYmax();
            double getZmin();
            double getZmax();
		
            double getAbsXDist(int i , double x);	// gets absolute distance between x-coord of node i and x
            double getAbsYDist(int i , double y);	//
            double getAbsZDist(int i , double z);	//
            double getAbsDistSqr(const unsigned int& i , double* coord);
            double* getPtrToNodeNormals();

            size_t getTotalSize(); // returns memory consumption

            void getTetBaryCentre(double* x, const unsigned int& it ); // calculates barycentre x,y,z components of tetrahedron it, x must be array of length 3
            void isValidNodeIndex(const unsigned int& i) const;
		


	// NodeNormal methods	
            double getNodeNormalsX(int i);
            double getNodeNormalsY(int i);
            double getNodeNormalsZ(int i);
            void PrintNodeNormals();
            void PrintNodes();
            void PrintNode(int i);


            bool checkForOverlapingNodes(); // Debug function that chaecks makes sure not nodes are overlapping. Returns TRUE if some are, false if everyting is OK
            void countNodeReferences(vector <int>& refc, Mesh& mesh); // counts the number of times each node is used in mesh. DEBUG

};


void prepareGeometry();




#endif

