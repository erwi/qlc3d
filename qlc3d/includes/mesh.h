#ifndef MESH_H
#define MESH_H

#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <list>
#include <set>
#include <algorithm>
#include <string.h>
#include <math.h>



using std::vector;
using std::set;
using std::list;
class Mesh 
{
private:

    int Dimension;	// number of dimensions of mesh - 2 for tris and 3 for tets
    int nElements;	//total number of elements
    int nNodes;		// number of nodes per element
    int *Elem;				// for tets and tris
    int *Mat;				// for tets and tris
    double *Determinant;	// for tets and tris
    double TotalSize;		// Total volume/area of the mesh

    int *ConnectedVolume; 	// for tris only - index to LC volume element
    double *SurfaceNormal;	// for tris only
    unsigned int MaxNodeNumber;



public:
	
    inline int getnElements()const{return nElements;}   // number of elements
    inline int getnNodes()const {return nNodes;}  	// nuber of nodes per element
    inline int getDimension()const {return Dimension;}	// number of dimesnions of mesh ( 3 / 2 for tets / tris)
    int getConnectedVolume(const int& e); 	// returns index to connected volume element, or -1 if not connected to LC1
    inline int getNode(const int &e, const int &n)const{	// returns node n of element e
	#ifdef DEBUG
	if ( (e<0) || (e>=nElements) || (n < 0) || (n >= nNodes))
	{
	    printf("error - Mesh::getNode(int,int) - index to node out of bounds, bye! ");
	    printf("requested elem %i, node %i\n", e, n);
	    exit(1);
	}
	#endif
	return Elem[e*nNodes + n];
    }
	
    int getMaterialNumber(int e);	// returns material number of element e

    int getFixLCNumber(int e);		// gets alignment layer number, i.e. FixLC 1, 2, 3...
    int getDielectricNumber(int e);     // gets dielectric materials number, i.e. Dielectric 1, 2, 3 ....
    inline unsigned int getMaxNodeNumber() { return MaxNodeNumber;}
    int* getPtrToElement(int e);    		// returns pointer to first node in element e.
    int* getPtrToMaterialNumber(int e); 	// returns pointer to material number of element e
    int* getPtrToConnectedVolume(int e); 	// returns pointer tp connected volume[e]
    double* getPtrToDeterminant(int e);
    double* getPtrToSurfaceNormal(int e);
	
    double Calculate4x4Determinant(double* M);
    double getDeterminant(const int& i) const; // returns value of determinant of element i
    inline double getTotalSize()const{return TotalSize;}
    void AllocateMemory();
    void setConnectedVolume(Mesh* vol);  	// sets indexes to connected LC volume elements
    void setDeterminant(int i, double det); 		// sets determinant i to value det
    void setAllNodes(int *nodes); 			// copies node data to array Elem
    void setAllMaterials(int *mat);			// copies material numbers to array Mat
    void setSurfaceNormal(int i, double norm[3]); 	// sets normal of elementt i to norm
    void setDimension(int i);						// set mesh dimension
    void setnElements(int nnelem);					// set numbero of elements
    void setnNodes(int nnodes);						// set number of nodes / element
    inline void setMaxNodeNumber( const unsigned int& mn) {MaxNodeNumber = mn;}
    void removeElements(std::set <unsigned int>& index);			// removes elements in index. index must be sorted in ascending order
    void ClearMesh();								// clears all data in mesh object
    void addElements(int* new_Elements, int* new_Materials , int num_new);	// adds num_new new elements. these must have same # of nodes/element as existing ones
    void addElements( vector <unsigned int> & m_new, vector <int>& mat_new ); //adds new elements and materials, assuming element types match existing elements (nodes/per element)

    // creates list of all elements of material numebr mat
    void listElementsOfMaterial(std::vector <unsigned int>& elems, const int& mat);


    bool ContainsAllNodes(int elem, int n, int* nodes); // checks if element elem, contains all n nodes in array nodes
    void ContainsNodes(list <int>* elems , list <int>* points ); // adds element number to elemes that contain any node in list points
    bool ContainsCoordinate(const unsigned int& elem, const double* p, const double* coord); // checks whether coordinate is within element elem
    void CompleteNodesSet(const unsigned int& elem, std::vector<unsigned int>& nodes) const; // completes nodes vector with those from element, if nodes is empty returns all elements
    void PrintElements();		// prints all elements
    void PrintElement(int e);	// prints element e
    //void PrintNormals();

    void CalculateDeterminants3D(double *p);
    void FindIndexToMaterialNodes(int mat, vector<int> *index);
    void CalculateSurfaceNormals(double *p, Mesh* tets = NULL);
    void CopySurfaceNormal(int i, double* norm); // copies value of surface normal of element i to norm (which must be array of size3)
    void CopyMesh(Mesh* rhs);	// makes this a copy of Mesh* rhs - why does operator= overloading not work???
	void ScaleDeterminants( const double& s); // scales all determinants by s, e.g. to go to microns
	void CalcLocCoords(const unsigned int& elem, double* p, double* coord, double* loc); // calculates 4 local coordinates of coordinate cord in element elem
	void CalcElemBary(const unsigned int& elem, double* p, double* bary); // calculates barycentre coords of element elem
	double CalcBaryDistSqr( double *p,const unsigned int& elem, double* coord);
	bool isOnXPlane(int e, double X, double* p);	// checks if all nodes of triangle e are on z= Z plane
    bool isOnYPlane(int e, double Y, double* p);	// checks if all nodes of triangle e are on z= Z plane
    bool isOnZPlane(int e, double Z, double* p);	// checks if all nodes of triangle e are on z= Z plane
    //bool isOnBackSurface(int e, Geometry* geom);	//
    bool isNeighbours(const unsigned int& el1, const unsigned int& el2); // checks whether el1 and el2 are neighbours
    void gen_p_to_elem(vector<set <unsigned int> >& p_to_elem); // generates index from points to mesh elements
    void gen_neighbour_list( vector < vector <unsigned int> >& neigh);	// neighbour list of elements. ordered according to local node order.
									// index number of non-existent neighbours elements equals total number of elements
									// i.e. 1 too large to use as an index
    Mesh(int n,int nNodes);
    Mesh();
    ~Mesh();
};

#endif



