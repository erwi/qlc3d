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
#include <globals.h>


using std::vector;
using std::set;
using std::list;
class Mesh 
{
private:

    idx Dimension;	// number of dimensions of mesh - 2 for tris and 3 for tets
    idx nElements;	//total number of elements
    idx nNodes;		// number of nodes per element
    idx *Elem;		// for tets and tris
    idx *Mat;		// for tets and tris
    double *Determinant;// for tets and tris
    double TotalSize;	// Total volume/area of the mesh

    idx *ConnectedVolume; 	// for tris only - index to LC volume element
    double *SurfaceNormal;	// for tris only
    idx MaxNodeNumber;

public:
	
    inline idx getnElements()const{return nElements;}   // number of elements

    inline idx getnNodes()const {return nNodes;}  	// nuber of nodes per element
    inline idx getDimension()const {return Dimension;}	// number of dimesnions of mesh ( 3 / 2 for tets / tris)
    idx getConnectedVolume(const idx e) const; 	// returns index to connected volume element, or -1 if not connected to LC1
    inline idx getNode(const idx e, const idx n)const	// returns node n of element e
    {
	#ifdef DEBUG
    if ( (e>=nElements) || (n >= nNodes))
	{
	    printf("error - Mesh::getNode(int,int) - index to node out of bounds, bye! ");
            printf("requested elem %u, node %u\n", e, n);
	    exit(1);
	}
	#endif
	return Elem[e*nNodes + n];
    }
	
    idx getMaterialNumber(const idx e) const;	// returns material number of element e
    idx getFixLCNumber(const idx e) const;	// gets alignment layer number, i.e. FixLC 1, 2, 3...
    idx getDielectricNumber(const idx e) const; // gets dielectric materials number, i.e. Dielectric 1, 2, 3 ....
    idx getMaxNodeNumber() const { return MaxNodeNumber;}
    idx* getPtrToElement(const idx e) const;    // returns pointer to first node in element e.
    idx* getPtrToMaterialNumber(const idx e) const; // returns pointer to material number of element e
    idx* getPtrToConnectedVolume(const idx e) const;// returns pointer tp connected volume[e]
    double* getPtrToDeterminant(const idx e) const;
    double* getPtrToSurfaceNormal(const idx e) const;
	
    double Calculate4x4Determinant(double* M) const;
    double getDeterminant(const idx i) const; // returns value of determinant of element i
    inline double getTotalSize()const{return TotalSize;}
    void AllocateMemory();
    void setConnectedVolume(Mesh* vol);  	// sets indexes to connected LC volume elements
    void setDeterminant(idx i, double det); 		// sets determinant i to value det
    void setAllNodes(idx *nodes); 			// copies node data to array Elem
    void setAllMaterials(idx *mat);			// copies material numbers to array Mat
    void setSurfaceNormal(idx i, double norm[3]); 	// sets normal of elementt i to norm
    void setDimension(idx i);						// set mesh dimension
    void setnElements(idx nnelem);					// set numbero of elements
    void setnNodes(idx nnodes);						// set number of nodes / element
    inline void setMaxNodeNumber( const idx mn) {MaxNodeNumber = mn;}
    void removeElements(std::set <idx>& index);			// removes elements in index. index must be sorted in ascending order
    void ClearMesh();								// clears all data in mesh object
    //void addElements(idx* new_Elements, idx* new_Materials , idx num_new);	// adds num_new new elements. these must have same # of nodes/element as existing ones
    void addElements( vector <idx> & m_new, vector <idx>& mat_new ); //adds new elements and materials, assuming element types match existing elements (nodes/per element)

    // creates list of all elements of material numebr mat
    void listElementsOfMaterial(std::vector <idx>& elems, const idx mat) const;
    // Creates list of all nodes belonging to elements of material mat
    void listNodesOfMaterial(std::vector <idx>& nodes, const idx mat) const;
    void listFixLCSurfaces(std::vector <idx>& nodes, const idx FixLCNum) const; // list all nodes of given FixLC surface number (FixLCNum = 1,2,3...)


    bool ContainsAllNodes(idx elem, idx n, idx* nodes) const; // checks if element elem, contains all n nodes in array nodes
    void ContainsNodes(list <idx>* elems , list <idx>* points ); // adds element number to elemes that contain any node in list points

    bool ContainsCoordinate(const idx elem, const double* p, const double* coord) const; // checks whether coordinate is within element elem
    bool ContainsCoordinate(const idx elem, const double* p, const double* coord, double& dist) const; // checks whether coordinate is within element element and also calcultes distance between coord and element centre

    void CompleteNodesSet(const idx elem, std::vector<idx>& nodes) const; // completes nodes vector with those from element, if nodes is empty returns all elements
    void PrintElements()const;       // prints all elements
    void PrintElement(idx e) const;  // prints element e
    void PrintNormals() const;       // prints element surface normals, if present

    void CalculateDeterminants3D(double *p);
   // void FindIndexToMaterialNodes(idx mat, vector<idx> *index) const;
    void CalculateSurfaceNormals(double *p, Mesh* tets = NULL);
    void CopySurfaceNormal(idx i, double* norm); // copies value of surface normal of element i to norm (which must be array of size3)
    void CopyMesh(Mesh* rhs);	// makes this a copy of Mesh* rhs - why does operator= overloading not work???
    void ScaleDeterminants( const double& s); // scales all determinants by s, e.g. to go to microns
    void CalcLocCoords(const idx elem, double* p, double* coord, double* loc); // calculates 4 local coordinates of coordinate cord in element elem
    void CalcElemBary(const idx elem, const double* p, double* bary) const; // calculates barycentre coords of element elem
    double CalcBaryDistSqr( const double *p,const idx elem, const double* coord) const;
    bool isOnXPlane(idx e, double X, double* p) const;	// checks if all nodes of triangle e are on x = X plane
    bool isOnYPlane(idx e, double Y, double* p) const;	// checks if all nodes of triangle e are on y = Y plane
    bool isOnZPlane(idx e, double Z, double* p) const;	// checks if all nodes of triangle e are on z = Z plane
    //bool isOnBackSurface(int e, Geometry* geom);	//
    bool isNeighbours(const idx el1, const idx el2) const; // checks whether el1 and el2 are neighbours
    void gen_p_to_elem(vector<set <idx> >& p_to_elem) const; // generates index from points to mesh elements
    void gen_neighbour_list( vector < vector <idx> >& neigh) const;	// neighbour list of elements. ordered according to local node order.
									// index number of non-existent neighbours elements equals total number of elements
									// i.e. 1 too large to use as an index
    Mesh(idx n,idx nNodes);
    Mesh();
    ~Mesh();
};

#endif



