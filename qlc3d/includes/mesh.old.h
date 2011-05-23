
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <algorithm>
using std::vector;

class Mesh 
{
public:
int Dimension;
int nElements;
int nNodes;

int *Elem;
int *Mat;


Mesh(int n,int nNodes);


//========================================
//
//
//
//
//
//
/*
void SetAllNodes(int *nodes, int npt)
{
	printf("setting nodes %i\n",npt);
	/*
	if ((nElements>0) && (Elem!=NULL)) // make sure mesh is initilised
	{
		for (int i = 0 ; i < nElements*nNodes ; i++)
				Elem[i] = nodes[i];
	}
	else
	{
		printf("error in Mesh::SetAllNodes - mesh not initilised\n");
		exit(1);
	}
	
}// end void SetAllNodes
/*
void SetAllMaterials(int *mat)
{
	int i;
	if (nElements>0)
	{
		for (i=0;i<nElements;i++)
		{
			Elements[i].Material = mat[i];
		}
	}
	else
	{
		printf("error in Mesh::SetAllMaterials - mesh not initilised\n");
		exit(1);
	}
}// end void SetAllMaterials

void PrintElements()
{
int i,j;
	
	for (i=0;i<nElements;i++)
	{
		printf("%i -- ",i);
		for(j=0;j<nNodes;j++)
			printf(" %i",Elem[i*nNodes +j]);
		
		printf("- mat = %i\n",Mat[i]);
	}
	
}

/*
void FindElementsMaterial(int Material, vector<int> *ind)
{
// Returns index to all elements of material Material, zero indexing
	int i;
	if (nElements>0)
	{
		for (i=0;i<nElements;i++)
		{
			if (Elements[i].Material==Material)
			{
			
				ind->push_back(i);
			}
		}
	
	}
	else
	{
		printf("mesh.h - elements not initialised");
		exit(1);
	}
}// end void FindElementsMaterial

void FindNodesMaterial(int Material, vector<int> *indp)
{
// Returns index to all nodes of material Material, same indexing as is stored in Mesh, i.e. not necessarily zero for first node
	if (nNodes>0)
	{
	// 1. find all elements of type Material
		vector<int> indt;
		FindElementsMaterial(Material, &indt);
	
	//2. extract all nodes from vector indt and store in vector indp
		int j;
		vector<int>::iterator i;
		vector<int> temp;
		for (i=indt.begin();i != indt.end();i++)
		{	
			for (j=0 ; j < Elements[*i].nNodes;j++)
				temp.push_back(Elements[*i].Nodes[j]);
		}
	// 3. remove duplicate entries
		std::sort(temp.begin(),temp.end()); // first need to sort 
		vector<int>::iterator new_end;
		new_end	= std::unique(temp.begin(),temp.end()); // reorder unique entries to start of vector 
		
		for (i=temp.begin() ; i != new_end; i++)	// copy only unique entries to output vector 
			indp->push_back(*i);
	
	}
	else
	{
		printf("mesh::FindNodesMaterial - mesh not initialised");
		exit(1);
	}
	
}// end void FindElementsMaterial



*/
// Create mesh consistin gof n elements, each of nNodes number of nodes
//Mesh(int n,int nNodes);
/*
{
	printf("Mesh.h\n");
	nElements = n;
	// allocate memory for mesh
	Elem = (int*)malloc(n * nNodes * sizeof(int));
	Mat  = (int*)malloc(n * sizeof(int));
	
	
}

~Mesh()//Destructor
{
		
	if (Elem != NULL) free(Elem);
	if (Mat  != NULL) free(Mat);

	}
/*
Mesh()
{
	nNodes=0; 
	nElements= 0; 
	Dimension = 0;
	Elem= NULL;
	Mat = NULL;
}
*/
};

