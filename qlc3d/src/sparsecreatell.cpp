//#include "qlc3d.h"
#include <stdio.h>
#include "sparsematrix.h"
#include "mesh.h"

struct spm{				// define sparce matrix link
		int row;
		spm * next;
		spm * prev;
	};
void sp_elem(spm *col, int r, int c, int *nnz);

void sparsecreatell(Mesh *mesh, SparseMatrix *K)//double **P, int **I, int **J , int *nnz)
{
	
	int maxp	= -10;// number of columns 
	//int *nnz	=(int*)malloc(sizeof(int));
	int i=0,j=0,n1=0,n2=0,x=0; // counters & temporary indicies
	int *tt;		// temporary element

	int dof_per_node =1; // assume this is same for all 
// determine matrix size -find largest node number
	int nnz[1]={0}; // this is stupid, but I am too lazy ...

	
	
	
	
// FIND LARGEST NODE NUMBER IN INPUT MESH	
	int nNodes = mesh->getnNodes(); // nodes per element
	for (int element = 0; element < mesh->getnElements(); element++) //loop number of elements
	{
		for(int node = 0; node < nNodes ; node ++) //loop nodes per element
		{
			//printf("%i %i\n",i,j);
			if (mesh->getNode(element,node) > maxp) 
				maxp = mesh->getNode(element,node);
		}
	}

	maxp++;
	spm *col = new spm[maxp];	// allocate column of pointers
	//for (i=0;i<maxp;i++) aa[i].row =i;
	//printf("maxp = %i\n",maxp);//, aa.?=%i\n",maxp,aa[maxp-1].row);


	for (i=0;i<maxp;i++)  
	{
		//printf("%i\n",i);
		col[i].next = NULL;  // initialize with no elements 
		col[i].prev = NULL;
		col[i].row  = 0;
	}
	//printf("maxp = %i\n",maxp);


// fill in first row of blocks... column blocks are added later
	int r,c;
	tt = (int*)malloc(dof_per_node * sizeof(int));
	// printf("aqda");
	//for (i=0;i<mesh->nElements;i++) //loop over all elements
	 for (int el = 0; el< mesh->getnElements() ; el ++ )
	{
		//tt = &mesh->Elem[i*nNodes]; // shortcut to element nodes 
		//printf("---element %i = [%i,%i,%i,%i]\n",i,tt[0],tt[1],tt[2],tt[3]);
		
		for (x=0;x<dof_per_node;x++){ // this loop is redundant in potential calculation
		//printf("nodes per element = %i \n",mesh->Elements[i].nNodes);
		for (n1=0 ; n1< mesh->getnNodes() ; n1++) // loop over each node in element "i" first loop
		{
			//r = tt[n1] + x*maxp; // sparse row index
			r = mesh->getNode(el,n1) + x*maxp; //sparse row index
			//printf("n1 = %i, r=%i\n",n1,r);
				for (n2=n1+1 ; n2< mesh->getnNodes(); n2++)// loop over each node in element "i" second loop
				{
					
					//c = tt[n2] + x*maxp; // sparse column index
				    c = mesh->getNode(el,n2) + x*maxp;		// column index
					
					
					
					sp_elem(col, mesh->getNode(el,n1), r, nnz); // add row,col n1,n2 nodes to LL
					sp_elem(col, mesh->getNode(el,n1), c, nnz);
					sp_elem(col, mesh->getNode(el,n2), r, nnz);
					sp_elem(col, mesh->getNode(el,n2), c, nnz);
				
				}//end for n2
		}//end for n1
		}//end  for x
	}//end for el
	//printf("maxp = %i\n",maxp);
// extend linked lists, column by column to make matrix square
	/*
	if (dof_per_node>1) // only needed if multiple variables/node
	{
		spm *head, *tail, *temp;
		int rows=0;
		nnz[0] = dof_per_node*nnz[0]; // must increase number of nonzeros
		for (i=0;i<maxp*dof_per_node;i++)
		{
		head = &col[i];
		if (head->next) // if column exists
		{
			head = head->next; // first valid node
			tail = head;
			rows=1;
			while (tail->next!=NULL) // find end of list 
			{
				tail=tail->next; 
				rows++;
			}
			
			for (j=1;j<dof_per_node;j++) // create of_per_node copies of existing LL
			{
				temp = head;
				for (x=0;x<rows;x++)	// make a copy of LL to the end 
				{
					tail->next = new spm;//(spm*)malloc(sizeof(spm));	
					tail->next->prev = tail;				
					tail=tail->next;						
					tail->next = NULL;
					tail->row = (temp->row) + j*maxp; 
					temp=temp->next;
				}//end for x
			}//end for j
	}//end if column exists
	}//end for i
	}//end if 
*/
// ALLOCATE SPARSE MATRIX ARRAYS
//printf("nnz= %i\n",nnz[0]);
//double *Pr = (double*)malloc(nnz[0]*sizeof(double));
int    *Ir = (int*)malloc(nnz[0]*sizeof(int));
int    *Jc = (int*)malloc((maxp+1)*sizeof(int));
//memset(Pr,0,nnz[0]*sizeof(double));
memset(Ir,0,nnz[0]*sizeof(int));
memset(Jc,0,(maxp+1)*sizeof(int));

// FILL IN CONNECTIVITY INFORMATION TO Ir AND Jc ARRAYS FROM LINKED LISTS
// Linked list is freed as it is copied into arrays
	spm *node;
	spm *temp;
	int ir = 0;
	int jc = 0;
	nnz[0]=0; // counts number of entries

	for (i=0;i<maxp*dof_per_node;i++)// loop through col[i];
	{
		node= &col[i];
		node = node->next; // first node in list
		Jc[jc]=nnz[0];		   // column index	
		while (node!=NULL) // follow list untill NULL
		{
			Ir[ir] = (node->row) ;	     // set row number in Ir 	
			nnz[0]++;
			ir++; //next position in Ir array

			//Cleanup, free linked list on the go...
			temp=node;
			node = node->next;
			if (!temp) delete temp;//free(temp); 
					
		}//end while
		
		if (!node) delete node;//free(node);		  // free linked list, so this doesnt need to be done later	
		jc++; // next position in Jc array
	}// end for i
	Jc[jc]=nnz[0]; // last value in Jc is nnz
	delete[] col;
	
	
	K->MakeSparseMatrix(maxp,maxp,nnz[0],Ir,Jc);
	//K->SPY();
	//K->PrintArrays();
	
	//*/
}//end void sprsecreatell


// void sp_elem inserts nodes into linked lists acording to row and column indicies
void sp_elem(spm *col, int r, int c, int *nnz)
{
	//printf("[%i,%i] -> ",r,c);
	spm *node = &col[c]; // was &col[c-1]
	
		
	if (node->next == NULL)//first link in list
	{
		//printf("first link, column %i\n",c);
		node->next = new spm;//(spm*)malloc(sizeof(spm));
		node->next->prev = node;
		node->next->row = r;
		node->next->next = NULL;
		node=node->next;
		nnz[0] ++;
	}

	else
	{
	node = node->next; // new db
	while ((node->row <r) && (node->next!=NULL))
	{
		node = node->next;
	}

	if (node->row==r) //existing entry nothing to be done
	{
		//do nothing 
		//printf("do nothing\n",r,c);
	}
	else if (node->row>r) // must insert new node in middle of list
	{
		//printf("insert in between\n",r,c);
		//spm *temp = new spm;		// new temp node
		spm *temp = new spm;//(spm*)malloc(sizeof(spm));

		node->prev->next = temp;	// insert temp node in between list << insertion fails
		temp->next=node;			//
		temp->prev=node->prev;
		node->prev=temp;			//
		temp->row = r;				// set row value
		nnz[0]++;						// add one to nonzero entries
	}

	else if(node->next==NULL) // found end of column, add new at end
	{
		//printf("add to end\n",r,c);
		//node->next = new spm;
		node->next=new spm;//(spm*)malloc(sizeof(spm));
		node->next->prev = node;	// link backwards
		node->next->next = NULL;	// link forwards = NULL, because end of list
		node->next->row = r;		// set row value
		nnz[0]++;					// add one to nonzero entries
	}
	
	else
	{printf("sparsecreatell.cpp -- !! unknown !!\n");}
	}//end else not first link*/
}

