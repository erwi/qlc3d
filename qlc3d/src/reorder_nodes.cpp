# include "mesh.h"

void Mesh::reorder_nodes()
{

    struct Node			//doubly linked list
	{
		Node *up;
		Node *down;
		int value;

	};

	int i;
	float *newp;		//reordered coordinates
	int *newt;			//reordered tetrahedra
	int *newe;			//reordered triangles
	int *newperi;		//reordered periodic nodes

	
	Node *head;
	Node *tail;
	Node *tn;
	
	

	tn = new Node;
	tn->value = 0;
	tn->up = NULL;
	tn->down = NULL;
	head = tn;
	tail = tn;
	

//GENERATE LIST OF MATERIAL NUMBERS FOR NODES. IN CASE OF DUAL VALUE LC PRECEDES	
	int *mat_index = (int*)malloc(np*sizeof(int)); // array of node materials
	memset(mat_index,36,np*sizeof(int)); // set all to DE to start		
	
	
	for (i=0;i<nt;i++)//iterate and re-set all LC nodes back to LC material
	{
		if (matt[i]==4)
		{
			mat_index[t[i*3+0]]=4;
			mat_index[t[i*3+1]]=4;
			mat_index[t[i*3+2]]=4;
			
		}
	}

//ASSEMBLE DOUBLY LINKED LIST, SORTED BY MATERIAL NUMBERS FROM NODE MATERIALS
int counter = 0;
	npLC=0;
	for (i=1;i<np;i++)
	{
		
		tn = new Node;
		tn->value = i;
		tn->up=NULL;
		tn->down=NULL;
		
			
			
			if (mat_index[i]==4) //if LC -> front of list
			{
				tn->down = head;	//point to previous first
				head->up = tn;		
				head = tn;
				counter ++;
				npLC++;
			}
			else // else DE -> back of list
			{
				tn->up = tail; //add to last position
				tail->down =tn;
				tail = tn;
				counter ++;			
			}
	//cout << tail->value << '\n';
	
	
	}//end for
	//cout <<"counter: "<< counter;

	
// CONVERT LINKED LIST TO ARRAY
	i=0;	
	while(head->down!=NULL)
	{
		mat_index[i] = head->value;
		i++;
		delete head->up;
		if (head->down!=NULL)head = head->down;
	}

	for (i = 0 ; i < np ; i ++ )
		printf("mat_index[%i] = %i\n", i , mat_index[i]);
	//cout <<"counter" << i << '\n';


	newp = (float*)malloc(3*np*sizeof(float)); // memory for reordered node coordinates

//REORDER NODES 
	for (i=0;i<np;i++)
	{
	newp[i*3+0]=p[mat_index[i]*3+0]; //x-coord
	newp[i*3+1]=p[mat_index[i]*3+1]; //y-coord
	newp[i*3+2]=p[mat_index[i]*3+2]; //z-coord
	}
	
	free(p); // make p = new reordered p
	p=newp;

//REORDER TETRAHEDRA
	newt = (int*)malloc(4*nt*sizeof(int));

	for (i=0;i<nt;i++)
	{
	newt[i*4+0]= mat_index[t[i*4+0]-1];
	newt[i*4+1]= mat_index[t[i*4+1]-1];
	newt[i*4+2]= mat_index[t[i*4+2]-1];
	newt[i*4+3]= mat_index[t[i*4+3]-1];
	}

	free(t);
	t=newt;

//REORDER TRIANGLES
	newe = (int*)malloc(3*ne*sizeof(int));

	for(i=0;i<ne;i++)
	{
		newe[i*3+0] = mat_index[e[i*3+0]-1];
		newe[i*3+1] = mat_index[e[i*3+1]-1];
		newe[i*3+2] = mat_index[e[i*3+2]-1];

	}
	free(e);
	e=newe;

//REORDER PERIODIC NODES
	if (nperi>0)
	{
		newperi = (int*)malloc(nperi*6*sizeof(int));
		for (i=0;i<nperi;i++)
		{
		newperi[i*6+0] = mat_index[peri[i*6+0]-1];
		newperi[i*6+1] = mat_index[peri[i*6+1]-1];
		newperi[i*6+2] = mat_index[peri[i*6+2]-1];
		newperi[i*6+3] = mat_index[peri[i*6+3]-1];
		newperi[i*6+4] = mat_index[peri[i*6+4]-1];
		newperi[i*6+5] = mat_index[peri[i*6+5]-1];
		}
	free(peri);
	peri = newperi;
	
	}



}
