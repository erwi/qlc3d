

//
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "../includes/material_numbers.h"
/*
#define MAT_DIELECTRIC1 32
#define MAT_DIELECTRIC2 36
#define MAT_DIELECTRIC3 40
#define MAT_DIELECTRIC4 44
#define MAT_DIELECTRIC5 48
#define MAT_DIELECTRIC6 52
#define MAT_DIELECTRIC7 56
*/




int ReorderNodes(double *p, int np,int *t, int nt, int *e, int ne,int *tmat,int *emat)
{// reorder mesh so that LC elements come first. Returns number of LC nodes npLC

    // check if dielectric elements exist
    bool DE_exist = false;
    int npLC = np;

    for (int i = 0 ; i <nt ; i++)
        if (tmat[i]>= MAT_DIELECTRIC1)	// if material nuber > LC number
        {
            DE_exist = true;
            break;
        }


    if (DE_exist) // dielectric materials exist, must reorder
    {
        printf("Reordering nodes...");
        struct Node			//Define linked list node
        {
            Node *up;
            Node *down;
            int value;
        };
        int i;
        double *newp;		//reordered coordinates
        int *newt;			//reordered tetrahedra
        int *newe;			//reordered triangles
        //int *newperi;		//reordered periodic nodes

        Node *head;			// linked list node pointers
        Node *tail;
        Node *tn;

        tn = new Node;		// initialise empty linked list
        tn->value = 0;
        tn->up = NULL;
        tn->down = NULL;
        head = tn;
        tail = tn;

        //GENERATE LIST OF MATERIAL NUMBERS FOR NODES. IN CASE OF DUAL VALUE LC PRECEDES
        int *mat_index = (int*)malloc(np*sizeof(int)); // array of node materials

        for (i = 0 ; i < np ; i++ )
            mat_index[i] = MAT_DIELECTRIC1;

        for (i=0;i<nt;i++)//iterate and re-set all LC nodes back to LC material
        {
            if (tmat[i]< MAT_DIELECTRIC1)
            {
                mat_index[t[i*4+0]]=tmat[i];
                mat_index[t[i*4+1]]=tmat[i];
                mat_index[t[i*4+2]]=tmat[i];
                mat_index[t[i*4+3]]=tmat[i];
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

            if (mat_index[i]< MAT_DIELECTRIC1) //if LC -> front of list
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
        }//end for linked list assembly
        npLC++;
        // CONVERT LINKED LIST TO ARRAY AND FREE LINKED LIST
        i=0;
        while(head->down!=NULL)
        {

            //mat_index[i] = head->value;
            mat_index[head->value] = i;
            i++;
            delete head->up;
            if (head->down!=NULL)head = head->down;
        }
        mat_index[head->value] = i;
        //mat_index[i] = head->value;
        delete head;
        printf("\n reordered mat_index:\n");
        for (i = 0 ; i < np ; i++ )	printf("mat_index[%i] = %i\n",i,mat_index[i]);
        // index is now formed. use this to rearrange mesh numbering

        //REORDER NODES
        newp = (double*)malloc(3*np*sizeof(double)); // memory for reordered node coordinates

        //reorder nodal coordinates
        for (i=0;i<np;i++)
        {
            newp[mat_index[i]*3+0] = p[i*3+0];
            newp[mat_index[i]*3+1] = p[i*3+1];
            newp[mat_index[i]*3+2] = p[i*3+2];
        }

        //copy back to original array
        for (i=0;i<3*np;i++)
            p[i]=newp[i];

        free(newp);

        /// REORDER TETRAHEDRA

        newt = (int*)malloc(4*nt*sizeof(int));

        //reorder t nodes
        for (i=0;i<nt*4;i++)
            newt[i]= mat_index[t[i]];

        //copy back to original array
        for (i = 0 ; i <4*nt ; i++)
            t[i] = newt[i];

        free(newt);


        /// REORDER TRIANGLES

        newe = (int*)malloc(3*ne*sizeof(int));

        //reorder e nodes
        for(i=0;i<3*ne;i++)
            newe[i] = mat_index[e[i]];

        // copy back to original array
        for(i=0;i<3*ne;i++)
            e[i] = newe[i];

        free(newe);
        free(mat_index);

        printf("OK\n");
        return npLC;
    }

    else	// no dielectric elements - no need for reordering
    {
        return npLC;
    }
}



