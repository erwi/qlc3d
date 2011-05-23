#include "mesh.h"

/*
%returns index to periodic nodes 
%
%       corn3   nright
%       *---------------*corn2
%       |               |
%       |               |
%       |nleft          | nright
%       |               |        
%       |               |
%       *---------------*corn1
%     corn0    nleft
%
%
*/




bool compare(int i, int j)
{
return (i==j);
}


void Mesh::sort_by_coordinate(vector<int> *index, int c)
{
// c = 0,1,2 for x,y or z coordinate
	int i,j,i_sm;
	float sm,te;
	float *tempp;
	
	//make temporary list of coordinttes to sort
	tempp =(float*) malloc(index->size()*sizeof(float));
	vector<int>::iterator itr;
	i=0;
//	cout <<"------\n";
	for(itr=index->begin();itr<index->end();itr++)
	{
		tempp[i]=p[(*itr)*3+c];
//		cout <<i<<" "<< tempp[i] <<'\n';
		i++;
	}
//SORT	
	for (i=0;i<index->size();i++)
	{
		//find smallest
		sm = 1e9;//tempp[i];
		i_sm = i;
		for (j=i;j<index->size();j++)	
		{
			if (tempp[j]<sm)
			{
				sm=tempp[j];
				i_sm=j; 
			}//find index to smallest
		}//end for j
		//swap
		if(i_sm>i)
		{
			te = tempp[i];
			tempp[i] = tempp[i_sm];
			tempp[i_sm] = te;
	
			j=*(index->begin()+i);
			*(index->begin()+i)=*(index->begin()+i_sm);
			*(index->begin()+i_sm)=j;
		}
	}//end for i

//	for (i=0;i<index->size();i++)
//	cout << *(index->begin()+i) <<" "<<tempp[i] <<"\n";	

	free(tempp);

}//end void sort_by_coordinate








void Mesh::find_periodic_nodes()
{
	int i;

	vector<int> i_perie;
	vector<int> i_perip;
	vector<int> temp;

	cout << "finding periodic nodes \n";
//FIND INDEX TO ALL PERIODIC TRIANGLES
	for (i=0;i<ne;i++){
		if (mate[i]==3) 
		{
			i_perie.push_back(i);
		}
		else
		{
	//	cout << i << " " <<mate[i] << "\n";
		}
	}
	cout << "periodic triangles #:" << i_perie.size() <<"\n"; 
	if (i_perie.size()>0)
	{
//FIND INDEX TO ALL PERIODIC NODES

	for (i=0;i<i_perie.size();i++)
	{
		i_perip.push_back(e[i_perie[i]*3+0]-1);
		i_perip.push_back(e[i_perie[i]*3+1]-1);
		i_perip.push_back(e[i_perie[i]*3+2]-1);
	}
	sort(i_perip.begin(),i_perip.end());
	vector<int>::iterator itr1; 
	itr1 = unique_copy(i_perip.begin(),i_perip.end(),i_perip.begin(),compare);
    i_perip.resize( itr1 - i_perip.begin() ); // <-- is now index to all periodic nodes

//---------------------------------------------------
//
//         FIND INDEX TO ALL CORNER NODES
//
//---------------------------------------------------
	

	for (i=0;i<i_perip.size();i++)
	{
	//determine corner based on coordinate
		//xmin and ymin -> corn0
    	if ((p[3*i_perip[i]+0]==xmin)&&(p[3*i_perip[i]+1]==ymin))
		{
			c0.push_back(i_perip[i]);
			i_perip[i]=-1;
		}
		
		//xmin and ymax -> corn3
		else if ((p[3*i_perip[i]+0]==xmin)&&(p[3*i_perip[i]+1]==ymax))
		{
			c3.push_back(i_perip[i]);
			i_perip[i]=-1;
		}
		//xmax and ymin -> corn1
		else if ((p[3*i_perip[i]+0]==xmax)&&(p[3*i_perip[i]+1]==ymin))
		{
			c1.push_back(i_perip[i]);
			i_perip[i]=-1;
		}
		//xmax and ymax -> corn2
		else if ((p[3*i_perip[i]+0]==xmax)&&(p[3*i_perip[i]+1]==ymax))
		{
			c2.push_back(i_perip[i]);
			i_perip[i]=-1;
		}

	  }//end for
	
		if ((c0.size()!=c1.size()) || (c0.size()!=c2.size()) || (c0.size() != c3.size()))
		{
			cout << "Mesh::find_periodic_nodes() - inequal number of periodic corner nodes !!\n";
			exit(1);
		}


		
//STILL NEED TO SORT CORNER NODE LISTS ACCORDING TO Z-COORDINATE 
		sort_by_coordinate(&c0,2);
		sort_by_coordinate(&c1,2);
		sort_by_coordinate(&c2,2);
		sort_by_coordinate(&c3,2);
		

		//remove corner nodes (set to -1) from index list
		sort(i_perip.begin(),i_perip.end()); // all -1's go first
		itr1 = unique_copy(i_perip.begin(),i_perip.end(),i_perip.begin(),compare); //remove repetitions (-1's)
		i_perip.resize( itr1 - i_perip.begin() ); 
		i_perip.erase(i_perip.begin()); // erase first entry (the remaining -1)

//---------------------------------------------------------------------------
//
//			THEN FIND NLEFT AND NRIGHT NODES
//
//---------------------------------------------------------------------------
	
		//separate remaining nodes to four planes : left, right, front and back

		//vector <int> nleft;
		//vector <int> right;
		vector <int> front;
		vector <int> back;

		vector <int>::iterator itr;
		for (itr=i_perip.begin();itr<i_perip.end();itr++)
		{
			if (p[(*itr)*3+0]==xmin) nleft.push_back(*itr);
			else if (p[(*itr)*3+0]==xmax) nright.push_back(*itr);
			else if (p[(*itr)*3+1]==ymin) front.push_back(*itr); 
			else if (p[(*itr)*3+1]==ymax) back.push_back(*itr);
			else
			{
				cout << "Mesh::find_periodic_nodes - incromulent periodic node (not on side)!!\n";
				exit(1);
			}
		}

		//make sure left/right and front/back have same number of nodes
		if (nleft.size()!=nright.size())
		{
			cout << "Mesh>>find_periodic_nodes - inequal number of nodes on left/right periodic boundaries\n";
			exit(1);
		}
		if (front.size()!=back.size())
		{
			cout << "Mesh>>find_periodic_nodes - inequal number of nodes on front/back periodic boundaries\n";
			exit(1);
		}

		//reorder index to right boundary to match left boundary
	
		int temp_int,j,l,r;
		for (i=0;i<nleft.size();i++)
		{
			l = nleft[i];
			for (j=0;j<nright.size();j++)
			{
				r = nright[j];
				//if match found, swap
				if ((p[l*3+1]==p[r*3+1]) && ((p[l*3+2]==p[r*3+2])))
				{
					temp_int = nright[j];
					nright[j] = nright[i];
					nright[i] = temp_int;
					break;
				}
				if (j==nright.size()-1) 
				{
					cout << "Mesh::find_periodic_nodes -did not find periodic left/right node \n";
					cout << "(change coordinate comparison tolerance)\n"	;
					exit(1);
				}
			}//end for j
		}//end for i


//reorder index to front boundary to match back boundary
	
		for (i=0;i<front.size();i++)
		{
			l = front[i];
			for (j=0;j<back.size();j++)
			{
				r = back[j];
				//if match found, swap
				if ((p[l*3+0]==p[r*3+0]) && ((p[l*3+2]==p[r*3+2])))
				{
					temp_int = back[j];
					back[j] = back[i];
					back[i] = temp_int;
					break;
				}
				if (j==back.size()-1) 
				{
					cout << "Mesh::find_periodic_nodes -did not find periodic front/back node \n";
					cout << "(change coordinate comparison tolerance)\n"	;
					exit(1);
				}
			}//end for j
		}//end for i

// COPY FRONT TO NLEFT AND BACK TO NRIGHT
for (i=0;i<front.size();i++)
{
	nleft.push_back(front[i]);
	nright.push_back(back[i]);
}

}//end if periodic sides exist

}//end void find_periodic nodes