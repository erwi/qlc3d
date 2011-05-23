#include "mesh.h"
//-----------------------------------------------------------------------------------------
void Mesh::show_nodes(int start, int end)
{
	if (np>0)
	{
		if ((start>=0)&&(end<np))
		{
			cout << "node coordinates x-y-z \n";
			for(int i=start;i<=end;i++)
			{
				cout << i<<"      "<< p[i*3+0]<<" "<< p[i*3+1]<<" "<< p[i*3+2] << "\n";
			}
		}
		else
		{
			cout << "Mesh::show_nodes(int start, int end) - check bounds!\n";
		}
	
	}
	else
	{
	cout << "Mesh::show_nodes - no coordinates, np = 0 \n";
	}
}//end void show_nodes;

void Mesh::show_nodes()
{
	show_nodes(0,np-1);
}//end void show nodes;




void Mesh::show_tetrahedra(int start, int end)
{
if (nt>0)
	{
		if ((start>=0)&&(end<nt))
		{
			cout << "tetrahedra nodes and material number\n";
			for(int i=start;i<=end;i++)
			{
				cout << i<<"      "<< t[i*4+0]<<" "<< t[i*4+1]<<" "<< t[i*4+2] <<" "<<t[i*4+3]<<"    "<<matt[i]<< "\n";
			}
		}
		else
		{
			cout << "Mesh::show_tetrahedra(int start, int end) - check bounds!\n";
		}
	
	}
	else
	{
	cout << "Mesh:show_tetrahedra - no tetrahedra, nt = 0 \n";
	}
}//end void show_tetrahedra
void Mesh::show_tetrahedra()
{
	show_tetrahedra(0,nt-1);
}//end void show_tetrahedra()



void Mesh::show_triangles(int start, int end)
{
if (ne>0)
	{
		if ((start>=0)&&(end<ne))
		{
			cout << "triangle nodes and material number \n";
			for(int i=start;i<=end;i++)
			{
				cout << i<<"      "<< e[i*3+0]<<" "<< e[i*3+1]<<" "<< e[i*3+2] <<"    "<<mate[i]<< "\n";
			}
		}
		else
		{
			cout << "Mesh::show_triangles(int start, int end) - check bounds!\n";
		}
	
	}
	else
	{
	cout << "Mesh:show_triangles - no triangles, ne = 0 \n";
	}
}//end void show_tetrahedra
void Mesh::show_triangles()
{
	show_triangles(0,ne-1);
}//end void show_tetrahedra()