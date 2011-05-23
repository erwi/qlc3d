
#include "mesh.h"

void Mesh::write_mesh(string filename)
{

	ofstream fout;
	char ch_filename[100], chstr[100];
	
	string str;
	int i,chlen;



	strcpy(ch_filename,filename.c_str());
	fout.open(ch_filename);

	if (fout.fail()) // could not open file
	{
	cout << "Mesh::write_mesh() - file open failed!\n";
	exit(1);
	}
//mesh file opened OK	
	else	
	{
		if((np>0)&&(nt>0)) // test if mesh is initialised
		{
			fout << "MESH    dimension 3 ElemType Tetrahedra  Nnode 4\n";
			fout << "Coordinates\n";
		
//write coordinates
			for (i=0;i<np;i++)
			{
				
				chlen = sprintf(chstr,"\t %d \t %f \t %f \t %f\n",i+1, p[i*3+0],p[i*3+1],p[i*3+2]);
				fout <<	chstr;
			}			
			fout <<"end coordinates\n\nElements\n";

//write tetrahedra elements
			for(i=0;i<nt;i++)
			{
			    chlen = sprintf(chstr,"\t %d \t %d \t %d \t %d \t %d \t %d \n", i+1, t[i*4+0], t[i*4+1], t[i*4+2], t[i*4+3], matt[i]);
				fout << chstr;
			}
			fout << "end elements\nMESH    dimension 3 ElemType Triangle  Nnode 3\nCoordinates\nend coordinates\n\nElements\n";
//write triangle elements
			for(i=0;i<ne;i++)
			{
				chlen = sprintf(chstr,"\t%d \t %d \t %d \t %d \t %d \n",i+1,e[i*3+0],e[i*3+1],e[i*3+2], mate[i]);
				fout << chstr;
			}
			fout << "end elements\n";



			fout.close();
		}
		else
		{
		cout << "Mesh::write_mesh() - no mesh to write, np or nt = 0 \n";
		fout.close();
		exit(1);
		}
			
		
	
	
	
	}//end if file opened OK



}//end void