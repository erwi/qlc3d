# include "mesh.h"
# include <string.h>

//		READ GID MESH
//
void Mesh::read_GID_mesh(string filename)
{
	
	char char_filename[30];
	char charray[200];
	int temp,i;

	strcpy(char_filename, filename.c_str());
	
	ifstream fin;
	fin.open(char_filename);
	
	if (fin.fail())
	{
	cout << "file open failed!\n";
	exit(1);
	}
	else
	{
	cout << "reading GID mesh: "<< filename << "\n";
	
	fin.getline(charray,200,'\n'); // READ 2 LINES OF JUNK
	fin.getline(charray,200,'\n');

//COUNT NODES
	while (strcmp(charray,"Coordinates")) fin.getline(charray,200,'\n')	;
	np = 0;
	while (strcmp(charray,"end coordinates"))
	{
			fin.getline(charray,200,'\n');
			//cout << charray<<"\n";
			np++;
	}
	np--;
	cout << "# nodes : "<< np <<"\n";

//COUNT TETRAHEDRA	
	fin.getline(charray,200,'\n'); // READ 2 LINES OF JUNK
	fin.getline(charray,200,'\n');
	nt = 0;
	while (strcmp(charray,"end elements"))
		{
		fin.getline(charray,200,'\n');
		nt++;
		}
		nt--;
		cout << "# tetrahedra : " << nt <<"\n";

//COUNT PERIODIC NODES
		nperi=0;
		bIsPeriodic = false;
		fin.getline(charray,200,'\n');
	
		//test if periodic nodes exist
		if(!strcmp(charray,"MESH    dimension 3 ElemType Prisma  Nnode 6")){
			bIsPeriodic = true;
			while (strcmp(charray,"Elements"))	{fin.getline(charray,200,'\n');	}
			while (strcmp(charray,"end elements"))
				{	fin.getline(charray,200,'\n');	
					nperi++;
				}
				nperi--;
				cout << "# periodic : " << nperi <<"\n";
		}//end if periodic nodes exist


//COUNT TRIANGLES
	while (strcmp(charray,"Elements"))	{fin.getline(charray,200,'\n');	}

	ne = 0;
	while (strcmp(charray,"end elements"))
		{
		fin.getline(charray,200,'\n');
		ne++;
		}
		ne--;
		cout << "# triangles : " << ne <<"\n";

//ALLOCATE MEMORY FOR MESH DATA
	p = (float*)malloc(3*np*sizeof(float));
	t = (int*)  malloc(4*nt*sizeof(int));
	e = (int*)  malloc(3*ne*sizeof(int));
	matt = (int*)malloc(nt*sizeof(int));
	mate = (int*)malloc(ne*sizeof(int));
	if (bIsPeriodic)								//allocate if periodc nodes exist
	{		peri = (int*)malloc(6*nperi*sizeof(int));	}
	else	{peri = NULL;}


//REWIND TO START OF FILE
	fin.clear();              // forget we hit the end of file
	fin.seekg(0, ios::beg);   // move to the start of the file


//READ NODE COORDINATES
	cout << "Reading node coordinates ...";
	fin.getline(charray,200,'\n'); // READ 2 LINES OF JUNK
	fin.getline(charray,200,'\n');
	for (i =0;i<np;i++)
	{
		fin >> temp; //node number - not needed
		fin >> p[i*3+0];
		fin >> p[i*3+1];
		fin >> p[i*3+2];
	}
	cout << "OK\n";

// READ TETRAHEDRA NODE NUMBERS AND MATERIAL NUMBER
	cout << "Readin tetrahedra ...";
	fin.getline(charray,200,'\n');
	while (strcmp(charray,"Elements"))	fin.getline(charray,200,'\n');	//find start of tet data	
	for (i =0;i<nt;i++)
	{
		fin >> temp; //tetrahedra number - not needed
		fin  >> t[i*4+0];
		fin  >> t[i*4+1];
		fin  >> t[i*4+2];
		fin  >> t[i*4+3];
		fin  >> matt[i];
	 //if (i<10) cout << i << t[i*4] << matt[i] <<"\n";
	}
	cout << "OK\n";



//READ PERIODIC NODES (IF EXIST)
	if (bIsPeriodic)
	{
		
	fin.getline(charray,200,'\n');
	
		while (strcmp(charray,"Elements")) fin.getline(charray,200,'\n'); //find start of periodic nodes
		cout << "Reading periodic nodes ...";
		
	
	for (i=0;i<nperi;i++)
	{
		fin >> temp;
		fin >> peri[i*6+0];
		fin >> peri[i*6+1];
		fin >> peri[i*6+2];
		fin >> peri[i*6+3];
		fin >> peri[i*6+4];
		fin >> peri[i*6+5];
	//	fin >> temp;
	
	
	}
  
	cout << "OK\n";
	}//end if read periodic nodes

//READ TRIANGLE NODE NUMBERS AND MATERIAL NUMBER

	fin.getline(charray,200,'\n');
//	cout << charray <<"\n";

	while (strcmp(charray,"Elements"))	fin.getline(charray,200,'\n');		//find start of triangle data	
	
		cout << "Reading triangles ...";
		for (i =0;i<ne;i++)
	{
		fin >> temp; //triangle number - not needed
		fin >> e[i*3+0];
		fin >> e[i*3+1];
		fin >> e[i*3+2];
		fin >> mate[i];
	}
	cout << "OK\n";





	fin.close(); // close file 

// ONLY POSITIVE COORDINATES ALLOWED - MOVE IF NECESSARY	
	xmin = 1e9 ; ymin = 1e9 ; zmin = 1e9;
	xmax = -1e9; ymax = -1e9; zmax = -1e9;
	for (i=0;i<np;i++)
	{
		if(p[i*3+0]<xmin) xmin = p[i*3+0]; // find min
		if(p[i*3+1]<ymin) ymin = p[i*3+1];
		if(p[i*3+2]<zmin) zmin = p[i*3+2];

		if(p[i*3+0]>xmax) xmax = p[i*3+0]; // find max
		if(p[i*3+1]>ymax) ymax = p[i*3+1];
		if(p[i*3+2]>zmax) zmax = p[i*3+2];
	}
	for (i=0;i<np;i++)
	{
		if(xmin<0) p[i*3+0]-=xmin; // move mesh 
		if(ymin<0) p[i*3+1]-=ymin;
		if(zmin<0) p[i*3+2]-=zmin;
	}
	if (xmin<0) {xmax-=xmin; xmin=0; cout << "shifting x\n";}
	if (ymin<0) {ymax-=ymin; ymin=0; cout << "shifting y\n";}
	if (zmin<0) {zmax-=zmin; zmin=0; cout << "shifting z\n";}

//	cout << xmin << xmax << ymin << ymax<< zmin << zmax<<'\n';
// REMOVE NUMERICAL NOISE FROM BOUNDARY NODES

	for (i=0;i<np;i++)
	{
		if( p[i*3+0]-xmin <=1e-9) p[i*3+0]=xmin;
		if( xmax-p[i*3+0] <=1e-9) p[i*3+0]=xmax;

		if( p[i*3+1]-ymin <=1e-9) p[i*3+1]=ymin;
		if( ymax-p[3*i+1] <=1e-9) p[i*3+1]=ymax;
		
		if( p[i*3+2]-zmin <=1e-9) p[i*3+2]=zmin;
		if( zmax-p[3*i+2] <=1e-9) p[i*3+2]=zmax;
	}




	set_ind_t_material(4,&ind_t_LC);
	set_ind_t_material(36,&ind_t_DE);
	cout << "LC tets: "<< ind_t_LC.size() << " DE tets: "<< ind_t_DE.size() <<"\n" ;

	if (ind_t_DE.size()>0)
	{
	cout << "reordering nodes ...";
	reorder_nodes();
	cout << "OK\n";
	}



	}//end if mesh file opened ok




}// end void read_GID_mesh
