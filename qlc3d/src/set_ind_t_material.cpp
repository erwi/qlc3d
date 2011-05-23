#include "mesh.h"

void Mesh::set_ind_t_material(int MATERIAL, vector<int> *index)
{
	int i=0;
	
	for (i=0;i<nt;i++)
	{
		if (matt[i]==MATERIAL) index->push_back(i);
		}


}