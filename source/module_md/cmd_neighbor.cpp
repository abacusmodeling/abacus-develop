#include "cmd_neighbor.h"

CMD_neighbor::CMD_neighbor()
{
    dim = 1;
    nlist = new int[dim];
    list = new int*[dim];
    for(int i=0; i<dim; i++)
	{
		list[i] = new int[1];
	}
}

CMD_neighbor::~CMD_neighbor()
{
    delete[] nlist;
    for(int i=0; i<dim; i++)
	{
		delete[] list[i];
	}
    delete[] list;
}

//------------------------------------------------------
// In big cell, for atom i,
// at most one mirror atom of atom j is its adjent atom
//------------------------------------------------------
Vector3<double> CMD_neighbor::Cell_periodic(const Vector3<double> a, const Vector3<double> b)
{
	Vector3<double> temp = a-b;
	while(temp.x<-0.5) temp.x+=1;
	while(temp.x>0.5) temp.x-=1;
	while(temp.y<-0.5) temp.y+=1;
    while(temp.y>0.5) temp.y-=1;
	while(temp.z<-0.5) temp.z+=1;
    while(temp.z>0.5) temp.z-=1;
	return temp;
}

//build the neighbor list
void CMD_neighbor::Neighbor(UnitCell_pseudo &ucell_c)
{
    TITLE("CMD_neighbor", "Neighbor");
    timer::tick("CMD_neighbor", "Neighbor");

    delete nlist;
    for(int i=0; i<dim; i++)
	{
		delete[] list[i];
	}
    delete list;

    dim = ucell_c.nat;
    nlist = new int[dim]();
    list = new int*[dim];
    int *list_temp = new int[dim];

	for(int i=0; i<dim; i++)
	{
        int T1 = ucell_c.iat2it[i];
        int I1 = ucell_c.iat2ia[i];
        Vector3<double> taud1 = ucell_c.atoms[T1].taud[I1];

		for(int j=0; j<dim; j++)
		{
            if(i==j) continue;

            int T2 = ucell_c.iat2it[j];
            int I2 = ucell_c.iat2ia[j];
            Vector3<double> taud2 = ucell_c.atoms[T2].taud[I2];

			Vector3<double> tempd = Cell_periodic(taud1,taud2);
            Vector3<double> temp;
            Mathzone::Direct_to_Cartesian(
			tempd.x, tempd.y, tempd.z,
			ucell_c.latvec.e11, ucell_c.latvec.e12, ucell_c.latvec.e13,
			ucell_c.latvec.e21, ucell_c.latvec.e22, ucell_c.latvec.e23,
			ucell_c.latvec.e31, ucell_c.latvec.e32, ucell_c.latvec.e33,
			temp.x, temp.y, temp.z);

			if(temp.norm()*ucell_c.lat0 <= GlobalV::SEARCH_RADIUS)
			{
				list_temp[nlist[i]]=j;
                nlist[i]++;
			}
		}
        list[i] = new int[nlist[i]];
        for(int index=0; index<nlist[i]; ++index)
        {
            list[i][index] = list_temp[index];
        }
	}

    delete[] list_temp;
    timer::tick("CMD_neighbor", "Neighbor");
    return;
}
