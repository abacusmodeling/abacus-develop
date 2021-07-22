#include "grid_meshcell.h"
#include "../src_pw/global.h"

Grid_MeshCell::Grid_MeshCell()
{
	allocate_pos = false;
	bx = by = bz = bxyz = 1;
}

Grid_MeshCell::~Grid_MeshCell()
{
	if(allocate_pos)
	{
		for(int ib=0; ib<bxyz; ib++)
		{
			delete[] meshcell_pos[ib];
		}
		delete[] meshcell_pos;
	}
}

void Grid_MeshCell::set_grid_dim(
    const int &ncx_in,
    const int &ncy_in,
    const int &ncz_in,
    const int &bx_in,
    const int &by_in,
    const int &bz_in,
    const int &nbx_in,
    const int &nby_in,
    const int &nbz_in,
    const int &nbxx_in,
    const int &nbzp_start_in,
    const int &nbzp_in
    )
{
    this->ncx = ncx_in;
    this->ncy = ncy_in;
    this->ncz = ncz_in;
    this->ncxyz = ncx * ncy * ncz;
    this->bx = bx_in;
    this->by = by_in;
    this->bz = bz_in;
    this->bxyz = bx*by*bz;
    this->nbx = nbx_in;
    this->nby = nby_in;
    this->nbz = nbz_in;
    this->nbxyz = nbx*nby*nbz;
    this->nbxx = nbxx_in;
    this->nbzp_start = nbzp_start_in;
    this->nbzp = nbzp_in;

	if(nbxx<=0)
	{
		cout << " nbx = " << nbx << endl;
		cout << " nby = " << nby << endl;
		cout << " nbz = " << nbz << endl;
		cout << " nbxx = " << nbxx << endl;
		WARNING_QUIT("Grid_MeshCell::set_grid_dim","nbxx is <=0 !");
	}

	//xiaohui add 'OUT_LEVEL' line, 2015-09-16
	if(OUT_LEVEL != "m") OUT(ofs_running,"real space grid",ncx,ncy,ncz); // real space uniform grid
	if(OUT_LEVEL != "m") OUT(ofs_running,"big cell numbers in grid",nbx,nby,nbz); // reduced by BIG_CELL
	if(OUT_LEVEL != "m") OUT(ofs_running,"meshcell numbers in big cell",bx,by,bz); // is small integer, typical number 2*2*2

	//cout << " bx=" << bx << " by=" << by << " bz=" << bz << endl;
	//cout << " nbx=" << nbx << " nby=" << nby << " nbz=" << nbz << endl; 
	//cout << " ncx=" << ncx << " ncy=" << ncy << " ncz=" << ncz << endl;

    return;
}



// (1)
void Grid_MeshCell::init_latvec(void)
{
	TITLE("Grid_MeshCell","init_latvec");
	// initialize the mesh cell vectors.
	assert(ncx>0);
	assert(ncy>0);
	assert(ncz>0);

	//size of each room (same shape with unitcell)
	this->meshcell_vec1[0]= ucell.a1.x / (double)ncx * ucell.lat0;
	this->meshcell_vec1[1]= ucell.a1.y / (double)ncx * ucell.lat0;
	this->meshcell_vec1[2]= ucell.a1.z / (double)ncx * ucell.lat0;

	this->meshcell_vec2[0]= ucell.a2.x / (double)ncy * ucell.lat0;
	this->meshcell_vec2[1]= ucell.a2.y / (double)ncy * ucell.lat0;
	this->meshcell_vec2[2]= ucell.a2.z / (double)ncy * ucell.lat0;

	this->meshcell_vec3[0]= ucell.a3.x / (double)ncz * ucell.lat0;
	this->meshcell_vec3[1]= ucell.a3.y / (double)ncz * ucell.lat0;
	this->meshcell_vec3[2]= ucell.a3.z / (double)ncz * ucell.lat0;

	this->meshcell_latvec0.e11 = this->meshcell_vec1[0];
	this->meshcell_latvec0.e12 = this->meshcell_vec1[1];
	this->meshcell_latvec0.e13 = this->meshcell_vec1[2];

	this->meshcell_latvec0.e21 = this->meshcell_vec2[0];
	this->meshcell_latvec0.e22 = this->meshcell_vec2[1];
	this->meshcell_latvec0.e23 = this->meshcell_vec2[2];

	this->meshcell_latvec0.e31 = this->meshcell_vec3[0];
	this->meshcell_latvec0.e32 = this->meshcell_vec3[1];
	this->meshcell_latvec0.e33 = this->meshcell_vec3[2];

	// why we need GT = meshcell_latvec0^(-1)?
	// note that (i,j,k) is a grid point.
	// (x,y,z) is the cartesian coordinates.
	// because
	// (x,y,z) = (i,j,k) * meshcell_latvec0
	// once we know (x,y,z) and meshcell_latvec0
	// we need to transform the formula to
	// (x,y,z) * meshcell_latvec0^(-1) = (i,j,k)
	this->meshcell_GT = this->meshcell_latvec0.Inverse();

	if(test_gridt)
	{
		ofs_running << " the VECTORS of MESHCELL are (Bohr): " << endl;
		ofs_running << " vec1( " 
			<< setw(15) << meshcell_vec1[0]
			<< setw(15) << meshcell_vec1[1]
			<< setw(15) << meshcell_vec1[2] 
			<< ")" << endl;

		ofs_running << " vec2( " 
			<< setw(15) << meshcell_vec2[0]
			<< setw(15) << meshcell_vec2[1]
			<< setw(15) << meshcell_vec2[2]
			<< ")" << endl;

		ofs_running << " vec3( " 
			<< setw(15) << meshcell_vec3[0]
			<< setw(15) << meshcell_vec3[1]
			<< setw(15) << meshcell_vec3[2]
			<< ")" << endl;
	}
	
	return;
}



void Grid_MeshCell::init_meshcell_pos(void)
{
	assert(bx>0);
	assert(by>0);
	assert(bz>0);
	assert(bxyz>0);

	if(!allocate_pos)
	{
		meshcell_pos = new double*[bxyz];
		for(int ib=0; ib<bxyz; ib++)
		{
			meshcell_pos[ib] = new double[3];
		}
	}

	int index=0;
	for(int i=0; i<bx; i++)
	{
		for(int j=0; j<by; j++)
		{
			for(int k=0; k<bz; k++)
			{
//				cout << setw(5) << i << setw(5) << j << setw(5) << k;
				for(int p=0; p<3; p++)
				{
					meshcell_pos[index][p] = i*meshcell_vec1[p] + j*meshcell_vec2[p] + k*meshcell_vec3[p];
					//cout << setw(15) << meshcell_pos[index][p]; 
				}
//				cout << endl;
				++index;
			}
		}
	}
	allocate_pos = true;

	return;
}
