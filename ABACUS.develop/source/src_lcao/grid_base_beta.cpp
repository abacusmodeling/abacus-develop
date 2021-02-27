#include "../src_pw/global.h"
#include "ylm.h"

#include "global_fp.h" // mohan add 2021-01-30

Grid_Base_Beta::Grid_Base_Beta()
{ 
	this->nnn = new int[1];
}

Grid_Base_Beta::~Grid_Base_Beta()
{
	delete[] nnn;
}

void Grid_Base_Beta::prepare(
    const Matrix3 &latvec_in,
    const double &lat0_in)
{
	TITLE("Grid_Base_Beta","prepare");

	this->latvec = latvec_in;
	this->lat0 = lat0_in;

	this->latvec0 = this->latvec;
	this->latvec0 *= this->lat0;

	delete[] this->nnn;
	this->nnn = new int[ucell.ntype];
	for(int T1=0; T1<ucell.ntype; T1++)
	{
		this->nnn[T1] = (ucell.atoms[T1].nwl+1) * (ucell.atoms[T1].nwl+1);
	//	cout << "\n nnn = " << nnn[T1];
	}	
	return;
}
