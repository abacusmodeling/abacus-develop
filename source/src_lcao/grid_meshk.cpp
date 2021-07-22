#include "grid_meshk.h"
#include "../src_pw/global.h"

Grid_MeshK::Grid_MeshK()
{
	ucell_index2x = new int[1];
	ucell_index2y = new int[1];
	ucell_index2z = new int[1];
}

Grid_MeshK::~Grid_MeshK()
{
	delete[] ucell_index2x;
	delete[] ucell_index2y;
	delete[] ucell_index2z;
}

int Grid_MeshK::cal_Rindex(const int &u1, const int &u2, const int &u3)const
{
	const int x1 = u1 - this->minu1;
	const int x2 = u2 - this->minu2;
	const int x3 = u3 - this->minu3;
	
	if(x1<0 || x2<0 || x3<0)
	{
		cout << " u1=" << u1 << " minu1=" << minu1 << endl;
		cout << " u2=" << u2 << " minu2=" << minu2 << endl;
		cout << " u3=" << u3 << " minu3=" << minu3 << endl;
		WARNING_QUIT("Grid_MeshK::cal_Rindex","x1<0 || x2<0 || x3<0 !");
	}

	assert(x1>=0);
	assert(x2>=0);
	assert(x3>=0);
	return (x3 + x2 * this->nu3 + x1 * this->nu2 * this->nu3);
}

void Grid_MeshK::cal_extended_cell(const int &dxe, const int &dye, const int &dze)
{
	TITLE("Grid_MeshK","cal_extended_cell");

	//--------------------------------------
	// max and min unitcell in expaned grid.
	//--------------------------------------
	this->maxu1 = dxe / pw.nbx + 1;
	this->maxu2 = dye / pw.nby + 1;
	this->maxu3 = dze / pw.nbz + 1;

	this->minu1 = (-dxe+1) / pw.nbx - 1; 
	this->minu2 = (-dye+1) / pw.nby - 1; 
	this->minu3 = (-dze+1) / pw.nbz - 1; 

	if(test_gridt)OUT(ofs_running,"MaxUnitcell",maxu1,maxu2,maxu3);
	if(test_gridt)OUT(ofs_running,"MinUnitcell",minu1,minu2,minu3);

	//--------------------------------------
	// number of unitcell in each direction.
	//--------------------------------------
	this->nu1 = maxu1 - minu1 + 1;
	this->nu2 = maxu2 - minu2 + 1;
	this->nu3 = maxu3 - minu3 + 1;
	this->nutot = nu1 * nu2 * nu3;

	if(test_gridt)OUT(ofs_running,"UnitCellNumber",nu1,nu2,nu3);
	//xiaohui add 'OUT_LEVEL' line, 2015-09-16
	if(OUT_LEVEL != "m") OUT(ofs_running,"UnitCellTotal",nutot);

//	cout << " nu1 = " << nu1 << " nu2 = " << nu2 << " nu3 = " << nu3 << endl;
//	cout << " nutot = " << nutot << endl;

	delete[] ucell_index2x;
	delete[] ucell_index2y;
	delete[] ucell_index2z;
	this->ucell_index2x = new int[nutot];
	this->ucell_index2y = new int[nutot];
	this->ucell_index2z = new int[nutot];
	ZEROS(ucell_index2x, nutot);
	ZEROS(ucell_index2y, nutot);
	ZEROS(ucell_index2z, nutot);

	this->nutot = nu1 * nu2 * nu3;

	for(int i=minu1; i<=maxu1; i++)
	{
		for(int j=minu2; j<=maxu2; j++)
		{
			for(int k=minu3; k<=maxu3; k++)
			{
				const int cell = cal_Rindex(i,j,k);	
				assert(cell<nutot);

				this->ucell_index2x[cell] = i-minu1;
				this->ucell_index2y[cell] = j-minu2;
				this->ucell_index2z[cell] = k-minu3;

			}
		}
	}

	return;
}
