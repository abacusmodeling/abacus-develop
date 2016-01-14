//=========================================================
//AUTHOR : liaochen
//DATE : 2008-03-04
//=========================================================
#include "numerical_nonlocal.h"

Numerical_Nonlocal::Numerical_Nonlocal()
{
	//make pair of new and delete
	//question remains
	this->type = 0;
	this->lmax = 0;
	this->LfromBeta = new int[1];
	this->Proj = new Numerical_Nonlocal_Lm[1];
	this->nproj = -1;
}

Numerical_Nonlocal::~Numerical_Nonlocal()
{
	delete[] Proj;
	delete[] LfromBeta;
}

void Numerical_Nonlocal::set_type_info
(
	const int& type_in,
	const string& label_in,
	const string& type_ps_in,
	const int& lmax_in,
	matrix& Coefficient_D_in,
	const int& nproj_in,
	int* LfromBeta_in,
	const Numerical_Nonlocal_Lm* Proj_in
)
{
	if (type_in < 0 || type_in > 2)
	{
		WARNING("Numerical_Nonlocal", "bad input of type_in: not ready yet for type >2");
	}

	this->type = type_in;
	this->label = label_in;
	this->type_ps = type_ps_in;

	if (lmax_in < -1 || lmax_in > 20)
	{
		WARNING_QUIT("Numerical_Nonlocal", "bad input of lmax : should be between -1 and 20");
	}

	this->lmax = lmax_in;
//----------------------------------------------------------
//EXPLAIN : Coefficient D used in calculate elements of NLps
//----------------------------------------------------------
	this->Coefficient_D.create( lmax_in+1, lmax_in+1);
	for (int L1 = 0; L1 < lmax + 1; L1++)
	{
		for (int L2 = 0; L2 < lmax + 1; L2++)
		{
			this->Coefficient_D(L1, L2) = Coefficient_D_in(L1, L2);
		}
	}

//----------------------------------------------------------
//EXPLAIN : LfromBeta
//----------------------------------------------------------
	this->nproj = nproj_in;
	assert(nproj <= lmax_in+1); //LiuXh 2016-01-13
	assert(nproj >= 0);

	delete[] LfromBeta;
	this->LfromBeta = new int[nproj];
	ZEROS(LfromBeta, nproj);
	for(int i=0; i<nproj; i++)
	{
		this->LfromBeta[i] = LfromBeta_in[i];
	}

//----------------------------------------------------------
// EXPLAIN : non_local pseudopotential projector for each l
//----------------------------------------------------------
	//only store radial function
	delete[] Proj;
	this->Proj = new Numerical_Nonlocal_Lm[this->nproj];

	for (int p1=0; p1<nproj; p1++)
	{
		this->Proj[p1] = Proj_in[p1];
	}

	this->rcut_max = 0.0;
	for(int p1=0; p1<nproj; p1++)
	{
		this->rcut_max = max( this->Proj[p1].getRcut(), rcut_max ); 
	}
	return;
}

