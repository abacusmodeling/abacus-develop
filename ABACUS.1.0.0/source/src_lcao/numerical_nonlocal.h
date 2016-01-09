//=========================================================
//AUTHOR : liaochen, mohan
//DATE : 2009-03-04
//=========================================================
#ifndef NUMERICAL_NONLOCAL_H
#define NUMERICAL_NONLOCAL_H

#include "../src_pw/tools.h"
#include "numerical_nonlocal_lm.h"
//=========================================================
//CLASS  Numerical_Nonlocal
//Note : contain nonlocal ps(:pseudopotential) information 
//about atoms
//Feature : set and store information about ps infomation 
//			related to atoms
//=========================================================

class Numerical_Nonlocal
{
public:

	Numerical_Nonlocal();
	~Numerical_Nonlocal();

//==========================================================
// EXPLAIN : get information from ps orbitals
// MEMBER FUNCTION :
//===========================================================
	
	// mohan add 2011-03-07
	const int& getL_Beta(const int &nb){ return this->LfromBeta[nb]; }
	const int& getLmax() const { return this->lmax; }
   	const int& getType() const { return this->type; }
	const string& getLabel() const { return this->label; }
	const string& getType_ps() const { return this->type_ps; }
	

	const double& getCoefficient_D(const int& L1, const int& L2) 
						const { return this->Coefficient_D(L1, L2); }
	
//==========================================================
// EXPLAIN : set information about Numerical Orbital
// why double** cannot have const in the front???
//==========================================================
	void set_type_info(
			const int& type_in,
			const string& label_in,
			const string& type_ps_in,
			const int& lmax_in,
			matrix& Coefficient_D_in,
			const int& nproj_in,
			int* LfromBeta_in,
			const Numerical_Nonlocal_Lm* ps_orbital_in);

	Numerical_Nonlocal_Lm* Proj;// length: nproj(only store radial function )

	const double& get_rcut_max(void) const { return rcut_max; }

	private:
	
//==========================================================
// MEMBER FUNCTION :
// NAME : label (atom type)
// NAME : lmax( max value of L angular momentum) 
//===========================================================
	string label;
	int type;
	int lmax;
	double rcut_max;
	string type_ps; //local or nonlocal
	matrix Coefficient_D; //Coefficient Matrix dimension : (lmax+1) * (lmax+1)
	
	// mohan add 2011-03-07
	// each Beta may have different L.
	int nproj;
	int *LfromBeta;
};

#endif
