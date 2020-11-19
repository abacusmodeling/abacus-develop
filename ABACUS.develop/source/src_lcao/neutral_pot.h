#ifndef NEUTRAL_POT_H
#define NEUTRAL_POT_H

#include <string>
#include "numerical_orbital.h"
#include "numerical_vna_lm.h"
using namespace std;

class Neutral_Pot
{
	public: 

	Neutral_Pot();
	~Neutral_Pot();

	void setup_Vna(const int &it);	
	void uniform_Vna(const int &it, const double &dr_in );

	double* vna;
	double* vna_u; // uniform mesh.
	double dr;
	double rcut;
	int nr;

	private:

	void output_1dv(const string &fn, const int &msh, const double* target, const double* r)const;


	//---------------------------------------------------------
	// the following is for the exansion of Neutral Potential
	//---------------------------------------------------------
	public:

	const int& getLmax(void) const {return lmax;}
	const int& getTotal_nchi(void) const { return total_nchi; }
	const int& getNchi(const int l) const { return this->nchi[l]; }


	// initial lmax, nmax, etc.
	void init_proj(Numerical_Orbital &Phi, const double &dr_uniform);


	Numerical_Vna_Lm* VnaLN;// length: total_nchi (only store radial function )

	private:

	int lmax;
	int total_nchi;
	int* nchi;
	int mesh_u;

	void set_vnaphi(double** VnaPhi, Numerical_Orbital &Phi, const double &dr_uniform);
	
};

#endif
