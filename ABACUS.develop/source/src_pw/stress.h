#ifndef STRESS_H
#define STRESS_H

#include "tools.h"

//-------------------------------------------------------------------
// mohan reconstruction note: 2021-02-07
// the stress code needs reconstructions (by Daye Zheng)
// 1) add explanations for each function, each variable, for 
// main procedures, make the code readable
// 2) divide the stress class into several files, each file
// deals with only one part of the stress, it is convenient for 
// the next-step reconstruction, for example, we want to make 
// pw as an external variable instead of a global variable
// 3) for PW and LCAO, keeps only one copy of the code, for example
// the ewald term needs only one copy, it will reduce the 
// time to maintain both copies of stress codes.
// 4) remain openning interfaces for others to contribute stress
// codes, for example, molecular dynamics will have an ionic stress
// term, +U? exx? may have other stress terms.
// 5) delete useless comments and tests, if you have a useless code,
// please explicitly explain why you want to keep the test
// 6) format should be beautiful! code should be readable like a 
// note (let readers be comfortable)
//-------------------------------------------------------------------

//----------------------------------------------------------------
// compute the stress terms in terms of the plane wave basis set
// the stress terms include: 
// 1) the stress from the electron kinetic energy 
// 2) the stress from the local pseudopotentials
// 3) the stress from the non-local pseudopotentials
// 4) the stress from the Hartree term
// 5) the stress from the non-linear core correction (if any)
// 6) the strees from the exchange-correlation functional term
// 7) the stress from the ewald term (ion-ion intraction under 
//		periodic boundary conditions). 
// 8) the stress from ionic contributions (for molecular dynamics)
//----------------------------------------------------------------
using namespace std;

class Stress
{
	public: 

	Stress(){};
	~Stress(){};

	void cal_stress();

	void dvloc_of_g (const int& msh,
			const double* rab,
			const double* r,
			const double* vloc_at,
			const double& zp,
			double*  dvloc);

	void dvloc_coul (const double& zp, double* dvloc);

	void deriv_drhoc (
			const bool &numeric,
			const int mesh,
			const double *r,
			const double *rab,
			const double *rhoc,
			double *drhocg);

	void print_stress(const string &name, double f[][3], const bool screen, bool ry)const;

	void printstress_total (bool ry);
	
	// total stress
	double sigmatot[3][3];

	private:

	void stres_knl();

	void stres_har();

	void stres_ewa();

	void stres_loc();

	void stres_cc();

	void stres_gradcorr();

	void stres_nl();

	void dylmr2 (
			const int nylm,
			const int ngy,
			Vector3<double> *gk,
			matrix &dylm,
			const int ipol);

	double Polynomial_Interpolation_nl(
			const realArray &table,
			const int &dim1,
			const int &dim2,
			const double &table_interval,
			const double &x);

	void get_dvnl2(
			ComplexMatrix &vkb,
			const int ik);

	void get_dvnl1(
			ComplexMatrix &vkb,
			const int ik,
			const int ipol);

	// stress for exchange-correlation functional term
    double sigmaxc[3][3];
	// stress for electron kinetic energy term
	double sigmakin[3][3];
	// stress for hartree term
	double sigmahar[3][3];
	// stress for local pseudopotential term
	double sigmaloc[3][3];
	// stress for nonlocal pseudopotential term
	double sigmanlc[3][3];
	// stress for Ewald term (ion-ion interaction)
	double sigmaewa[3][3];
	// stress for non-linear core correction term
	double sigmaxcc[3][3];
}; 







#endif
