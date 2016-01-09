//=========================================================
//AUTHOR : mohan
//DATE : 2009-09-16
//=========================================================
#ifndef GINT_GAMMA_H
#define GINT_GAMMA_H

#include "../src_pw/tools.h"
#include "grid_base_beta.h"
#include "grid_technique.h"
#include "lcao_matrix.h"
//=========================================================
//CLASS  Grid_Integral_Beta
//Note : Integral On 3D Grids, different from Grid_Integral
//Feature : Matrix Elements Of Local Potential For 
//Numerical Orbitals
//=========================================================

class Gint_Gamma : public Grid_Base_Beta
{
private:

	double* transformer;
	double psiv1;
	double psiv2;
	double* ylm1;
	double* ylm2;

	int grid_index;
	int max_size;
	
	// these parameters are for interpolation.
	// we store these parameters at first to speed
	// up the calculation.
	double *x0;
	double *x1;
	double *x2;
	double *x3;
	double* x12;
	double* x03;
	int *iq;

	void save_atoms_on_grid(const Grid_Technique &gt);
	
	
	// for calculation of < phi_i | Vna | phi_j >
	// on normal real space FFT grid.
	void gamma_vna(void);


	// for calculation of < phi_i | Vna | phi_j >
	// on dense real space FFT grid.
	void gamma_vna_d(const Grid_Technique &gt, const char &matrix_type);// d stands for 'dense'

	// for calculation of < phi_i | Vlocal | phi_j >
	void gamma_vlocal(void);  


	// for calculation of pVp under B field.
	void gamma_vlocal_B(void);  


	// for calculation of charege 
	double gamma_charge(void);

	
	// for calculatin of charge under B field
	double gamma_charge_B(void); //mohan add 2012-04-19


	// for calculation of Mulliken charge.
	void gamma_mulliken(double** mulliken);


	// for calculation of envelope functions.
	void gamma_envelope(const double* wfc, double* rho);// mohan add 2011-07-01


	// for calculation of < beta_i | phi_j> under magnetic B field.
	void grid_integration_vnl(complex<double> **Tab);  


	// for calculatin of < dphi_i | Vlocal | phi_j > for foce calculation
	// on regular FFT real space grid.
	void gamma_force(void);


	// for calculation of < dphi_i | Vlocal | phi_j> for force calculation
	// on dense FFT grid.
	void gamma_force_vna(const Grid_Technique &gt, LCAO_Matrix &lm);


	// for claculation of terms under B field,
	void gamma_S_T_AP(char type, const Grid_Technique &gt);// sun zhiyuan add 2012-01-03//


    //To be used in cal_S_T_AP(), zhiyuan add 2011-01-06//
    complex<double> Add_S(complex<double> phase, double psi1, double psi2);

    complex<double> Add_T(complex<double> phase, double A1[3], double A2[3], 
			double dphi1x,double dphi1y,double dphi1z,
			double dphi2x,double dphi2y,double dphi2z, 
			double psi1, double psi2);

    complex<double> Add_AP(complex<double> phase, double Aldr3[3],double A2[3],
			double dphi2x,double dphi2y,double dphi2z,
			double psi1, double psi2);

    complex<double> Add_AA(complex<double> phase, double AAdr,
			double psi1, double psi2);


public:

	Gint_Gamma();
	~Gint_Gamma();

	double cal_rho(void);
	void cal_mulliken(double** mulliken);
	void cal_env(const double* wfc, double* rho);
	void cal_vna( const double* vlocal_in);
	void cal_vna_d( const Grid_Technique &gt, const double* vlocal_in, const char &matrix_type);// d stands for 'dense'
	void cal_vlocal( const double* vlocal_in);
	void cal_force( const double* vlocal_in);
	void cal_force_vna( const double* vlocal_in, const Grid_Technique &gt, LCAO_Matrix &lm );
	void cal_vnl_B(complex<double> **Tab); // mohan add 2011-04-08
	
	
	
	// needed while for gauge-invarient B-field.
	void cal_S_T_AP(char type, const Grid_Technique &gt);//for calculating S,T,and A*P. sun zhiyuan add 2012-01-03


};

#endif
