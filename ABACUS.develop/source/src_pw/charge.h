//==========================================================
// AUTHOR : Lixin He, mohan , fangwei
// DATE : 2008-11-10
//==========================================================
#ifndef CHARGE_H
#define CHARGE_H

#include "tools.h"
#include "../src_parallel/parallel_global.h"

//==========================================================
// Electron Charge Density 
//==========================================================

class Charge
{

public:

    Charge();
    ~Charge();

	// mohan update 2021-02-20
	void allocate(const int &nspin_in, const int &nrxx_in, const int &ngmc_in);

//==========================================================
// MEMBER VARIABLES :
// NAME : rho (nspin,ncxyz), the charge density in real space
// NAME : rho_save (nspin,ncxyz), another charge density in
// NAME : rho_ion(extra_pot, nspin, nxyz), mohan add 2011-03-22
// rho_ion save the previous ion step charge density.
// for charge extrapolation purpose.
// real space
// NAME : rho_core [nrxx], the core charge in real space
// NAME : rhog_core [ngm], the core charge in reciprocal space
//==========================================================

    double** rho;
    double** rho_save;

    complex<double>** rhog;
    complex<double>** rhog_save;

    double *rho_core;
	complex<double> *rhog_core; // mohan add 2009-12-15

    int out_charge; // output charge if out_charge > 0, and output every "out_charge" elec step.

    double *start_mag_type;
    double *start_mag_atom;


    void atomic_rho(const int spin_number_need, double **rho_in)const;

    void set_rho_core(const ComplexMatrix &structure_factor);

    void write_rho(const int &is, const int &iter, const string &fn, 
		const int &precision = 11, const bool for_plot = false);//mohan add 2007-10-17

    void write_rho_dipole(const int &is, const int &iter, const string &fn, 
		const int &precision = 11, const bool for_plot = false);//fuxiang add 2017-3-15    

    bool read_rho(const int &is, const string &fn);//mohan add 2007-10-17

    void sum_band(void);

    void renormalize_rho(void);

    void save_rho_before_sum_band(void);

    void non_linear_core_correction// drhoc
    (
        const bool &numeric,
        const int mesh,
        const double *r,
        const double *rab,
        const double *rhoc,
        double *rhocg
    );

	double check_ne(const double *rho_in) const;

    void init_final_scf(); //LiuXh add 20180619

	private:

	// mohan add 2021-02-20
	int nrxx; // number of r vectors in this processor
	int ngmc; // number of g vectors in this processor
	int nspin; // number of spins

    void sum_band_k();

    void rho_mpi(void);

    double sum_rho(void) const;

	bool allocate_rho;

	bool allocate_rho_final_scf; //LiuXh add 20180606

};

#endif //charge

