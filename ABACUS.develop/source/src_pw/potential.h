#ifndef POTENTIAL_H
#define POTENTIAL_H
#include "tools.h"

class potential
{
	public:

	friend class Hamilt_PW; 
	friend class Electrons;
	friend class ELEC_scf;

    // constructor and deconstructor
    potential();
    ~potential();

    //==========================================================
    // start_pot : "atomic" or "file"
	// extra_pot : extrapolation methods for potential
    // vr(nspin,ncxyz) : Hartree + xc potentials in real space
    // vrs(nspin,ncxyz) : total potential in real space
    // vnew(nspin,ncxyz) : V_out - V_in, needed in scf
	// vltot: the local potential in real space
	// out_potential: options to print out potentials 
    //==========================================================
    string start_pot;
	string extra_pot;
    matrix vr;
    matrix vrs;
    matrix vnew;
    double *vrs1;
    double *vltot;
	int out_potential; // mohan add 2011-02-28

    void allocate(const int nrxx);

	void init_pot(const int &istep);

    void v_of_rho(double** rho_in, matrix &v);

    void set_vrs(void);


    void newd(void);

	public:

	// mohan add 2011-02-28
	// here vh is complex because the array is got after complex FFT.
	void write_potential(const int &is, const int &iter, const string &fn,
		const matrix &v, const int &precision, const int &hartree = 0)const;

    void write_elecstat_pot(const string &fn, const string &fn_ave);
	
	private:

    void set_local(double *vl_pseudo)const;



	private:

	// TDDFT related, fuxiang add
    double *vext;

    double *vextold;

    void set_vrs_tddft(const int istep);

};

#endif	// POTENTIAL_H
