#ifndef POTENTIAL_H
#define POTENTIAL_H
#include "tools.h"

class potential
{
	public:

    // constructor and deconstructor
    potential();
    ~potential();

    //==================================================================
    // start_pot : "atomic" or "file"
    // vr(nspin,ncxyz) : The Hartree + xc potential in real space
    // vrs(nspin,ncxyz): The total potential in real space(smooth grid)
    // vnew(nspin,ncxyz) : V_out - V_in , needed in scf
    //==================================================================
    string start_pot;
	string extra_pot; //mohan add 2011-03-13, xiaohui modify 2015-02-01
    matrix vr;
    matrix vrs;
    matrix vnew;
    double *vrs1;	// mohan add 2007-11-12
    double *vext;		// fuxiang add 2017-05
    double *vextold;		//fuxiang add 2018-01-15


    // member functions
    void init(const int nrxx);
    
	// if delta_vh==true, vl_pseudo not calculated,
	// V(rho-rho_atom) is calculated,
	// V(xc) is calculated.

	// if vna==true, vl_pseudo is calculated,
	// V(rho_atom) is calclated,
	// V(xc) is not calculated.
	void init_pot(const int &istep, const bool delta_vh=false, const bool vna=false);
    void newd(void);

    void v_of_rho( double** rho_in, double &ehart, double &etxc,
        double &vtxc, matrix &v, const bool delta_vh=false, const bool vna=false);

    // mix input and output potentials.
    //	void mix_potential (int ndim, double *vout, double *vin, double alphamix,
    //		double &dr2,double tr2,int iter,int n_iter,string filename,bool &conv);

    // vrs = vr + vltotK
    // vr : Hartree potential + V_xc potential .
    // vltot : From pseudopotential .
    void set_vrs(const bool doublegrid);

	// I guess this is done by Fuxiang He, -- mohan 2021-02-01
    void set_vrs_tddft(const bool doublegrid, const int istep);

    void print_pot(ofstream &ofs)const;
    double vr_ave(const int, const int, const double *) ;

    void v_xc(double** rho_in, double &etxc, double &vtxc, matrix &v);

    double *vltot; 	// = new [ncxyz],the local potential in real space

	int out_potential; // mohan add 2011-02-28

	// mohan add 2011-02-28
	// here vh is complex because the array is got after complex FFT.
	void write_potential(const int &is, const int &iter, const string &fn,
		const matrix &v, const int &precision, const int &hartree = 0)const;
    void write_elecstat_pot(const string &fn, const string &fn_ave);
	
	private:

    // use fft to set vltot.
    // do once in demo.cpp
    void set_local(double *vl_pseudo)const;
    void v_h( int nspin, double &ehart, matrix &v, double** rho);

    int test;
};

#endif	// POTENTIAL_H
