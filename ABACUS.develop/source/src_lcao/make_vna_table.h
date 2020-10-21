#ifndef MAKE_VNA_TABLE_H
#define MAKE_VNA_TABLE_H

#include "../src_pw/tools.h"
#include "neutral_pot.h"
#include "numerical_orbital_lm.h"

class Make_Vna_Table 
{
	public:
		
	Make_Vna_Table();
	~Make_Vna_Table();

    void allocate (
        const int &ntype,
        const int &lmax_in,
        const int &kmesh_in,
        const double &Rmax_in,
        const double &dR_in,
        const double &dk_in);

	
    static double dr;
    int Rmesh;

	static int get_rmesh( const double &R1, const double &R2);

	void init_Table_Chi(void);
	void Destroy_Table_Chi(void);


    // Five dimension:
    // (1) 0: normal (S(R)) ; 1: derivative( dS/dR )
    // (2) pairs type number.
    // (3) pairs chi.
    // (4) Max angular momentum: L.
    // (5) Distance between atoms: R.
	double***** Table_VnaR;
    
	
	void init_Vna_Tpair(void);
	void init_Vna_Opair(void);

	
	// number of pairs
    int Vna_nTpairs;

	// index for atom type pairs
    IntArray Vna_Tpair;

	IntArray Vna_Opair;
	
    IntArray Vna_L2plus1;

	bool destroy_Vna;

	private:


	// phi is the atomic orbital, chi is the projector of Vna
	void cal_Vna_PhiChi_R(
		const int &l,
		const Numerical_Orbital_Lm &n1,
		const Neutral_Pot &n2,
		const int &rmesh,
		double* rs,
		double* drs);


	//variables
    int ntype;
    int lmax;
    double Rmax;
    double dk;
    int nlm;
    int kmesh;
    double *kpoint;
    double *r;
    double *rab;
    double *kab;


};

#endif	
