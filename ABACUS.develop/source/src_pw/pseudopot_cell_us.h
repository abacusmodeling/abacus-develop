#ifndef pseudopot_cell_US_H
#define pseudopot_cell_US_H

#include "tools.h"
#include "pseudopot_cell_vnl.h"

class pseudopot_cell_us: public pseudopot_cell_vnl
{
public:
	// PUBLIC :: nlx, lpx, lpl, ap, aainit, indv, nhtol, nhtolm, nkb, nkbus,
	//   vkb, dvan, deeq, qq, nhtoj, beta, becsum, deallocate_uspp
	//  PUBLIC :: qq_so, dvan_so, deeq_nc

//	int nlx;	// = (lmaxx+1)**2,
	// maximum number of combined angular momentum
//	int	mx;		// = 2*lqmax-1, maximum magnetic angular momentum of Q
//	matrix lpx;	// (nlx,nlx); for each pair of combined momenta lm(1),lm(2):
	// maximum combined angular momentum LM
//	realArray lpl;	//(nlx,nlx,mx)  list of combined angular momenta  LM

//	realArray ap;	//(lqmax*lqmax,nlx,nlx)
	// Clebsch-Gordan coefficients for spherical harmonics

//	int nkbus;        // as above, for US-PP only

	realArray becsum;	//(:,:,:), \sum_i f(i) <psi(i)|beta_l><beta_m|psi(i)>
	realArray beta;		//(:,:,:), beta functions for CP (without struct.factor) //not sure!

	// from spin_orb
//	bool lspinorb,
//	domag;			// if .TRUE. this is a spin-robit calculation
	// from us
//	ComplexMatrix rot_ylm;	//(2*lmaxx+1,2*lmaxx+1), transform real
	// spherical harmonics into complex ones
//	ComplexArray fcoef;		// (:,:,:,:,:), function needed to
	// account for spinors.

	// from uspp
//	ComplexArray dvan_so,	// (:,:,:,:), D_{nm}
//	qq_so;				// (:,:,:,:), Q_{nm}

	//MODULE us
	//. These parameters are needed with the US pseudopotentials

//	::real dq = 0.01;   // space between points in the pseudopotential tab.
	// moved to constants.h

//	realArray qrad;		//(:,:,:,:),radial FT of Q functions

	// END MODULE us

	pseudopot_cell_us();
	~pseudopot_cell_us();
	/*
		// member functions:
	//	void init_us_1();  // generate  qrad(:,:,:,:) ?
		void init_us_1(PW_Basis pw);	// generate tab(:,:,:), indv, nhtol, nhtolm
		void aainit(int);     //Clebsch-Gordan coefficients "ap"
		::real compute_ap(int l,int li,int lj,int llx,matrix ylm,matrix mly);
		void gen_rndm_r(int llx, Vector3 <::real> *r, ::real *rr);
	//	void gen_rndm_r(int llx,Vector3 *r,::real *rr);
		void qvan2 (int, int ih, int jh, int nt, ::real *gg,
			complex<::real> qgm, matrix ylmk0);
	//	void newd();  //computes the fourier transform of the Q function, qq

		void printus(ofstream &ofs);

	//	template <class T>
	//	void copyarr(T *, T *, int len);
	//	void copyarr(::real *, ::real *, int );

	*/
};

#endif // PSEUDOPOT_CELL_US_H
