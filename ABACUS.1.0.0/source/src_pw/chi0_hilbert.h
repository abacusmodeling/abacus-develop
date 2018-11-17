//==========================================================
//   AUTHOR : Pengfei Li
//   DATE : 2016-12-14
//==========================================================
#ifndef CHI0_HILBERT_H
#define CHI0_HILBERT_H
#include "../input_conv.h"
#include "wavefunc.h"
#include "../src_parallel/parallel_global.h"
#if defined __FFTW2
#include "../src_parallel/fftw.h"
#elif defined __FFTW3
#include "../src_parallel/fftw3.h"
#endif

class Chi0_hilbert 
{
public:
	Chi0_hilbert();
	~Chi0_hilbert();
	
	bool epsilon;           // calculate epsilon or not
	string kernel_type;      // the type of kernel
	string system;          // bulk or surface
	double eta;             // unit(Ry)
	double domega;          // unit(Ry)
	int nomega;
	int dim;                // the dimension of G 
	int oband;              // the number of "occupied" bands 
	double q_start[3];      // the position of the first q point in direct coordinate
	double direct[3];       // the q direction
	int start_q;            // the serial number of the start qpoint
	int interval_q;         // the interval of the qpoints
	int nq;                 // the total number of qpoints for calculation
	bool out_epsilon;       // output epsilon or not
	bool out_chi;           // output chi or not
	bool out_chi0;          // output chi0 or not
	double fermi_level;     // the change the fermi level(Ry)
	bool coulomb_cutoff;    // turn on or off the Coulomb_cutoff 0/1  (The Z axis must be set to be a vacuum layer)
	bool kmesh_interpolation;       // calculting <i,0|j,R>
	double qcar[100][3];    // the Cartesian position of q points(unit: 2*PI/lat0) 
	int lcao_box[3];        // the scale for searching the existence of the overlap <i,0|j,R>
	
	int dim_para;

        
	void Init();
	void Delete();
	void Chi();
	void Parallel_G();
	void Cal_Psi(int iq, complex<double> **psi_r);
	void Cal_Psi_down(int iq, complex<double> **psi_r);
	void Cal_lcao_psi();
	void Cal_b(int iq, int ik, int iqk, int ispin);
	void Cal_b_lcao( int iq, int ik, int iqk);
	void Cal_Chi0s(int iq);
	void Cal_T();
	void Cal_Chi0();
	complex<double> f(int k, int j);
	void Cal_Chi(int iq);
	void Cal_Rpa(int iq);
	double qg2( int iq, int g0); 
	int Cinv(int n, complex<double>** a);
	void CMatrixMul(int m, int n, int l, complex<double>** A, complex<double>** B, complex<double>** C); 
	int Cal_iq(int ik, int iq, int a, int b, int c);
	
	void Cal_Chi_surface(int iq);
	int parallel_g();
	void chi0_para_g();
	void Cal_kernel(int iq);
	void Cal_kernel_2D(int iq);
	double qg(int iq, int g0);
	double qG(int iq, int g0);
	void Cal_Rpa1();
	void chi_para_g();
	complex<double> Cal_g(int iq);
	void C_inverse(int n, complex<double>** a);
	void plot_chi0(int iq);
	
	Vector3<double> *all_gcar;
	complex<double> **chi0;
	complex<double> **chi;

private:


	bool init_finish;

        complex<double> *b_core;
        int *num_G_core;
        int *num_G_dis;
        complex<double> *b_summary;
        complex<double> *b_order;
        double *G_r_core;
        int *num_Gvector_core;
        int *num_Gvector_dis;
        int *flag;
        double *G_r;
        double *Gvec_core;
        double *Gvec;
        double **cweight;
        complex<double> **psi_r1;
        complex<double> **psi_r2;
        complex<double> ***b;
        complex<double> **chi0s;
        complex<double> **chi0_gg;
        complex<double> **T;
		complex<double> **kernel;
        complex<double> **rpa;

        int *flag1;
        double **para_g;
        complex<double> **chi0_para;
        complex<double> **rpa1;
        complex<double> **chi_para;
 
        Vector3<int> ***R; 
        Vector3<double> ***Rcar; 
        int ** Rmax;
        int NR;
        complex<double> ****overlap;
};

extern Chi0_hilbert chi0_hilbert;

#endif
