#ifndef CHI0_STANDARD_H
#define CHI0_STANDARD_H

#include "wavefunc.h"
#include "../src_parallel/parallel_global.h"
#if defined __FFTW2
#include "../src_parallel/fftw.h"
#elif defined __FFTW3
#include "../src_parallel/fftw3.h"
#endif

class Chi0_standard
{
public:

      Chi0_standard();
	  ~Chi0_standard();
	  
	  bool epsilon;
	  string system;
	  double eta;
	  double domega;
	  int nomega;
	  int dim;
	  int oband;
	  double q_start[3];
	  double direct[3];
	  int start_q;
	  int interval_q;
	  int nq;
	  bool out_epsilon;
	  
	  void Chi();
	  void Parallel_G();
	  void Init();
	  void Delete();
	  void Cal_Psi(int iq, complex<double> **psi_r);
	  void Cal_b(int iq, int ik, int iqk);
	  void Cal_weight(int iq, int ik, double omega);
	  void Cal_last();
	  void Cal_first();
	  void Cal_chi0(int iq, double omega);
	  void Cal_rpa(int iq);
	  void Cal_chi();
	  double qg2( int iq, int g0);
	  int Cinv(int n, complex<double>** a);
	  void CMatrixMul(int m, int n, int l, complex<double>** A, complex<double>** B, complex<double>** C);
	  int Cal_iq(int ik, int iq, int a, int b, int c);
	  int parallel_g();
	  void chi0_para_g();
	  
	  Vector3<double> *all_gcar;	  
	  complex<double> **chi0;
	  complex<double> **chi;
	  
	  
private:
	  
	  bool init_finish;
	  complex<double> **psi_r1;
	  complex<double> **psi_r2;
	  complex<double> **b;
	  complex<double> **A;
	  complex<double> **B;
	  complex<double> *weight;
	  complex<double> **rpa;
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
	  int *flag1;
	  double **para_g;
	  int dim_para;
	  complex<double> **chi0_para;
	  	
};

extern Chi0_standard chi0_standard;

#endif