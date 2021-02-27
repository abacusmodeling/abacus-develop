#ifndef EPSILON0_VASP_H
#define EPSILON0_VASP_H

#include "src_pw/wavefunc.h"
#include "../src_parallel/parallel_global.h"

class Epsilon0_vasp 
{

public:
		
		Epsilon0_vasp();
		~Epsilon0_vasp();
		
		bool epsilon;
		double domega;
		int nomega;
		double eta;
		
		int oband;
		int uband;

		void Cal_psi_nu(int ik);
		
		void cal_epsilon0();
		void Init();
		void Delete();
		void Cal_b(int ik);
		void Cal_psi(int ik);
		void Cal_psi_nabla(int ik);
		void Cal_epsilon0s();
		void Cal_T();
		void Cal_epsilon0();
		void CMatrixMul(int m, int n, int l, complex<double>** A, complex<double>** B, complex<double>** C);
		
		complex<double> **eps0;

private:
		
		bool init_finish;
		complex<double> ***b;
		complex<double> **psi;
		complex<double> ***psi_nabla;
		complex<double> ***psi_nu;
		complex<double> **eps0s;
		complex<double> **T;
};

extern Epsilon0_vasp epsilon0_vasp;

#endif	
