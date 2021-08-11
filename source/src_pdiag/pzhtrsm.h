#ifndef PZHTRSM_H
#define PZHTRSM_H

#include "../src_pw/tools.h"
#include "pdiag_common.h"

void pzhtrsm(char isuplo,int b_n,MPI_Comm comm_2D,int NB,int N_A,
	    std::complex<double> *A,std::complex<double> *B,LocalMatrix loc_A,LocalMatrix loc_B);

#endif
