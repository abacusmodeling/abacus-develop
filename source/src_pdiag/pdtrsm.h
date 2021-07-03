#ifndef PDTRSM_H
#define PDTRSM_H

#include "../src_pw/tools.h"
#include "pdiag_common.h"

void pdtrsm(char isuplo,int b_n,MPI_Comm comm_2D,int NB,int N_A,
	    double *A,double *B,LocalMatrix loc_A,LocalMatrix loc_B);

#endif
