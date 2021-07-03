#ifndef PDT2S_H
#define PDT2S_H

#include "../src_pw/tools.h"
#include "pdiag_common.h"

void pdt2s(MPI_Comm comm2D,int N_A,int NB,double *A,double *X,
	     LocalMatrix loc_A,int loc_size,double *norm,char uplo);

#endif
