#ifndef PDST2G_H
#define PDST2G_H

#include "../src_pw/tools.h"
#include "pdiag_common.h"

void pdst2g(MPI_Comm comm_2D,int NB,int N_A,double *A,double *B,
	     LocalMatrix loc_A,int loc_size,char isuplo);

#endif
