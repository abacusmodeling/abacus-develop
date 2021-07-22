#ifndef PZST2G_H
#define PZST2G_H

#include "../src_pw/tools.h"
#include "pdiag_common.h"

void pzst2g(MPI_Comm comm_2D,int NB,int N_A,complex<double> *A,complex<double> *B,
	     LocalMatrix loc_A,int loc_size,char isuplo);

#endif
