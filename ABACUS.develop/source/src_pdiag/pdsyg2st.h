#ifndef PDSYG2ST_H
#define PDSYG2ST_H

#include "src_pw/tools.h"
#include "pdiag_common.h"

void pdsyg2st(
	char isuplo,
	MPI_Comm comm_2D,
	double *A, 
	LocalMatrix loc_A,
	double *B,
	LocalMatrix loc_B,
	int NB,
	int N_A);

#endif
