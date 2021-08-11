#ifndef PZHEG2ST_H
#define PZHEG2ST_H

#include "../src_pw/tools.h"
#include "pdiag_common.h"

void pzheg2st(
	char isuplo,
	MPI_Comm comm_2D,
	std::complex<double> *A, 
	LocalMatrix loc_A,
	std::complex<double> *B,
	LocalMatrix loc_B,
	int NB,
	int N_A);

#endif
