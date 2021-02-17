#ifndef PZGSEPS_H
#define PZGSEPS_H

#include "src_pw/tools.h"
#include "pdiag_common.h"

void pzgseps(
		MPI_Comm comm_2D,
		int n,
		int nb,
		int &egnum,
		complex<double> *A, 
		complex<double> *B,
		complex<double> *Z,
		double *eigen,
		LocalMatrix loc_A,
		char uplo, 
		int &loc_size,
		int &loc_pos);

#endif
