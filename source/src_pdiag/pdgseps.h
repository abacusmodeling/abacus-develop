#ifndef PDGSEPS_H
#define PDGSEPS_H

#include "../src_pw/tools.h"
#include "pdiag_common.h"

void pdgseps(
		MPI_Comm comm_2D,
		int n,
		int nb,
		double *A, 
		double *B,
		double *Z,
		double *eigen,
		LocalMatrix loc_A,
		char uplo, 
		int &loc_size,
		int &loc_pos);

#endif
