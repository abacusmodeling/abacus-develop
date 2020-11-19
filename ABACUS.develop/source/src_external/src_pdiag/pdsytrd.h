#ifndef PDSYTRD_H
#define PDSYTRD_H

#include "../../src_pw/tools.h"
#include "pdiag_common.h"
#include "pdsytrd.h"

void pdsytrd(MPI_Comm comm2D, LocalMatrix loc_A,int N,int NB,
	    double *a, double *diag, double *off_diag, double *norm,char uplo);

#endif
