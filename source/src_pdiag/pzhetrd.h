#ifndef PZHETRD_H
#define PZHETRD_H

#include "../src_pw/tools.h"
#include "pdiag_common.h"
#include "pzhetrd.h"

void pzhetrd(MPI_Comm comm2D, LocalMatrix loc_A,int N,int NB,
	    complex<double> *a,double *diag,double *off_diag,complex<double>  *norm,char uplo);

#endif
