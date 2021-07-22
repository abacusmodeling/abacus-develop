#ifndef PDSTEBZ_H
#define PDSTEBZ_H

#include "../src_pw/tools.h"
#include "pdiag_common.h"

void pdstebz(MPI_Comm comm_2D, double *D,double *E,double *eigen,int N);

#endif
