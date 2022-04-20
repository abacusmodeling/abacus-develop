#ifndef PARALLEL_COMMON_H
#define PARALLEL_COMMON_H

#ifdef __MPI
#include "mpi.h"
#endif
#include <complex>
#include <string>

namespace Parallel_Common
{
	//(1) bcast array
	void bcast_complex_double( std::complex<double> *object, const int n);
	void bcast_string(std::string *object,const int n);
	void bcast_double(double *object,const int n);
	void bcast_int(int *object,const int n);
	void bcast_char(char *object,const int n);
	
	//(2) bcast single
	void bcast_complex_double( std::complex<double> &object);
	void bcast_string(std::string &object);
	void bcast_double(double &object);
	void bcast_int(int &object);
	void bcast_bool(bool &object);

}

#endif
