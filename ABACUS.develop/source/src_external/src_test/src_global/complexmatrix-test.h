//==========================================================
// AUTHOR : Peize Lin
// DATE : 2016-09-08
//==========================================================

#ifndef COMPLEXMATRIX_TEST_H
#define COMPLEXMATRIX_TEST_H

#include "src_global/complexmatrix.h"
#include<limits>

static std::ostream & operator<<( std::ostream & os, const ComplexMatrix & m )
{
	for( int ir=0; ir!=m.nr; ++ir )
	{
		for( int ic=0; ic!=m.nc; ++ic )
			os<<m(ir,ic)<<"\t";
		os<<std::endl;
	}	
	return os;
}

#endif