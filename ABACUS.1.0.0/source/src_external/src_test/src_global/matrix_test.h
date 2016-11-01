//==========================================================
// AUTHOR : Peize Lin
// DATE : 2016-09-08
//==========================================================

#ifndef MATRIX_TEST_H
#define MATRIX_TEST_H

#include "../../../src_global/matrix.h"
#include<limits>

static std::ostream & operator<<( std::ostream & os, const matrix & m )
{
	for( int ir=0; ir!=m.nr; ++ir )
	{
		for( int ic=0; ic!=m.nc; ++ic )
			os<<m(ir,ic)<<"\t";
		os<<std::endl;
	}	
	return os;
}

static double sum( const matrix & m )
{
	double value = 0.0;
	for( int ir=0; ir!=m.nr; ++ir )
		for( int ic=0; ic!=m.nc; ++ic )
			value += m(ir,ic);
	return value;
}

static double average( const matrix & m )
{
	return sum(m) / (m.nr*m.nc);
}

#endif	// MATRIX_TEST_H