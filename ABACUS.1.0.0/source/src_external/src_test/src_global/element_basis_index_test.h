//==========================================================
// AUTHOR : Peize Lin
// DATE : 2016-06-02
//==========================================================

#ifndef ELEMENT_BASIS_INDEX_TEST_H
#define ELEMENT_BASIS_INDEX_TEST_H

#include "../../..//src_global/element_basis_index.h"

#include <iostream>
using namespace std;

ostream & operator << (ostream & os, const Element_Basis_Index::Range & range)
{
	for( size_t T=0; T!=range.size(); ++T )
	{	
		os<<T<<endl;
		for( size_t L=0; L!=range[T].size(); ++L )
		{
			os<<"\t"<<L<<"\t"<<range[T][L].M<<"\t"<<range[T][L].N<<endl;
		}
	}
	return os;
}

ostream & operator << (ostream & os, const Element_Basis_Index::Index & index)
{
	for( size_t T=0; T!=index.size(); ++T )
	{
		os<<T<<":\t"<<index[T].count_size<<endl;;
		for( size_t L=0; L!=index[T].size(); ++L )
		{
			os<<"\t"<<L<<endl;
			for( size_t M=0; M!=index[T][L].size(); ++M )
			{
				os<<"\t\t"<<M<<endl;
				for( size_t N=0; N!=index[T][L][M].size(); ++N )
				{
					os<<"\t\t\t"<<N<<"\t"<<index[T][L][M][N]<<endl;
				}
			}
		}
	}
	return os;
}
#endif		// ELEMENT_BASIS_INDEX_TEST_H