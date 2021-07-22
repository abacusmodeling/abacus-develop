//==========================================================
// AUTHOR : Peize Lin
// DATE : 2016-06-02
//==========================================================

#pragma once


#include <iostream>
using namespace std;

static ostream & operator << (ostream & os, const Element_Basis_Index::Range & range)
{
	for( size_t T=0; T!=range.size(); ++T )
	{	
		os<<T<<endl;
		for( size_t L=0; L!=range[T].size(); ++L )
		{
			os<<"\t"<<L<<"\t"<<range[T][L].N<<"\t"<<range[T][L].M<<endl;
		}
	}
	return os;
}

static ostream & operator << (ostream & os, const Element_Basis_Index::IndexLNM & index)
{
	for( size_t T=0; T!=index.size(); ++T )
	{
		os<<T<<":\t"<<index[T].count_size<<endl;;
		for( size_t L=0; L!=index[T].size(); ++L )
		{
			os<<"\t"<<L<<endl;
			for( size_t N=0; N!=index[T][L].size(); ++N )
			{
				os<<"\t\t"<<N<<endl;
				for( size_t M=0; M!=index[T][L][N].size(); ++M )
				{
					os<<"\t\t\t"<<M<<"\t"<<index[T][L][N][M]<<endl;
				}
			}
		}
	}
	return os;
}