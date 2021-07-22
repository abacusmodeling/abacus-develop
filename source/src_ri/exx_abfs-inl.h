#ifndef EXX_ABFS_INL_H
#define EXX_ABFS_INL_H

#include "exx_abfs.h"

template<typename T1, typename T2> 
void Exx_Abfs::minus_matrixes(
	map<T1,T2> &A,
	const map<T1,T2> &B) const
{
	for( auto & a : A )
	{
		const T1 i = a.first;
		minus_matrixes( a.second, B.at(i) );
	}
	return A;
}

template<typename T>
void Exx_Abfs::minus_matrixes(
	vector<T> &A,
	const vector<T> &B) const
{
	for( size_t i=0; i!=A.size(); ++i )
		minus_matrixes( A[i], B[i] );
	return A;
}

template<>
void Exx_Abfs::minus_matrixes(
	matrix &A,
	const matrix &B) const
{
	A -= B;
}


/*
map<size_t,map<size_t,map<size_t,map<size_t,matrix>>>> &Exx_Abfs::minus_matrixes(
	map<size_t,map<size_t,map<size_t,map<size_t,matrix>>>> &mAs, 
	const map<size_t,map<size_t,map<size_t,map<size_t,matrix>>>> & mBs) const
{
	for( const auto m1 : mBs )
	{
		const size_t i1 = m1.first;
		for( const auto m2 : m1.second )
		{
			const size_t i2 = m2.first;
			for( const auto m3 : m2.second )
			{
				const size_t i3 = m3.first;
				for( const auto m4 : m3.second )
				{
					const size_t i4 = m4.first;
					
					// Peize Lin test
					cout<<i1<<i2<<i3<<i4<<endl;
					cout<<mAs[i1][i2][i3][i4]<<endl;
					cout<<m4.second<<endl;
					
					mAs[i1][i2][i3][i4] -= m4.second;
					
					// Peize Lin test
					cout<<mAs[i1][i2][i3][i4]<<endl;
				}
			}
		}	
	}
	return mAs;
}
*/

#endif // EXX_ABFS_INL_H