#ifndef EXX_ABFS_INL_H
#define EXX_ABFS_INL_H

#include "exx_abfs.h"

template<typename T1, typename T2> 
void Exx_Abfs::minus_matrixes(
	std::map<T1,T2> &A,
	const std::map<T1,T2> &B) const
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
	std::vector<T> &A,
	const std::vector<T> &B) const
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
std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,matrix>>>> &Exx_Abfs::minus_matrixes(
	std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,matrix>>>> &mAs, 
	const std::map<size_t,std::map<size_t,std::map<size_t,std::map<size_t,matrix>>>> & mBs) const
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
					std::cout<<i1<<i2<<i3<<i4<<std::endl;
					std::cout<<mAs[i1][i2][i3][i4]<<std::endl;
					std::cout<<m4.second<<std::endl;
					
					mAs[i1][i2][i3][i4] -= m4.second;
					
					// Peize Lin test
					std::cout<<mAs[i1][i2][i3][i4]<<std::endl;
				}
			}
		}	
	}
	return mAs;
}
*/

#endif // EXX_ABFS_INL_H