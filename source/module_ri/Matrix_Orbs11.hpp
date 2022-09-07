//=======================
// AUTHOR : Peize Lin
// DATE :   2022-08-17
//=======================

#ifndef MATRIX_ORB11_HPP
#define MATRIX_ORB11_HPP

#include "Matrix_Orbs11.h"

template<typename Tdata>
Tensor<Tdata> Matrix_Orbs11::cal_overlap_matrix( 
	const size_t TA, 
	const size_t TB, 
	const ModuleBase::Vector3<double> &tauA,
	const ModuleBase::Vector3<double> &tauB, 
	const ModuleBase::Element_Basis_Index::IndexLNM &index_A, 
	const ModuleBase::Element_Basis_Index::IndexLNM &index_B,
	const Matrix_Order &matrix_order ) const
{
	Tensor<Tdata> m;
	const size_t sizeA = index_A[TA].count_size;
	const size_t sizeB = index_B[TB].count_size;
	switch(matrix_order)
	{
		case Matrix_Order::AB: m = Tensor<Tdata>({sizeA, sizeB});	break;
		case Matrix_Order::BA: m = Tensor<Tdata>({sizeB, sizeA});	break;
		default:	throw std::invalid_argument(std::string(__FILE__)+" line "+std::to_string(__LINE__));
	}
	
	for( const auto &co3 : center2_orb11_s.at(TA).at(TB) )
	{
		const int LA = co3.first;
		for( const auto &co4 : co3.second )
		{
			const size_t NA = co4.first;
			for( size_t MA=0; MA!=2*LA+1; ++MA )
			{
				for( const auto &co5 : co4.second )
				{
					const int LB = co5.first;
					for( const auto &co6 : co5.second )
					{
						const size_t NB = co6.first;	
						for( size_t MB=0; MB!=2*LB+1; ++MB )
						{
							const Tdata overlap = co6.second.cal_overlap( tauA*GlobalC::ucell.lat0, tauB*GlobalC::ucell.lat0, MA, MB );
							const size_t iA = index_A[TA][LA][NA][MA];
							const size_t iB = index_B[TB][LB][NB][MB];
							switch(matrix_order)
							{
								case Matrix_Order::AB:	m(iA,iB) = overlap;		break;
								case Matrix_Order::BA:	m(iB,iA) = overlap;		break;
								default:	throw std::invalid_argument(std::string(__FILE__)+" line "+std::to_string(__LINE__));
							}
						}
					}
				}
			}
		}
	}
	return m;
}

#endif