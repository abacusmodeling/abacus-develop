#include "exx_opt_orb.h"
#include "src_pw/global.h"
#include "src_ri/exx_abfs-jle.h"

void Exx_Opt_Orb::print_matrix(
	const string &file_name,
	const vector<matrix> &matrix_Q, 
	const vector<vector<matrix>> &matrix_S,
	const matrix &matrix_V,
	const size_t TA, const size_t IA, const size_t TB, const size_t IB,
	const Element_Basis_Index::Range &range_jles, 
	const Element_Basis_Index::IndexLNM &index_jles) const
{
	auto print_header = [&]( ofstream &ofs )
	{
		ofs << ucell.lat0 << endl;

		ofs << ucell.latvec.e11 << " " << ucell.latvec.e12 << " " << ucell.latvec.e13 << endl;
		ofs << ucell.latvec.e21 << " " << ucell.latvec.e22 << " " << ucell.latvec.e23 << endl;
		ofs << ucell.latvec.e31 << " " << ucell.latvec.e32 << " " << ucell.latvec.e33 << endl;
		
		if( TA==TB )
		{
			ofs << 1 << " ntype" << endl;
			ofs << ucell.atoms[TA].label << " label" << endl;
			if( IA==IB )
			{
				ofs << 1 << " na" << endl;
				ofs << ucell.atoms[TA].tau[IA].x << " " 
					<< ucell.atoms[TA].tau[IA].y << " " 
					<< ucell.atoms[TA].tau[IA].z << endl;
			}
			else
			{
				ofs << 2 << " na" << endl;
				ofs << ucell.atoms[TA].tau[IA].x << " " 
					<< ucell.atoms[TA].tau[IA].y << " "
					<< ucell.atoms[TA].tau[IA].z << endl;
				ofs << ucell.atoms[TB].tau[IB].x << " " 
					<< ucell.atoms[TB].tau[IB].y << " " 
					<< ucell.atoms[TB].tau[IB].z << endl;
			}
		}
		else
		{
			ofs << 2 << " ntype" << endl;
			ofs << ucell.atoms[TA].label << " label" << endl;
			ofs << 1 << " na" << endl;
			ofs << ucell.atoms[TA].tau[IA].x << " " 
				<< ucell.atoms[TA].tau[IA].y << " " 
				<< ucell.atoms[TA].tau[IA].z << endl;
			ofs << ucell.atoms[TB].label << " label" << endl;
			ofs << 1 << " na" << endl;
			ofs << ucell.atoms[TB].tau[IB].x << " " 
				<< ucell.atoms[TB].tau[IB].y << " " 
				<< ucell.atoms[TB].tau[IB].z << endl;
		}
		
		// ecutwfc_jlq determine the jlq corresponding to plane wave calculation.
		ofs << Exx_Abfs::Jle::Ecut_exx << " ecutwfc" << endl; // mohan add 2009-09-08

		// this parameter determine the total number of jlq.
		ofs << Exx_Abfs::Jle::Ecut_exx << " ecutwfc_jlq" << endl;//mohan modify 2009-09-08

		if(TA==TB)
			ofs << ORB.Phi[TA].getRcut() << " rcut_Jlq" << endl;
		else
			ofs << ORB.Phi[TA].getRcut() << " " << ORB.Phi[TB].getRcut() << " rcut_Jlq" << endl;

		// mohan add 'smooth' and 'sigma' 2009-08-28
		ofs << 0 << " smooth" << endl;
		ofs << 0 << " sigma" << endl;

		ofs << Exx_Abfs::Jle::tolerence << " tolerence" << endl;

		ofs << Exx_Abfs::Jle::Lmax << " lmax" << endl;

		ofs << kv.nkstot << " nks" << endl;
		assert( matrix_V.nr == matrix_V.nc );
		ofs	<< matrix_V.nr << " nbands" << endl;
		
		auto cal_sum_M = [&range_jles](size_t T) -> size_t
		{
			size_t sum_M = 0;
			for( size_t L = 0; L!=range_jles[T].size(); ++L )
				sum_M += range_jles[T][L].M;
			return sum_M;
		};
		const size_t nwfc = (TA==TB && IA==IB) ? cal_sum_M(TA) : cal_sum_M(TA)+cal_sum_M(TB);
		ofs	<< nwfc << " nwfc" << endl;
		
		const size_t ecut_numberA = static_cast<size_t>( sqrt( Exx_Abfs::Jle::Ecut_exx ) * ORB.Phi[TA].getRcut() / PI ); // Rydberg Unit
		const size_t ecut_numberB = static_cast<size_t>( sqrt( Exx_Abfs::Jle::Ecut_exx ) * ORB.Phi[TB].getRcut() / PI ); // Rydberg Unit
		if(TA==TB)
			ofs	<< ecut_numberA << " ne" << endl;
		else
			ofs	<< ecut_numberA << " " << ecut_numberB << " ne" << endl;
		
		ofs << "<WEIGHT_OF_KPOINTS>" << endl;
		for( int ik=0; ik!=kv.nkstot; ++ik )		
		{
			ofs << kv.kvec_c[ik].x << " " << kv.kvec_c[ik].y << " " << kv.kvec_c[ik].z;
			ofs << " " << kv.wk[ik] * 0.5 << endl;
		}
		ofs << "</WEIGHT_OF_KPOINTS>" << endl;

		ofs << endl;
	};
	
	
	auto print_Q = [&]( ofstream &ofs )
	{
		//---------------------
		//  < Psi | jY >
		//---------------------
		ofs<< "<OVERLAP_Q>" << std::endl;
		
		for( size_t ib=0; ib!=matrix_V.nr; ++ib )
		{
			for( size_t iat=0; iat!=matrix_Q.size(); ++iat )
			{
				const size_t it = (iat==0) ? TA : TB;
				for( size_t il=0; il!=range_jles[it].size(); ++il )
				{
					for( size_t im=0; im!=range_jles[it][il].M; ++im )
					{
						for( size_t iq=0; iq!=range_jles[it][il].N; ++iq )
						{
							ofs<<matrix_Q[iat]( ib, index_jles[it][il][iq][im] )<<"\t"<<0<<endl;
						}
					}
				}
			}
		}
		
		ofs<< "</OVERLAP_Q>" << std::endl << std::endl;
	};
	
	
	auto print_S = [&]( ofstream &ofs, const double scale=1 )
	{
		//---------------------
		//  < jY | jY >
		//---------------------
		ofs<< "<OVERLAP_Sq>" <<endl;
		
		for( size_t iat1=0; iat1!=matrix_S.size(); ++iat1 )
		{
			const size_t it1 = (iat1==0) ? TA : TB;
			for( size_t il1=0; il1!=range_jles[it1].size(); ++il1 )
			{
				for( size_t im1=0; im1!=range_jles[it1][il1].M; ++im1 )
				{
					for( size_t iat2=0; iat2!=matrix_S[iat1].size(); ++iat2 )
					{
						const size_t it2 = (iat2==0) ? TA : TB;
						for( size_t il2=0; il2!=range_jles[it2].size(); ++il2 )
						{
							for( size_t im2=0; im2!=range_jles[it2][il2].M; ++im2 )
							{
								for( size_t iq1=0; iq1!=range_jles[it1][il1].N; ++iq1 )
								{
									for( size_t iq2=0; iq2!=range_jles[it2][il2].N; ++iq2 )
									{
										ofs<<matrix_S[iat1][iat2]( index_jles[it1][il1][iq1][im1], index_jles[it2][il2][iq2][im2] )*scale<<"\t"<<0<<endl;
									}
								}
							}
						}
					}
				}
			}
		}
		
		ofs<< "</OVERLAP_Sq>" << std::endl << std::endl;
	};
	
	
	auto print_V = [&]( ofstream &ofs, const double scale=1 )
	{
		//---------------------
		//  < Psi | Psi >
		//---------------------	
		ofs << "<OVERLAP_V>" << std::endl;
		
		for( size_t ib1=0; ib1!=matrix_V.nr; ++ib1 )
		{
			for( size_t ib2=0; ib2!=matrix_V.nc; ++ib2 )
			{
				ofs<<matrix_V(ib1,ib2)*scale<<"\t";
			}
			ofs<<endl;
		}
		
		ofs << "</OVERLAP_V>" << std::endl << std::endl;
	};
	
	ofstream ofs(file_name+"_"+TO_STRING(TA)+"_"+TO_STRING(IA)+"_"+TO_STRING(TB)+"_"+TO_STRING(IB));
	print_header(ofs);
	print_Q(ofs);
	print_S(ofs);
	print_V(ofs);
	ofs.close();
}