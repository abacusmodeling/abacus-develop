#ifndef EXX_ABFS_PARALLEL_COMMUNICATE_HEXX_TEMPLATE_H
#define EXX_ABFS_PARALLEL_COMMUNICATE_HEXX_TEMPLATE_H

#include "exx_abfs-parallel-communicate-hexx.h"
#include "src_pw/global.h"
#include "src_lcao/global_fp.h"
#include <cassert>

/*
template< typename T >
T Exx_Abfs::Parallel::Communicate::Hexx::a2D_to_m2D( const map<size_t,map<size_t,T>> & H_a2D ) const
{
	TITLE("Exx_Abfs::Parallel::Communicate::Hexx::a2D_to_m2D");
	
	T H_m2D;
	if(KS_SOLVER=="genelpa")
		H_m2D.create( ParaO.ncol, ParaO.nrow );
	else
		H_m2D.create( ParaO.nrow, ParaO.ncol );
	
	for( const auto H_a2DA : H_a2D )
	{
		const size_t iat1 = H_a2DA.first;
		for( const auto H_a2DB : H_a2DA.second )
		{
			const size_t iat2 = H_a2DB.first;
			const T & H = H_a2DB.second;
			
			for( int iw1=0; iw1!=H.nr; ++iw1 )
			{
				const int iwt1 = ucell.itiaiw2iwt( ucell.iat2it[iat1], ucell.iat2ia[iat1], iw1 );
				const int iwt1_m2D = ParaO.trace_loc_row[iwt1];
				if( iwt1_m2D == -1 )	continue;
				
				for( int iw2=0; iw2!=H.nc; ++iw2 )
				{
					const int iwt2 = ucell.itiaiw2iwt( ucell.iat2it[iat2], ucell.iat2ia[iat2], iw2 );
					const int iwt2_m2D = ParaO.trace_loc_col[iwt2];
					if( iwt2_m2D == -1 )	continue;
					
					if(KS_SOLVER=="genelpa")
						H_m2D( iwt2_m2D, iwt1_m2D ) = H(iw1,iw2);
					else
						H_m2D( iwt1_m2D, iwt2_m2D ) = H(iw1,iw2);
				}
			}
		}
	}
	return H_m2D;
}
*/

template<>
matrix Exx_Abfs::Parallel::Communicate::Hexx::H_phase<matrix>(
	matrix &&HR, const int ik, const Abfs::Vector3_Order<int> &box2) const
{
	assert(box2 == Abfs::Vector3_Order<int>(0,0,0));
	return std::forward<matrix>(HR);
}

template<>
ComplexMatrix Exx_Abfs::Parallel::Communicate::Hexx::H_phase<ComplexMatrix>(
	matrix &&HR, const int ik, const Abfs::Vector3_Order<int> &box2) const
{
	return ComplexMatrix(HR) * exp( TWO_PI*IMAG_UNIT * (kv.kvec_c[ik] * (box2*ucell.latvec)) );
}

template<typename Tmatrix>
Tmatrix Exx_Abfs::Parallel::Communicate::Hexx::Ra2D_to_Km2D(
	vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> &HR_a2D, const int ik) const
{
	TITLE("Exx_Abfs::Parallel::Communicate::Hexx::Ra2D_to_Km2D");

	Tmatrix HK_m2D;
	if(KS_SOLVER=="genelpa")
		HK_m2D.create( ParaO.ncol, ParaO.nrow );
	else
		HK_m2D.create( ParaO.nrow, ParaO.ncol );

	const int is_begin = (NSPIN==4) ? 0 : kv.isk[ik];
	const int is_end = (NSPIN==4) ? 4 : kv.isk[ik]+1;
	for(int is=is_begin; is!=is_end; ++is)
	{
		for(auto & HR_a2D_A : HR_a2D[is])
		{
			const size_t iat1 = HR_a2D_A.first;
			const size_t it1 = ucell.iat2it[iat1];
			const size_t ia1 = ucell.iat2ia[iat1];
			for(auto & HR_a2D_B : HR_a2D_A.second)
			{
				const size_t iat2 = HR_a2D_B.first;
				const size_t it2 = ucell.iat2it[iat2];
				const size_t ia2 = ucell.iat2ia[iat2];

				Tmatrix HK_a2D;
				for(auto & HR_a2D_C : HR_a2D_B.second)
				{
					const Abfs::Vector3_Order<int> &box2 = HR_a2D_C.first;
					if(HK_a2D.c)
						HK_a2D += H_phase<Tmatrix>(std::forward<matrix>(HR_a2D_C.second), ik, box2);
					else
						HK_a2D = H_phase<Tmatrix>(std::forward<matrix>(HR_a2D_C.second), ik, box2);
				}

				for(int iw1=0; iw1!=HK_a2D.nr; ++iw1)
				{
					const int iw1_tmp = (NSPIN==4)
						? (iw1*2+is/2)
						: iw1;
					const int iwt1 = ucell.itiaiw2iwt(it1,ia1,iw1_tmp);
					const int iwt1_m2D = ParaO.trace_loc_row[iwt1];
					if(iwt1_m2D<0) continue;
					for(int iw2=0; iw2!=HK_a2D.nc; ++iw2)
					{
						const int iw2_tmp = (NSPIN==4)
							? (iw2*2+is%2)
							: iw2;
						const int iwt2 = ucell.itiaiw2iwt(it2,ia2,iw2_tmp);
						const int iwt2_m2D = ParaO.trace_loc_col[iwt2];
						if(iwt2_m2D<0) continue;

						if(KS_SOLVER=="genelpa" || KS_SOLVER=="scalapack_gvx")
							HK_m2D(iwt2_m2D,iwt1_m2D) = HK_a2D(iw1,iw2);
						else
							HK_m2D(iwt1_m2D,iwt2_m2D) = HK_a2D(iw1,iw2);
					}
				}
			}
		}
	}
	return HK_m2D;
}

template<typename Tmatrix>
void Exx_Abfs::Parallel::Communicate::Hexx::Ra2D_to_Km2D_mixing(
	vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> &HR_a2D,
	vector<Tmatrix> &HK_m2D,
	vector<deque<Tmatrix>> &HK_m2D_pulay_seq) const
{
	TITLE("Exx_Abfs::Parallel::Communicate::Hexx::Ra2D_to_Km2D_mixing");

	HK_m2D.resize(kv.nks);
	HK_m2D_pulay_seq.resize(kv.nks);
	for( int ik=0; ik!=kv.nks; ++ik )
	{
//gettimeofday( &t_start, NULL);
//			const map<size_t,map<size_t,matrix>> HK_a2D = R_to_K(HR_a2D[ik]);
//ofs_time<<"TIME@ Exx_Abfs::Parallel::Communicate::Hexx::R_to_K\t"<<time_during(t_start)<<endl;
//ofs_matrixes( exx_lcao.test_dir+"test-HK_a2D_"+TO_STRING(ik)+"_"+TO_STRING(MY_RANK), HK_a2D );
//gettimeofday( &t_start, NULL);
		switch(mixing_mode)
		{
			case Mixing_Mode::No:
				HK_m2D[ik] = Ra2D_to_Km2D<Tmatrix>(HR_a2D, ik);
				break;
			case Mixing_Mode::Plain:
				if( HK_m2D[ik].nr && HK_m2D[ik].nc )
					HK_m2D[ik] = (1-mixing_beta) * HK_m2D[ik] + mixing_beta * Ra2D_to_Km2D<Tmatrix>(HR_a2D, ik);
				else
					HK_m2D[ik] = Ra2D_to_Km2D<Tmatrix>(HR_a2D, ik);
				break;
			case Mixing_Mode::Pulay:
				HK_m2D[ik] = pulay_mixing( HK_m2D[ik], HK_m2D_pulay_seq[ik], Ra2D_to_Km2D<Tmatrix>(HR_a2D, ik) );
				break;
			default:
				throw domain_error(TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));	break;
		}
//ofs_time<<"TIME@ Exx_Abfs::Parallel::Communicate::Hexx::a2D_to_m2D\t"<<time_during(t_start)<<endl;
	}
//ofs_matrixes( exx_lcao.test_dir+"test-HK_m2D_"+TO_STRING(MY_RANK), HK_m2D );
}

template<typename Tmatrix>
Tmatrix Exx_Abfs::Parallel::Communicate::Hexx::pulay_mixing(
	const Tmatrix &H_pulay_old, deque<Tmatrix> &H_seq, const Tmatrix &H_new ) const
{
	TITLE("Exx_Abfs::Parallel::Communicate::Hexx::pulay_mixing");

	Tmatrix H_pulay;
	if( 0==CHR.totstep )
	{
		H_seq.clear();
		H_pulay = H_new;
	}
	else
	{
		H_seq.push_back( 
			(1-CHR.mixing_beta) * H_pulay_old + CHR.mixing_beta * H_new );
		if( H_seq.size() > CHR.dstep+1 )
			H_seq.pop_front();
		
		if( 1==H_seq.size() )
		{
			H_pulay = H_seq[0];
		}
		else
		{
			auto pos_mod = [](const int i, const int N){ return (i%N+N)%N; };
			auto index = [&](const int i) -> int
			{
				const int alpha_size = H_seq.size()-1;
				const int alpha_begin = CHR.idstep - alpha_size;
				return pos_mod( alpha_begin+i, CHR.dstep );
			};
			H_pulay = (1+CHR.alpha[index(H_seq.size()-2)]) * H_seq[H_seq.size()-1];
			for( int i=1; i<=H_seq.size()-2; ++i )
				H_pulay += ( CHR.alpha[index(i-1)] - CHR.alpha[index(i)] ) * H_seq[i];
			H_pulay -= CHR.alpha[index(0)] * H_seq[0];
		}
	}
	return H_pulay;
}

#endif