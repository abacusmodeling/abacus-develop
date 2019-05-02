#ifndef EXX_ABFS_PARALLEL_COMMUNICATE_HEXX_TEMPLATE_H
#define EXX_ABFS_PARALLEL_COMMUNICATE_HEXX_TEMPLATE_H

#include "exx_abfs-parallel-communicate-hexx.h"

#include "src_pw/global.h"
#include "src_lcao/global_fp.h"

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


template<typename T>
T Exx_Abfs::Parallel::Communicate::Hexx::pulay_mixing( const T &H_pulay_old, deque<T> &H_seq, const T &H_new )
{
	TITLE("Exx_Abfs::Parallel::Communicate::Hexx::pulay_mixing");

	
	T H_pulay;
	if( 0==chr.totstep )
	{
		H_seq.clear();
		H_pulay = H_new;
	}
	else
	{
		H_seq.push_back( 
			(1-chr.mixing_beta) * H_pulay_old + chr.mixing_beta * H_new );
		if( H_seq.size() > chr.dstep+1 )
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
				const int alpha_begin = chr.idstep - alpha_size;
				return pos_mod( alpha_begin+i, chr.dstep );
			};
			H_pulay = (1+chr.alpha[index(H_seq.size()-2)]) * H_seq[H_seq.size()-1];
			for( int i=1; i<=H_seq.size()-2; ++i )
				H_pulay += ( chr.alpha[index(i-1)] - chr.alpha[index(i)] ) * H_seq[i];
			H_pulay -= chr.alpha[index(0)] * H_seq[0];
		}
	}
	return H_pulay;
}

#endif