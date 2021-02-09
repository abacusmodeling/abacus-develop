#include "exx_abfs.h"

#include "src_pw/global.h"
#include "src_global/mathzone.h"
#include "src_pw/charge_pulay.h"
#include "src_global/container_operator.h"

#include "src_external/src_test/src_global/matrix-test.h"
#include "src_external/src_test/src_global/complexmatrix-test.h"
#include "src_external/src_test/src_lcao/exx_lcao-test.h"
#include "src_external/src_test/src_lcao/exx_abfs-dm-test.h"

void Exx_Abfs::DM::cal_DM(
	const set<pair<size_t,size_t>> &atom_pairs,
	const vector<Abfs::Vector3_Order<int>> &Born_von_Karman_boxes)
{
	TITLE("Exx_Abfs::DM::cal_DM");
	
	cal_DMk_mixing( CHR, atom_pairs );

	for( const pair<size_t,size_t> & atom_pair : atom_pairs )
	{
		const size_t iat1 = atom_pair.first;
		const size_t iat2 = atom_pair.second;
		const size_t it1 = ucell.iat2it[iat1];
		const size_t it2 = ucell.iat2it[iat2];

		for( const Vector3<int> &box : Born_von_Karman_boxes )
		{
			DMr[iat1][iat2][box] = vector<matrix>( NSPIN, {ucell.atoms[it1].nw,ucell.atoms[it2].nw} );
			for( size_t ik=0; ik!=kv.nks; ++ik )
			{
				DMr[iat1][iat2][box][kv.isk[ik]] += ( DMk[iat1][iat2][ik] * exp( -TWO_PI*IMAG_UNIT* (kv.kvec_c[ik]* (box*ucell.latvec)) ) ).real();
			}
		}
	}
}



void Exx_Abfs::DM::cal_DMk_mixing(
	const Charge_Broyden &charge,
	const set<pair<size_t,size_t>> &atom_pairs )
{
	TITLE("Exx_Abfs::DM::cal_DMk_mixing");

	if(flag_mix)
	{
		if ( charge.mixing_mode == "plain" )
		{
			plain_mixing( charge, atom_pairs );
		}
		else if ( charge.mixing_mode == "kerker" )
		{
			throw invalid_argument("mixing density matrix can't be kerker. In "+TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));
		}
		else if ( charge.mixing_mode == "pulay" )
		{
			pulay_mixing( charge, atom_pairs );
		}
		else if ( charge.mixing_mode == "pulay-kerker" )
		{
			throw invalid_argument("mixing density matrix can't be pulay-kerker. In "+TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));
		}
		else
		{
			throw invalid_argument("mixing density matrix error. In "+TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));
		}
	}
	else
	{
		DMk = cal_DMk_raw(atom_pairs);
//		DMk = Exx_Abfs_DM_Test::cal_DMk_raw_readfile(atom_pairs);			// Peize Lin test 2018-03-27

		#if TEST_EXX_LCAO==1
			static int istep=0;
			ofs_matrixes("DMk_"+TO_STRING(istep++)+".dat(@Exx_Abfs::DM::cal_DMk_mixing)",DMk);
		#elif TEST_EXX_LCAO==-1
			#error
		#endif
	}
}



map<size_t,map<size_t,vector<ComplexMatrix>>> Exx_Abfs::DM::cal_DMk_raw( const set<pair<size_t,size_t>> &atom_pairs ) const
{
	TITLE("Exx_Abfs::DM::cal_DMk_raw");

	const double SPIN_multiple = 0.5*NSPIN;
	
	map<size_t,map<size_t,vector<ComplexMatrix>>> DMk_raw;
	for( const pair<size_t,size_t> & atom_pair : atom_pairs )
	{
		const size_t iat1 = atom_pair.first;
		const size_t iat2 = atom_pair.second;
		const size_t it1 = ucell.iat2it[iat1];
		const size_t it2 = ucell.iat2it[iat2];
		const size_t ia1 = ucell.iat2ia[iat1];
		const size_t ia2 = ucell.iat2ia[iat2];

		DMk_raw[iat1][iat2] = vector<ComplexMatrix>( kv.nks, {ucell.atoms[it1].nw,ucell.atoms[it2].nw} );
		for( size_t ik=0; ik!=kv.nks; ++ik )
		{
			for( size_t iw1=0; iw1!=ucell.atoms[it1].nw; ++iw1 )
			{
				for( size_t iw2=0; iw2!=ucell.atoms[it2].nw; ++iw2 )
				{
					for( size_t ib=0; ib!=NBANDS; ++ib )
					{
						if( GAMMA_ONLY_LOCAL )
						{
							//---------------------------------------------------------
							// LOWF.WFC_GAMMA has been replaced by wfc_dm_2d.cpp 
							// we need to fix this function in near future.
							// -- mohan add 2021-02-09
							//---------------------------------------------------------
							WARNING_QUIT("Exx_Abfs::DM::cal_DMk_raw","need to update LOWF.WFC_GAMMA");
//							DMk_raw[iat1][iat2][ik](iw1,iw2) += wf.wg(ik,ib) 
//								* LOWF.WFC_GAMMA[ik][ib][ucell.itiaiw2iwt(it1,ia1,iw1)] 
//								* LOWF.WFC_GAMMA[ik][ib][ucell.itiaiw2iwt(it2,ia2,iw2)];
						}
						else
						{
							DMk_raw[iat1][iat2][ik](iw1,iw2) += wf.wg(ik,ib) 
								* LOWF.WFC_K[ik][ib][ucell.itiaiw2iwt(it1,ia1,iw1)] 
								* conj(LOWF.WFC_K[ik][ib][ucell.itiaiw2iwt(it2,ia2,iw2)]);
						}
					}
				}
			}
			DMk_raw[iat1][iat2][ik] *= SPIN_multiple;
		}
	}
	
	#if TEST_EXX_LCAO==1
		static int istep=0;
		ofs_matrixes("DMk_raw_"+TO_STRING(istep++)+".dat(@Exx_Abfs::DM::cal_DMk_raw)",DMk_raw);
	#elif TEST_EXX_LCAO==-1
		#error
	#endif
	return DMk_raw;
}



void Exx_Abfs::DM::plain_mixing(
	const Charge_Broyden &charge,
	const set<pair<size_t,size_t>> &atom_pairs)
{
	TITLE("Exx_Abfs::DM::plain_mixing");

	if(DMk.empty())
		DMk = cal_DMk_raw(atom_pairs);
	else
		DMk = charge.mixing_beta * cal_DMk_raw(atom_pairs) + (1-charge.mixing_beta) * DMk;
}



void Exx_Abfs::DM::pulay_mixing(
	const Charge_Broyden &charge,
	const set<pair<size_t,size_t>> &atom_pairs)
{
	if( 1==charge.totstep )
	{
		DMk_pulay_seq.clear();
	}
	
	DMk_pulay_seq.push_back( charge.mixing_beta * cal_DMk_raw(atom_pairs) + (1-charge.mixing_beta) * DMk );
	if( charge.totstep > charge.rstep )
		DMk_pulay_seq.pop_front();
	
	if( 1==charge.totstep )
	{
		DMk = DMk_pulay_seq.front();
	}
	else
	{
		const int alpha_size = DMk_pulay_seq.size()-1;
		const int alpha_end = charge.idstep;
		const int alpha_begin = alpha_end - alpha_size;
		auto alpha_index = [&](const int i){ return ((alpha_begin+i)%charge.dstep+charge.dstep)%charge.dstep; };
		
		DMk = (1+charge.alpha[alpha_index(alpha_size-1)]) * DMk_pulay_seq.back();
		for( size_t i=1; i<DMk_pulay_seq.size()-1; ++i )
			DMk = DMk + ( charge.alpha[alpha_index(i-1)] - charge.alpha[alpha_index(i)] ) * DMk_pulay_seq[i];
		DMk = DMk - charge.alpha[alpha_index(0)] * DMk_pulay_seq.front();
				
		#if TEST_EXX_LCAO==1
		{
			cout<<"charge.alpha"<<endl;
			cout<<alpha_begin<<"\t"<<alpha_end<<"\t"<<alpha_size<<endl;
			cout<<charge.alpha[alpha_index(alpha_size-1)]<<endl;
			for( size_t i=1; i<DMk_pulay_seq.size()-1; ++i )
				cout<<charge.alpha[alpha_index(i-1)]<<"\t"<<charge.alpha[alpha_index(i)]<<endl;
			cout<<charge.alpha[alpha_index(0)]<<endl;
		}
		#elif TEST_EXX_LCAO==-1
			#error
		#endif
	}
	
	#if TEST_EXX_LCAO==1
	{
		cout<<"charge.alpha_all"<<endl;
		for( size_t i=0; i!=charge.dstep; ++i )
			cout<<charge.alpha[i]<<"\t";
		cout<<endl;
	}
	#elif TEST_EXX_LCAO==-1
		#error
	#endif	
	
	#if TEST_EXX_LCAO==1
//	{
//		static int istep=0;
//		for( size_t i=0; i!=DMk_pulay_seq.size(); ++i )
//			ofs_matrixes( "DMk_pulay_seq_"+TO_STRING(istep)+"_"+TO_STRING(i), DMk_pulay_seq[i] );
//		++istep;
//	}
	#elif TEST_EXX_LCAO==-1
		#error
	#endif		
}
