#include "exx_abfs-parallel-communicate-hexx.h"
#include "exx_abfs-parallel-communicate-hexx-template.h"
#include "exx_abfs-parallel-communicate-allreduce-template.h"
#include "src_pw/global.h"
#include "src_global/global_function.h"
#include "exx_abfs-io.h"
#include "exx_abfs-io-template.h"

#include "src_external/src_test/src_lcao/exx_lcao-test.h"
//#include <gperftools/profiler.h>

void Exx_Abfs::Parallel::Communicate::Hexx::Rexx_to_Km2D( 
	vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> &HR_exx,
	const pair<bool,bool> &io_HR_a2D )
{
	TITLE("Exx_Abfs::Parallel::Communicate::Hexx::Rexx_to_Km2D");
	
//ofstream ofs_time("time_"+TO_STRING(MY_RANK),ofstream::app);
//timeval t_start;	

//gettimeofday( &t_start, NULL);	
	MPI_Barrier(MPI_COMM_WORLD);		// Peize Lin test
//ofs_time<<"TIME@ MPI_Barrier\t"<<time_during(t_start)<<endl;
//ofs_matrixes( exx_lcao.test_dir+"test-HR_exx_"+TO_STRING(MY_RANK), HR_exx );
	
//gettimeofday( &t_start, NULL);
	vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> HR_a2D;
	if(io_HR_a2D.first)
		HR_a2D = Exx_Abfs::IO::input_binary<vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>>>(
			global_out_dir+"HR_exx_"+TO_STRING(MY_RANK));
	else
	{
		Allreduce allreduce(MPI_COMM_WORLD,HR_exx);
		HR_a2D = allreduce.exx_to_a2D();
	}
	if(io_HR_a2D.second)
		Exx_Abfs::IO::output_binary( HR_a2D, global_out_dir+"HR_exx_"+TO_STRING(MY_RANK) );
//ofs_time<<"TIME@ Exx_Abfs::Parallel::Communicate::Hexx::Allreduce::exx_to_a2D\t"<<time_during(t_start)<<endl;
//ofs_matrixes( exx_lcao.test_dir+"test-HR_a2D_"+TO_STRING(MY_RANK), HR_a2D );
	if(GAMMA_ONLY_LOCAL)
	{
		HK_Gamma_m2D.resize(NSPIN);
		HK_Gamma_m2D_pulay_seq.resize(NSPIN);
		for( size_t is=0; is!=NSPIN; ++is )
		{
//gettimeofday( &t_start, NULL);
			const map<size_t,map<size_t,matrix>> HK_a2D = R_to_K(HR_a2D[is]);
//ofs_time<<"TIME@ Exx_Abfs::Parallel::Communicate::Hexx::R_to_K\t"<<time_during(t_start)<<endl;
//ofs_matrixes( exx_lcao.test_dir+"test-HK_a2D_"+TO_STRING(is)+"_"+TO_STRING(MY_RANK), HK_a2D );
//gettimeofday( &t_start, NULL);
			switch(mixing_mode)
			{
				case Mixing_Mode::No:
					HK_Gamma_m2D[is] = a2D_to_m2D(HK_a2D);	break;
				case Mixing_Mode::Plain:
					if( HK_Gamma_m2D[is].nr && HK_Gamma_m2D[is].nc )
						HK_Gamma_m2D[is] = (1-mixing_beta) * HK_Gamma_m2D[is] + mixing_beta * a2D_to_m2D(HK_a2D);
					else
						HK_Gamma_m2D[is] = a2D_to_m2D(HK_a2D);
					break;
				case Mixing_Mode::Pulay:
					HK_Gamma_m2D[is] = pulay_mixing( HK_Gamma_m2D[is], HK_Gamma_m2D_pulay_seq[is], a2D_to_m2D(HK_a2D) );
					break;
				default:
					throw domain_error(TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));	break;
			}
//ofs_time<<"TIME@ Exx_Abfs::Parallel::Communicate::Hexx::a2D_to_m2D\t"<<time_during(t_start)<<endl;
		}
//ofs_matrixes( exx_lcao.test_dir+"test-HK_Gamma_m2D_"+TO_STRING(MY_RANK), HK_Gamma_m2D );
	}
	else
	{
		HK_K_m2D.resize(kv.nks);
		HK_K_m2D_pulay_seq.resize(kv.nks);
		for( size_t ik=0; ik!=kv.nks; ++ik )
		{
//gettimeofday( &t_start, NULL);	
			const map<size_t,map<size_t,ComplexMatrix>> HK_a2D = R_to_K(HR_a2D[kv.isk[ik]],ik);
//ofs_time<<"TIME@ Exx_Abfs::Parallel::Communicate::Hexx::R_to_K\t"<<time_during(t_start)<<endl;
//ofs_matrixes( exx_lcao.test_dir+"test-HK_a2D_"+TO_STRING(ik)+"_"+TO_STRING(MY_RANK), HK_a2D );
//gettimeofday( &t_start, NULL);	
			switch(mixing_mode)
			{
				case Mixing_Mode::No:
					HK_K_m2D[ik] = a2D_to_m2D(HK_a2D);	break;
				case Mixing_Mode::Plain:
					if( HK_K_m2D[ik].nr && HK_K_m2D[ik].nc )
						HK_K_m2D[ik] = (1-mixing_beta) * HK_K_m2D[ik] + mixing_beta * a2D_to_m2D(HK_a2D);
					else
						HK_K_m2D[ik] = a2D_to_m2D(HK_a2D);
					break;
				case Mixing_Mode::Pulay:
					HK_K_m2D[ik] = pulay_mixing( HK_K_m2D[ik], HK_K_m2D_pulay_seq[ik], a2D_to_m2D(HK_a2D) );
					break;	
				default:
					throw domain_error(TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));	break;
			}
//ofs_time<<"TIME@ Exx_Abfs::Parallel::Communicate::Hexx::a2D_to_m2D\t"<<time_during(t_start)<<endl;
		}
//ofs_matrixes( exx_lcao.test_dir+"test-HK_K_m2D_"+TO_STRING(MY_RANK), HK_K_m2D );
	}
//ofs_time.close();
}


map<size_t,map<size_t,matrix>> Exx_Abfs::Parallel::Communicate::Hexx::R_to_K( 
	map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>> & HR) const
{
	TITLE("Exx_Abfs::Parallel::Communicate::Hexx::R_to_K");
	
	map<size_t,map<size_t,matrix>> HK;
	for( auto & HR1 : HR )
	{
		const size_t iat1 = HR1.first;
		for( auto & HR2 : HR1.second )
		{
			const size_t iat2 = HR2.first;
			assert(HR2.second.size()==1);
			HK[iat1][iat2] = std::move( HR2.second.at(Vector3<int>{0,0,0}) );
		}
	}
	return HK;
}


map<size_t,map<size_t,ComplexMatrix>> Exx_Abfs::Parallel::Communicate::Hexx::R_to_K( 
	const map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>> & HR,
	const size_t ik) const
{
	TITLE("Exx_Abfs::Parallel::Communicate::Hexx::R_to_K");
	
	map<size_t,map<size_t,ComplexMatrix>> HK;
	for( auto & HRA : HR )
	{
		const size_t iat1 = HRA.first;
		for( auto & HRB : HRA.second )
		{
			const size_t iat2 = HRB.first;
			
			ComplexMatrix HK_tmp(
				ucell.atoms[ucell.iat2it[iat1]].nw,
				ucell.atoms[ucell.iat2it[iat2]].nw);
			for( auto & HRC : HRB.second )
			{
				const Abfs::Vector3_Order<int> & box2 = HRC.first;
				HK_tmp += ComplexMatrix(HRC.second) * exp( TWO_PI*IMAG_UNIT * (kv.kvec_c[ik] * (box2*ucell.latvec)) );
			}
			HK[iat1][iat2] = std::move(HK_tmp);
		}
	}
	return HK;
}



