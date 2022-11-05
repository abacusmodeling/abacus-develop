//=======================
// AUTHOR : Peize Lin
// DATE :   2022-08-17
//=======================

#ifndef EXX_LRI_HPP
#define EXX_LRI_HPP

#include "Exx_LRI.h"
#include "RI_2D_Comm.h"
#include "RI_Util.h"
#include "src_ri/exx_abfs-construct_orbs.h"
#include "src_ri/exx_abfs-util.h"
#include "src_ri/exx_abfs-io.h"
#include "src_ri/conv_coulomb_pot_k.h"
#include "module_base/tool_title.h"
#include "module_base/timer.h"
#include "src_lcao/serialization_cereal.h"
#include "src_lcao/local_orbital_charge.h"
#include "module_orbital/parallel_orbitals.h"

#include <RI/distribute/Distribute_Equally.h>
#include <RI/global/Map_Operator-3.h>

#include <fstream>
#include <string>

//test
#include "/share/home/linpz/ABACUS/LibRI/unittests/global/Tensor-test.h"
#include "/share/home/linpz/ABACUS/abacus-develop_force/source/src_external/src_test/test_function.h"
#include <iomanip>

template<typename Tdata>
void Exx_LRI<Tdata>::init(const MPI_Comm &mpi_comm_in)
{
	ModuleBase::TITLE("Exx_LRI","init");
	ModuleBase::timer::tick("Exx_LRI", "init");

//	if(GlobalC::exx_info.info_global.separate_loop)
//	{
//		Hexx_para.mixing_mode = Exx_Abfs::Parallel::Communicate::Hexx::Mixing_Mode::No;
//		Hexx_para.mixing_beta = 0;
//	}
//	else
//	{
//		if("plain"==GlobalC::CHR.mixing_mode)
//			Hexx_para.mixing_mode = Exx_Abfs::Parallel::Communicate::Hexx::Mixing_Mode::Plain;
//		else if("pulay"==GlobalC::CHR.mixing_mode)
//			Hexx_para.mixing_mode = Exx_Abfs::Parallel::Communicate::Hexx::Mixing_Mode::Pulay;
//		else
//			throw std::invalid_argument("exx mixing error. exx_separate_loop==false, mixing_mode!=plain or pulay");
//		Hexx_para.mixing_beta = GlobalC::CHR.mixing_beta;
//	}


	this->mpi_comm = mpi_comm_in;

	this->lcaos = Exx_Abfs::Construct_Orbs::change_orbs( GlobalC::ORB, this->info.kmesh_times );

//	#ifdef __MPI
//	Exx_Abfs::Util::bcast( this->info.files_abfs, 0, this->mpi_comm );
//	#endif

	const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>>
		abfs_same_atom = Exx_Abfs::Construct_Orbs::abfs_same_atom( this->lcaos, this->info.kmesh_times, this->info.pca_threshold );		// Peize Lin test
	if(this->info.files_abfs.empty())
		this->abfs = abfs_same_atom;
	else
		this->abfs = Exx_Abfs::IO::construct_abfs( abfs_same_atom, GlobalC::ORB, this->info.files_abfs, this->info.kmesh_times );

	switch(this->info.hybrid_type)
	{
		case Exx_Info::Hybrid_Type::HF:
		case Exx_Info::Hybrid_Type::PBE0:
			this->abfs_ccp = Conv_Coulomb_Pot_K::cal_orbs_ccp( this->abfs, Conv_Coulomb_Pot_K::Ccp_Type::Ccp, {}, this->info.ccp_rmesh_times );		break;
		case Exx_Info::Hybrid_Type::HSE:
			this->abfs_ccp = Conv_Coulomb_Pot_K::cal_orbs_ccp( this->abfs, Conv_Coulomb_Pot_K::Ccp_Type::Hse, {{"hse_omega",this->info.hse_omega}}, this->info.ccp_rmesh_times );	break;
		default:
			throw std::domain_error(ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));	break;
	}

	for( size_t T=0; T!=this->abfs.size(); ++T )
		Exx_Abfs::Lmax = std::max( Exx_Abfs::Lmax, static_cast<int>(this->abfs[T].size())-1 );

	this->cv.set_orbitals(
		this->lcaos, this->abfs, this->abfs_ccp,
		this->info.kmesh_times, this->info.ccp_rmesh_times );

	this->exx_lri.set_csm_threshold(this->info.cauchy_threshold);

//	this->m_abfs_abfs.init( 2, this->kmesh_times, (1+this->info.ccp_rmesh_times)/2.0 );
//	this->m_abfs_abfs.init_radial( abfs_ccp, abfs );
//
//	this->m_abfslcaos_lcaos.init( 1, this->kmesh_times, 1 );
//	this->m_abfslcaos_lcaos.init_radial( abfs_ccp, lcaos, lcaos );
	ModuleBase::timer::tick("Exx_LRI", "init");
}

template<typename Tdata>
void Exx_LRI<Tdata>::cal_exx_ions()
{
	ModuleBase::TITLE("Exx_LRI","cal_exx_ions");
	ModuleBase::timer::tick("Exx_LRI", "cal_exx_ions");

//	init_radial_table_ions( cal_atom_centres_core(atom_pairs_core_origin), atom_pairs_core_origin );

//	this->m_abfsabfs.init_radial_table(Rradial);
//	this->m_abfslcaos_lcaos.init_radial_table(Rradial);

	std::vector<TA> atoms(GlobalC::ucell.nat);
	for(int iat=0; iat<GlobalC::ucell.nat; ++iat)
		atoms[iat] = iat;
	std::map<TA,TatomR> atomsR;
	for(int iat=0; iat<GlobalC::ucell.nat; ++iat)
		atomsR[iat] = RI_Util::Vector3_to_array3( GlobalC::ucell.atoms[ GlobalC::ucell.iat2it[iat] ].tau[ GlobalC::ucell.iat2ia[iat] ] );
	const std::array<TatomR,Ndim> latvec
		= {RI_Util::Vector3_to_array3(GlobalC::ucell.a1),
		   RI_Util::Vector3_to_array3(GlobalC::ucell.a2),
		   RI_Util::Vector3_to_array3(GlobalC::ucell.a3)};
	const std::array<Tcell,Ndim> period = {GlobalC::kv.nmp[0], GlobalC::kv.nmp[1], GlobalC::kv.nmp[2]};

	this->exx_lri.set_parallel(this->mpi_comm, atomsR, latvec, period);

	// std::max(3) for gamma_only, list_A2 should contain cell {-1,0,1}. In the future distribute will be neighbour.
	const std::array<Tcell,Ndim> period_tmp = {std::max(3,GlobalC::kv.nmp[0]), std::max(3,GlobalC::kv.nmp[1]), std::max(3,GlobalC::kv.nmp[2])};
	std::pair<std::vector<TA>, std::vector<std::vector<std::pair<TA,std::array<Tcell,Ndim>>>>>
		list_As = Distribute_Equally::distribute_atoms_periods(this->mpi_comm, atoms, period_tmp, 2);
	const std::vector<TA> list_A1 = std::move(list_As.first);
	const std::vector<TAC> list_A2 = std::move(list_As.second[0]);

	std::map<TA,std::map<TAC,RI::Tensor<Tdata>>> Vs = this->cv.cal_Vs(list_A1, list_A2, {{"writable_Vws",true}});
	this->exx_lri.set_Vs(std::move(Vs), this->info.V_threshold);

	if(this->test_Fexx)
	{
		std::array<std::map<TA,std::map<TAC,RI::Tensor<Tdata>>>,3> dVs = this->cv.cal_dVs(list_A1, list_A2, {{"writable_dVws",true}});
		this->exx_lri.set_dVs(std::move(dVs), this->info.V_threshold);
	}

	std::pair<std::map<TA,std::map<TAC,RI::Tensor<Tdata>>>, std::array<std::map<TA,std::map<TAC,RI::Tensor<Tdata>>>,3>>
		Cs_dCs = this->cv.cal_Cs_dCs(list_A1, list_A2,
			{{"cal_dC",this->test_Fexx}, {"writable_Cws",true}, {"writable_dCws",true}, {"writable_Vws",false}, {"writable_dVws",false}});
	std::map<TA,std::map<TAC,RI::Tensor<Tdata>>> &Cs = std::get<0>(Cs_dCs);
	this->exx_lri.set_Cs(std::move(Cs), this->info.C_threshold);

	if(this->test_Fexx)
	{
		std::array<std::map<TA,std::map<TAC,RI::Tensor<Tdata>>>,3> &dCs = std::get<1>(Cs_dCs);
		this->exx_lri.set_dCs(std::move(dCs), this->info.C_threshold);
	}

	ModuleBase::timer::tick("Exx_LRI", "cal_exx_ions");
}

template<typename Tdata>
void Exx_LRI<Tdata>::cal_exx_elec(const Local_Orbital_Charge &loc, const Parallel_Orbitals &pv)
{
	const bool flag_cal_dCV = true;			// tmp
	ModuleBase::TITLE("Exx_LRI","cal_exx_elec");
	ModuleBase::timer::tick("Exx_LRI", "cal_exx_elec");

	const std::vector<std::tuple<std::set<TA>, std::set<TA>>> judge = RI_2D_Comm::get_2D_judge(pv);

	std::vector<std::map<TA,std::map<TAC,RI::Tensor<Tdata>>>> Ds =
		GlobalV::GAMMA_ONLY_LOCAL
		? RI_2D_Comm::split_m2D_ktoR<Tdata>(loc.dm_gamma, pv)
		: RI_2D_Comm::split_m2D_ktoR<Tdata>(loc.dm_k, pv);

	this->Hexxs.resize(GlobalV::NSPIN);
	this->Eexx = 0;
	for(int is=0; is<GlobalV::NSPIN; ++is)
	{
		if(!this->test_Fexx)
		{
			this->exx_lri.set_Ds(std::move(Ds[is]), this->info.dm_threshold);
			this->exx_lri.cal_Hs();
		}
		else
		{
			this->exx_lri.set_Ds(std::move(Ds[is]), this->info.dm_threshold, std::to_string(is));
			this->exx_lri.cal_Hs({"","",std::to_string(is)});
		}
		this->Hexxs[is] = Communicate_Tensors_Map_Judge::comm_map2_first(this->mpi_comm, std::move(this->exx_lri.Hs), std::get<0>(judge[is]), std::get<1>(judge[is]));
		this->Eexx += this->exx_lri.energy;
		post_process_Hexx(this->Hexxs[is]);
	}
	this->Eexx = post_process_Eexx(this->Eexx);
	ModuleBase::timer::tick("Exx_LRI", "cal_exx_elec");

	for(int is=0; is<GlobalV::NSPIN; ++is)
		if(this->test_Fexx)
		{
std::cout<<__FILE__<<__LINE__<<std::endl;
			this->exx_lri.cal_Fs({"","",std::to_string(is),"",""});
			GlobalV::ofs_running<<"Fexx\t"<<this->exx_lri.force<<std::endl;
		}
}

template<typename Tdata>
void Exx_LRI<Tdata>::post_process_Hexx( std::map<TA, std::map<TAC, RI::Tensor<Tdata>>> &Hexxs_io ) const
{
	ModuleBase::TITLE("Exx_LRI","post_process_Hexx");
	constexpr Tdata frac = -1 * 2;								// why?	Hartree to Ry?
	const std::function<void(RI::Tensor<Tdata>&)>
		multiply_frac = [&frac](RI::Tensor<Tdata> &t)
		{ t = t*frac; };
	Map_Operator::for_each( Hexxs_io, multiply_frac );
}

template<typename Tdata>
Tdata Exx_LRI<Tdata>::post_process_Eexx( const Tdata &Eexx_in ) const
{
	ModuleBase::TITLE("Exx_LRI","post_process_Eexx");
	const Tdata SPIN_multiple = std::map<int,Tdata>{{1,2}, {2,1}, {4,1}}.at(GlobalV::NSPIN);				// why?
	const Tdata frac = - SPIN_multiple;
	return frac * Eexx_in;
}

/*
post_process_old
{
	// D
	const std::map<int,double> SPIN_multiple = {{1,0.5}, {2,1}, {4,1}};							// ???
	DR *= SPIN_multiple.at(NSPIN);

	// H
	HR *= -2;

	// E
	const std::map<int,double> SPIN_multiple = {{1,2}, {2,1}, {4,1}};							// ???
	energy *= SPIN_multiple.at(GlobalV::NSPIN);			// ?
	energy /= 2;					// /2 for Ry
}
*/

template<typename Tdata>
void Exx_LRI<Tdata>::write_Hexxs(const std::string &file_name) const
{
	ModuleBase::TITLE("Exx_LRI","write_Hexxs");
	ModuleBase::timer::tick("Exx_LRI", "write_Hexxs");
	std::ofstream ofs(file_name, std::ofstream::binary);
	cereal::BinaryOutputArchive oar(ofs);
	oar(this->Hexxs);
	ModuleBase::timer::tick("Exx_LRI", "write_Hexxs");
}

template<typename Tdata>
void Exx_LRI<Tdata>::read_Hexxs(const std::string &file_name)
{
	ModuleBase::TITLE("Exx_LRI","read_Hexxs");
	ModuleBase::timer::tick("Exx_LRI", "read_Hexxs");
	std::ifstream ifs(file_name, std::ofstream::binary);
	cereal::BinaryInputArchive iar(ifs);
	iar(this->Hexxs);
	ModuleBase::timer::tick("Exx_LRI", "read_Hexxs");
}



#endif
