//=======================
// AUTHOR : Peize Lin
// DATE :   2022-08-17
//=======================

#ifndef EXX_LRI_H
#define EXX_LRI_H

#include "module_xc/exx_info.h"
#include "module_orbital/ORB_atomic_lm.h"
#include "LRI_CV.h"
#include <RI/physics/Exx.h>

#include <vector>
#include <array>
#include <map>

template<typename Tdata>
class Exx_LRI
{
private:
	using TA = int;
	using Tcell = int;
	static const size_t Ndim = 3;
	using TC = std::array<Tcell,Ndim>;
	using TAC = std::pair<TA,TC>;
	using TatomR = std::array<double,Ndim>;		// tmp
	
public:
	Exx_LRI( const Exx_Info::Exx_Info_RI &info_in ) :info(info_in){}

	void init(const MPI_Comm &mpi_comm_in);
	void cal_exx_ions();
	void cal_exx_elec(const Local_Orbital_Charge &loc, const Parallel_Orbitals &pv);

	std::vector< std::map<TA, std::map<TAC, Tensor<Tdata>>>> Hexxs;
	Tdata get_energy() const;

private:
	const Exx_Info::Exx_Info_RI &info;
	MPI_Comm mpi_comm;

	std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> lcaos;
	std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> abfs;
	std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> abfs_ccp;

	LRI_CV<Tdata> cv;
	Exx<TA,Tcell,Ndim,Tdata> exx_lri;

	void post_process_Hexx(std::map<TA, std::map<TAC, Tensor<Tdata>>> &Hexxs) const;	
};

#include "Exx_LRI.hpp"

#endif