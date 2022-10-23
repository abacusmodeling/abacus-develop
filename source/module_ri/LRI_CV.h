//=======================
// AUTHOR : Peize Lin
// DATE :   2022-08-17
//=======================

#ifndef LRI_CV_H
#define LRI_CV_H

#include "Matrix_Orbs11.h"
#include "Matrix_Orbs21.h"
#include "module_orbital/ORB_atomic_lm.h"
#include "src_ri/abfs-vector3_order.h"
#include "module_base/element_basis_index.h"

#include <RI/global/Tensor.h>
#include <RI/global/Global_Func-2.h>

#include <vector>
#include <map>
#include <functional>
#include <pthread.h>

template<typename Tdata>
class LRI_CV
{
private:
	using TA = int;
	using TC = std::array<int,3>;
	using TAC = std::pair<TA,TC>;
	using Tdata_real = Global_Func::To_Real_t<Tdata>;

public:
	LRI_CV();
	~LRI_CV();

	void set_orbitals(
		const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &lcaos_in,
		const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &abfs_in,
		const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &abfs_ccp_in,
		const double &kmesh_times,
		const double &ccp_rmesh_times_in);

	std::map<TA,std::map<TAC,Tensor<Tdata>>>
	cal_Vs(
		const std::vector<TA> &list_A0,
		const std::vector<TAC> &list_A1,
		const double threshold_V,
		const bool flag_writable);
	std::map<TA,std::map<TAC,Tensor<Tdata>>>
	cal_Cs(
		const std::vector<TA> &list_A0,
		const std::vector<TAC> &list_A1,
		const double threshold_C,
		const bool flag_writable);			

private:
	std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> lcaos;
	std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> abfs;
	std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> abfs_ccp;
	ModuleBase::Element_Basis_Index::IndexLNM index_lcaos;
	ModuleBase::Element_Basis_Index::IndexLNM index_abfs;
	double ccp_rmesh_times;

	std::map<int,std::map<int,std::map<Abfs::Vector3_Order<double>,Tensor<Tdata>>>> Vws;
	std::map<int,std::map<int,std::map<Abfs::Vector3_Order<double>,Tensor<Tdata>>>> Cws;
	std::map<int,std::map<int,std::map<Abfs::Vector3_Order<double>,std::array<Tensor<Tdata>,3>>>> dVws;
	pthread_rwlock_t rwlock_Vw;
	pthread_rwlock_t rwlock_Cw;
	pthread_rwlock_t rwlock_dVw;

	Matrix_Orbs11 m_abfs_abfs;
	Matrix_Orbs21 m_abfslcaos_lcaos;

	using T_func_DPcal_data = std::function<Tensor<Tdata>(
		const int it0,
		const int it1,
		const Abfs::Vector3_Order<double> &R,
		const bool flag_writable)>;
	std::map<TA,std::map<TAC,Tensor<Tdata>>>
	cal_datas(
		const std::vector<TA> &list_A0,
		const std::vector<TAC> &list_A1,
		const double threshold,
		const bool flag_writable,
		const T_func_DPcal_data &func_DPcal_data);

	inline Tensor<Tdata>
	DPcal_V(
		const int it0,
		const int it1,
		const Abfs::Vector3_Order<double> &R,
		const bool flag_writable);	
	Tensor<Tdata>
	DPcal_C(
		const int it0,
		const int it1,
		const Abfs::Vector3_Order<double> &R,
		const bool flag_writable);	
	inline std::array<Tensor<Tdata>,3>
	DPcal_dV(
		const int it0,
		const int it1,
		const Abfs::Vector3_Order<double> &R,
		const bool flag_writable);

	template<typename To11, typename Tfunc>
	To11 DPcal_o11(
		const int it0,
		const int it1,
		const Abfs::Vector3_Order<double> &R,
		const bool flag_writable,
		pthread_rwlock_t &rwlock_o11,
		std::map<int,std::map<int,std::map<Abfs::Vector3_Order<double>,To11>>> &o11_ws,
		const Tfunc &func_cal_o11);

	Tensor<Tdata>                           cal_I( const Tensor<Tdata>                           &m  );	
	std::vector<std::vector<Tensor<Tdata>>>	cal_I( const std::vector<std::vector<Tensor<Tdata>>> &ms );	

	inline Tensor<Tdata>               transform(const Tensor<Tdata>               &V ) const;
	inline std::array<Tensor<Tdata>,3> transform(const std::array<Tensor<Tdata>,3> &dV) const;
	
	inline bool exist(const Tensor<Tdata> &V) const;
	inline bool exist(const std::array<Tensor<Tdata>,3> &dV) const;
};

#include "LRI_CV.hpp"

#endif
