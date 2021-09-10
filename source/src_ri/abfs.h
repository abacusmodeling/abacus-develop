#ifndef ABFS_H
#define ABFS_H	

#include "../module_base/vector3.h"
#include "../module_base/matrix.h"
#include "exx_abfs.h"

#include <set>
#include <vector>
#include <map>
#include <memory>
#include <pthread.h>

class Abfs
{
public:
	
	template<typename T> class Vector3_Order;

	// Cs[iat1][iat2][box2]
	// Cws[it1,it1,0][it2,R]
	// Vws[it1,0][it2,R]
	static std::map<size_t,std::map<size_t,std::map<Vector3_Order<int>,std::shared_ptr<ModuleBase::matrix>>>> cal_Cs(
		const std::set<size_t> &atom_centres,
		const Exx_Abfs::Matrix_Orbs11 &m_abfs_abfs,
		const Exx_Abfs::Matrix_Orbs21 &m_abfslcaos_lcaos,
		const ModuleBase::Element_Basis_Index::IndexLNM &index_abfs,
		const ModuleBase::Element_Basis_Index::IndexLNM &index_lcaos,
		const double threshold,
		std::map<size_t,std::map<size_t,std::map<Vector3_Order<double>,std::weak_ptr<ModuleBase::matrix>>>> &Cws,
		std::map<size_t,std::map<size_t,std::map<Vector3_Order<double>,std::weak_ptr<ModuleBase::matrix>>>> &Vws );

	// Vs[iat1][iat2][box2]
	// Vws[it1,0][it2,R]
	static std::map<size_t,std::map<size_t,std::map<Vector3_Order<int>,std::shared_ptr<ModuleBase::matrix>>>> cal_Vs(
		const std::vector<std::pair<size_t,size_t>> &atom_pairs,
		const Exx_Abfs::Matrix_Orbs11 &m_abfs_abfs,
		const ModuleBase::Element_Basis_Index::IndexLNM &index_abfs,
		const double rmesh_times,
		const double threshold,
		std::map<size_t,std::map<size_t,std::map<Vector3_Order<double>,std::weak_ptr<ModuleBase::matrix>>>> &Vws );	
		
	// use matrixes ms to get period mps
	static std::map<Abfs::Vector3_Order<int>,std::shared_ptr<ModuleBase::matrix>> 
		cal_mps(
			const Abfs::Vector3_Order<int> &Born_von_Karman_period,
			const std::map<Vector3_Order<int>,std::shared_ptr<ModuleBase::matrix>> &ms );
	// mps[iat1][iat2][box2]
	static std::map<size_t,std::map<size_t,std::map<Abfs::Vector3_Order<int>,std::shared_ptr<ModuleBase::matrix>>>> 
		cal_mps(
			const Abfs::Vector3_Order<int> &Born_von_Karman_period,
			const std::map<size_t,std::map<size_t,std::map<Vector3_Order<int>,std::shared_ptr<ModuleBase::matrix>>>> &ms );

	// C[it1,it1,0][it2,R]
	static std::shared_ptr<ModuleBase::matrix> DPcal_C( 
		const size_t &it1, 
		const size_t &it2, 
		const Vector3_Order<double> &R, 												// unit: ucell.lat0
		const Exx_Abfs::Matrix_Orbs11 &m_abfs_abfs,
		const Exx_Abfs::Matrix_Orbs21 &m_abfslcaos_lcaos,
		const ModuleBase::Element_Basis_Index::IndexLNM &index_abfs,
		const ModuleBase::Element_Basis_Index::IndexLNM &index_lcaos,
		const double threshold,
		const bool writable,
		pthread_rwlock_t &rwlock_Cw,
		pthread_rwlock_t &rwlock_Vw,
		std::map<size_t,std::map<size_t,std::map<Vector3_Order<double>,std::weak_ptr<ModuleBase::matrix>>>> &Cws,
		std::map<size_t,std::map<size_t,std::map<Vector3_Order<double>,std::weak_ptr<ModuleBase::matrix>>>> &Vws );

	// V[it1,0][it2,R]
	static std::shared_ptr<ModuleBase::matrix> DPcal_V( 
		const size_t &it1, 
		const size_t &it2, 
		const Vector3_Order<double> &R, 												// unit: ucell.lat0
		const Exx_Abfs::Matrix_Orbs11 &m_abfs_abfs,
		const ModuleBase::Element_Basis_Index::IndexLNM &index_abfs,
		const double threshold,
		const bool writable,
		pthread_rwlock_t &rwlock_Vw,
		std::map<size_t,std::map<size_t,std::map<Vector3_Order<double>,std::weak_ptr<ModuleBase::matrix>>>> &Vws );		
		
//	static std::map<size_t,std::vector<Abfs::Vector3_Order<int>>> get_adjs( const size_t &T, const ModuleBase::Vector3<double> &tau_cartesian );
	static std::map<size_t,std::vector<Abfs::Vector3_Order<int>>> get_adjs( const size_t &iat );
	static std::vector<std::map<size_t,std::vector<Abfs::Vector3_Order<int>>>> get_adjs();
	static std::map<std::set<size_t>,std::set<size_t>> get_H_pairs_core_group( const std::vector<std::pair<size_t,size_t>> &atom_pairs );
	static std::set<std::pair<size_t,size_t>> get_H_pairs_core( const std::vector<std::pair<size_t,size_t>> &atom_pairs );
	static std::vector<Vector3_Order<int>> get_Coulomb_potential_boxes( const double rmesh_times );
	static std::vector<Vector3_Order<int>> get_Born_von_Karmen_boxes( const Abfs::Vector3_Order<int> &Born_von_Karman_period );

	static std::shared_ptr<ModuleBase::matrix> cal_I( const std::shared_ptr<ModuleBase::matrix> &m );
	static std::vector<std::vector<std::shared_ptr<ModuleBase::matrix>>> cal_I( const std::vector<std::vector<std::shared_ptr<ModuleBase::matrix>>> &ms );
	
	template<typename T1,typename T2,typename T3,typename T4>
	static void delete_empty_ptrs( std::map<T1,std::map<T2,std::map<T3,std::weak_ptr<T4>>>> &ptrs );
//	template<typename T1,typename T2,typename T3,typename Tmatrix>
//	static void delete_threshold_ptrs( std::map<T1,std::map<T2,std::map<T3,std::shared_ptr<Tmatrix>>>> &ptrs, const double threshold);
//	template<typename T1,typename T2,typename T3,typename Tmatrix>
//	static void delete_threshold_ptrs( std::map<T1,std::map<T2,std::map<T3,Tmatrix>>> &ptrs, const double threshold);
	template<typename Tkey,typename Tmatrix>
	static void delete_threshold_ptrs( std::map<Tkey,Tmatrix> &ptrs, const double threshold );
	template<typename Tkey,typename Tmatrix>
	static void delete_threshold_ptrs( std::map<Tkey,std::shared_ptr<Tmatrix>> &ptrs, const double threshold );
	template<typename Tkey1,typename Tkey2,typename Tvalue>
	static void delete_threshold_ptrs( std::map<Tkey1,std::map<Tkey2,Tvalue>> &ptrs, const double threshold );
//	template<typename Tkey,typename Tvalue>
//	static void delete_threshold_ptrs( std::map<Tkey,Tvalue> &ptrs, const double threshold );

	template<typename T1, typename T2, typename Tother>
	static std::vector<std::pair<T1,T2>> get_atom_pair(const std::map<T1,std::map<T2,Tother>> &m);
};

#endif	// ABFS_H
