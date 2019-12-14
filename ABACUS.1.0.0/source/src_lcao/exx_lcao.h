#ifndef EXX_LCAO_H
#define EXX_LCAO_H

#include "src_lcao/abfs.h"
#include "src_lcao/abfs-vector3_order.h"
#include "src_lcao/exx_abfs-matrix_orbs11.h" 
#include "src_lcao/exx_abfs-matrix_orbs21.h" 
#include "src_lcao/exx_abfs-dm.h"
#include "src_lcao/exx_abfs-parallel-communicate-hexx.h"
#include "src_lcao/exx_abfs-parallel-communicate-dm.h"
#include "src_lcao/exx_abfs-parallel-communicate-dm2.h"
#include "src_lcao/exx_abfs-screen-schwarz.h"
#include "src_lcao/exx_abfs-screen-cauchy.h"
#include "src_global/element_basis_index.h"
#include "src_pw/functional.h"
#include "src_pw/exx_global.h"

#include<set>
#include<vector>
#include<map>
#include<deque>
#include<memory>
#include<limits>
#include<atomic>

class Exx_Lcao
{
public:
	struct{ string process; string thread; string matrix; } test_dir;	// Peize Lin test
	Exx_Lcao( const Exx_Global::Exx_Info &info_global );				// Peize Lin test
public:
	void init();
	void cal_exx_ions();
	void cal_exx_elec();
	void cal_exx_elec_nscf();
	void add_Hexx(const size_t ik, const double alpha) const;
private:
	vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> cal_Hexx() const;
	void cal_Hexx_thread(
		atomic<size_t> & i_atom_pair,
		atomic_flag & Hexx_lock,
		vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> & HexxR) const;
	double cal_energy(
		const vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> &HexxR ) const;
	void init_radial_table_ions( const set<size_t> &atom_centres_core, const vector<pair<size_t,size_t>> &atom_pairs_core );
		
public:
	enum class Distribute_Type {Htime,Kmeans2,Kmeans1};
	
public:
	double get_energy() const { return energy; }
	
private:
	
	int kmesh_times = 4;				// Peize Lin test
	
#if EXX_DM==1
	Exx_Abfs::Parallel::Communicate::DM DM_para;
#elif EXX_DM==2
	Exx_Abfs::Parallel::Communicate::DM2 DM_para;
#endif
	Exx_Abfs::Parallel::Communicate::Hexx Hexx_para;
	double energy = 0.0;
	
	vector<vector<vector<Numerical_Orbital_Lm>>> lcaos;
	vector<vector<vector<Numerical_Orbital_Lm>>> abfs;
	vector<vector<vector<Numerical_Orbital_Lm>>> abfs_ccp;
	Element_Basis_Index::IndexLNM index_lcaos;
	Element_Basis_Index::IndexLNM index_abfs;

	Exx_Abfs::Matrix_Orbs11 m_abfs_abfs;
	Exx_Abfs::Matrix_Orbs21 m_abfslcaos_lcaos;	
	
	map<size_t,map<size_t,map<Abfs::Vector3_Order<double>,weak_ptr<matrix>>>> Cws;		// Cws[it1][it2][R]
	map<size_t,map<size_t,map<Abfs::Vector3_Order<double>,weak_ptr<matrix>>>> Vws;		// Vws[it1][it2][R]
	map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,shared_ptr<matrix>>>> Cs;		// Cs[iat1][iat2][box2]
	map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,shared_ptr<matrix>>>> Cps;		// Cs[iat1][iat2][box2]
	map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,shared_ptr<matrix>>>> Vs;		// Vs[iat1][iat2][box2]
	map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,shared_ptr<matrix>>>> Vps;		// Vps[iat1][iat2][box2]

	vector<pair<size_t,size_t>> atom_pairs_core_origin;
	vector<pair<size_t,size_t>> atom_pairs_core;
	set<pair<size_t,size_t>> H_atom_pairs_core;
	Abfs::Vector3_Order<int> Born_von_Karman_period;

	Exx_Abfs::Screen::Schwarz schwarz;
	Exx_Abfs::Screen::Cauchy cauchy;
		
public:
	struct Exx_Info
	{
		const Exx_Global::Hybrid_Type &hybrid_type;
		
		const double &hse_omega;
		
		double pca_threshold = 0;
		vector<string> files_abfs;
		double c_threshold  = 0;
		double v_threshold  = 0;
		double dm_threshold = 0;
		double schwarz_threshold = 0;
		double cauchy_threshold = 0;
		double ccp_threshold = 0;
		double ccp_rmesh_times = 10;
		
		Exx_Lcao::Distribute_Type distribute_type;

		Exx_Info( const Exx_Global::Exx_Info &info_global );
	};
	Exx_Info info;
	
	friend class Local_Orbital_Charge;
};

#endif	// EXX_LCAO_H