#ifndef EXX_ABFS_DM_H
#define EXX_ABFS_DM_H

#include "exx_abfs.h"
#include "src_lcao/abfs-vector3_order.h"

#include<vector>
#include<map>
#include<deque>

class Charge_Broyden;

class Exx_Abfs::DM
{
public:
	map<size_t,map<size_t,vector<ComplexMatrix>>> DMk;								// DMk[iat1][iat2][ik](iw1,iw2)
	map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,vector<matrix>>>> DMr;		// DMr[iat1][iat2][box][is](iw1,iw2)
	deque<map<size_t,map<size_t,vector<ComplexMatrix>>>> DMk_pulay_seq;				// DMk_pulay_seq[istep][iat1][iat2][ik](iw1,iw2)
//	double DM_threshold = std::numeric_limits<double>::min();
	bool flag_mix;

public:
	void cal_DM(
		const set<pair<size_t,size_t>> &atom_pairs,
		const vector<Abfs::Vector3_Order<int>> &Born_von_Karman_boxes);
	map<size_t,map<size_t,vector<ComplexMatrix>>> cal_DMk_raw( const set<pair<size_t,size_t>> &atom_pairs ) const;
		
private:
	void cal_DMk_mixing(
		const Charge_Broyden &charge,
		const set<pair<size_t,size_t>> &atom_pairs );
	void plain_mixing(
		const Charge_Broyden &charge,
		const set<pair<size_t,size_t>> &atom_pairs);
	void pulay_mixing(
		const Charge_Broyden &charge,
		const set<pair<size_t,size_t>> &atom_pairs);
		
//	double cal_DM_delta();

//	void   set_DM_threshold(double DM_threshold_in) { DM_threshold = DM_threshold_in; }
//	double get_DM_threshold() const { return DM_threshold; }
};

#endif	// EXX_ABFS_DM_H