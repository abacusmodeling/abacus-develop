//=========================================================
//AUTHOR : Peize Lin 
//DATE : 2016-01-24
//=========================================================

#ifndef CENTER2_ORB_ORB11_H
#define CENTER2_ORB_ORB11_H

#include<vector>
using std::vector;
#include<map>
using std::map;
#include<set>
using std::set;

#include "center2_orb.h"

#include "src_lcao/make_overlap_table.h"
#include "src_lcao/make_gaunt_table.h"
#include "src_lcao/ORB_atomic_lm.h"
#include "src_global/vector3.h"

class Center2_Orb::Orb11
{
public:
	Orb11(
		const Numerical_Orbital_Lm &nA_in,
		const Numerical_Orbital_Lm &nB_in,
		const Make_Overlap_Table &MOT_in,
		const Make_Gaunt_Table &MGT_in	);
	void init_radial_table();
	void init_radial_table( const set<size_t> &radials );				// unit: Bohr/MOT.dr
	double cal_overlap(
		const Vector3<double> &RA, const Vector3<double> &RB,			// unit: Bohr
		const int &mA, const int &mB) const;
		
private:
	const Numerical_Orbital_Lm &nA;
	const Numerical_Orbital_Lm &nB;
	
	const Make_Overlap_Table &MOT;
	const Make_Gaunt_Table &MGT;	

	map<int,vector<double>> Table_r;									// unit: Bohr/MOT.dr
	map<int,vector<double>> Table_dr;
};

// this->Table_r[LAB][ir]

#endif	// CENTER2_ORB_ORB11_H