//=========================================================
//AUTHOR : Peize Lin 
//DATE : 2016-01-24
//=========================================================

#ifndef CENTER2_ORB_ORB22_H
#define CENTER2_ORB_ORB22_H

#include<vector>
using std::vector;
#include<map>
using std::map;
#include<set>
using std::set;

#include "center2_orb.h"
#include "center2_orb-orb21.h"

#include "../module_orbital/ORB_table_phi.h"
#include "../module_orbital/ORB_gaunt_table.h"
#include "../module_orbital/ORB_atomic_lm.h"
#include "../module_base/vector3.h"

class Center2_Orb::Orb22
{
public:
	Orb22(
		const Numerical_Orbital_Lm &nA1_in,
		const Numerical_Orbital_Lm &nA2_in,
		const Numerical_Orbital_Lm &nB1_in,
		const Numerical_Orbital_Lm &nB2_in,		
		const ORB_table_phi &MOT_in,
		const ORB_gaunt_table &MGT_in	);
	void init_radial_table();
	void init_radial_table( const set<size_t> &radials );						// unit: Bohr/MOT.dr
	double cal_overlap(
		const Vector3<double> &RA, const Vector3<double> &RB,					// unit: Bohr
		const int &mA1, const int &mA2, const int &mB1, const int &mB2) const;
		
protected:								// Peize Lin test 2016-10-07
	const Numerical_Orbital_Lm &nA1;
	const Numerical_Orbital_Lm &nA2;
	const Numerical_Orbital_Lm &nB1;
	const Numerical_Orbital_Lm &nB2;
	
	const ORB_table_phi &MOT;
	const ORB_gaunt_table &MGT;	

	map<int,Numerical_Orbital_Lm> nB;
	map<int,Center2_Orb::Orb21> orb21s;
};

// this->orb21s[L34].psi2_center2[L12].Table_r[L1234][ir]

#endif	// CENTER2_ORB_ORB22_H
