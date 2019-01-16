//==========================================================
// AUTHOR : Peize Lin
// DATE : 2014-04-25
//==========================================================
#ifndef VDWD2_H
#define VDWD2_H
#include"../input_conv.h"
#include"unitcell_pseudo.h"
#include"../src_global/vector3.h"
#include<cmath>
#include<string>
#include<vector>

class VdwD2
{
public:

	VdwD2( const UnitCell_pseudo &unitcell);

	static bool vdwD2;	
	
	static double energy_result;
	double energy();
	
	static vector< vector<double> > force_result;
	vector< vector<double> > force(matrix &stress_result, const bool stress_for_vdw);
	
private:
	
	static double scaling;
	static double damping;
	
	static string model;
	static double radius;
	static int period[3];
	
	static const double C6_tmp[];
	static vector<double> C6;
	static void C6_input(string &file, string &unit);
	static const double R0_tmp[];
	static vector<double> R0;
	static void R0_input(string &file, string &unit);	

	const UnitCell_pseudo &ucell;	
	static vector<int> atom_kind;	
	static void atomkind (const UnitCell_pseudo &unitcell);

	static bool init_set;	
	void initset();
	
	inline Vector3<double> period_atom_tau( Vector3<double> &tau1, int *lat_loop);
	
	friend void Input_Conv::Convert();
};

inline Vector3<double> VdwD2::period_atom_tau ( Vector3<double> &tau, int *lat_loop )
{
	Vector3<double> true_tau;
	true_tau.x = tau.x + lat_loop[0] * ucell.latvec.e11 + lat_loop[1] * ucell.latvec.e21 + lat_loop[2] * ucell.latvec.e31 ;
	true_tau.y = tau.y + lat_loop[0] * ucell.latvec.e12 + lat_loop[1] * ucell.latvec.e22 + lat_loop[2] * ucell.latvec.e32 ;
	true_tau.z = tau.z + lat_loop[0] * ucell.latvec.e13 + lat_loop[1] * ucell.latvec.e23 + lat_loop[2] * ucell.latvec.e33 ;
	return true_tau;
}

#endif
