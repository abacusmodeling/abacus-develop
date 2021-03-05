//==========================================================
// AUTHOR : Peize Lin
// DATE : 2014-04-25
// UPDATE : 2019-04-26
//==========================================================
#ifndef VDWD2_H
#define VDWD2_H
#include"input_conv.h"
#include"src_pw/unitcell_pseudo.h"
#include"src_global/vector3.h"
#include<string>
#include<vector>

class Vdwd2
{
public:

	Vdwd2( const UnitCell_pseudo &unitcell);

	bool vdwD2;	
	
	double energy_result;
	double energy();
	
	std::vector<Vector3<double>> force_result;
	const std::vector<Vector3<double>> &force(
		const bool force_for_vdw, 
		const bool stress_for_vdw, 
		matrix &force_vdw, 
		matrix &stress_result );
	
private:
	
	double scaling;
	double damping;
	
	std::string model;
	double radius;
	Vector3<int> period;
	
	map<std::string,double> C6;
	map<std::string,double> R0;
	void init_C6();
	void init_R0();
	void C6_input(const std::string &file, const std::string &unit);
	void R0_input(const std::string &file, const std::string &unit);	

	const UnitCell_pseudo &ucell;

	bool init_set;	
	void initset();
	
	friend void Input_Conv::Convert();
};

#endif
