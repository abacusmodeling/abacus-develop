//==========================================================
// AUTHOR:	Peize Lin
// DATE:	2021-03-09
//==========================================================
#ifndef VDWD2_PARAMETERS_H
#define VDWD2_PARAMETERS_H

#include <map>
#include <string>
#include "../module_base/vector3.h"
#include "../module_cell/unitcell_pseudo.h"
#include "../input.h"

class Vdwd2_Parameters
{
public:
	bool flag_vdwd2=false;

	Vdwd2_Parameters();

	std::map<std::string,double> C6;
	std::map<std::string,double> R0;

	void C6_input(const std::string &file, const std::string &unit);
	void R0_input(const std::string &file, const std::string &unit);	
		
	double scaling;
	double damping;
	
	std::string model;
	double radius;
	ModuleBase::Vector3<int> period;

	void initset(const UnitCell_pseudo &ucell); //init sets of vdwd2 once this correction is called
	void initial_parameters(const Input &input); //initial parameters of Vdwd2 with INPUT file 

private:

	std::map<std::string,double> C6_default;
	std::map<std::string,double> R0_default;

	void default_C6();
	void default_R0();
};

#endif // define VDWD2_PARAMETERS_H
