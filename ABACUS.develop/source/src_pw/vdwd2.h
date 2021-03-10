//==========================================================
// AUTHOR : Peize Lin
// DATE : 2014-04-25
// UPDATE : 2019-04-26
//==========================================================
#ifndef VDWD2_H
#define VDWD2_H
#include"src_pw/vdwd2_parameters.h"
#include"src_pw/unitcell_pseudo.h"
#include"src_global/vector3.h"
#include<vector>

class Vdwd2
{
public:

	Vdwd2(const UnitCell_pseudo &unit_in, Vdwd2_Parameters &para_in);
	
	double energy_result = 0;
	void cal_energy();
	
	std::vector<Vector3<double>> force_result;
	void cal_force();

	Matrix3 stress_result;
	void cal_stress();

private:

	const UnitCell_pseudo &ucell;
	Vdwd2_Parameters &para;

public:

	bool flag_vdwd2()const{ return para.flag_vdwd2; }
};

#endif



/*
	for(int iat=0;iat<ucell.nat;iat++)
	{
		switch(ipol)
		{
			case 0: force_vdw(iat,ipol) += force_result[iat].x;        break;
			case 1: force_vdw(iat,ipol) += force_result[iat].y;        break;
			case 2: force_vdw(iat,ipol) += force_result[iat].z;        break;
		}
	}
*/