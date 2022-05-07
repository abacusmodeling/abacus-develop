#include "../src_pw/global.h"
#include "../module_base/ylm.h"
#include "global_fp.h" // mohan add 2021-01-30

Grid_Base_Beta::Grid_Base_Beta()
{ 
}

Grid_Base_Beta::~Grid_Base_Beta()
{
}

void Grid_Base_Beta::prepare(
    const ModuleBase::Matrix3 &latvec_in,
    const double& lat0_in)
{
	ModuleBase::TITLE("Grid_Base_Beta","prepare");

	this->lat0 = lat0_in;

	this->latvec0 = latvec_in;
	this->latvec0 *= this->lat0;
	
	return;
}
