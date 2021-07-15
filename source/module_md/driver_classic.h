#ifndef DRIVER_CLASSIC_H
#define DRIVER_CLASSIC_H

#include "../module_cell/unitcell_pseudo.h"

class Driver_classic
{
	public:
	
	Driver_classic();
	~Driver_classic();

	void init();

	private:

	// reading the parameters
	void reading();

    // convert INPUT parameters for classic MD
    void convert(UnitCell_pseudo &ucell_c);

	// classic MD
	void classic_world();


};

#endif