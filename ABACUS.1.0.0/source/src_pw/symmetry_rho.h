#ifndef SYMMETRY_RHO_H
#define SYMMETRY_RHO_H

#include "symmetry.h"

class Symmetry_rho
{
	public:
	Symmetry_rho();
	~Symmetry_rho();
	
	void begin(const int &spin_now);

	private:

	void psymm(double *rho_part);
	
};

#endif
