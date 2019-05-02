#ifndef CONV_COULOMB_POT_K_TEMPLATE_H
#define CONV_COULOMB_POT_K_TEMPLATE_H

#include "conv_coulomb_pot_k.h"
#include <vector>
#include <cmath>

template< typename T >
T Conv_Coulomb_Pot_K::cal_orbs_ccp(
	const T & orbs,
	const Ccp_Type &ccp_type,
	const double rmesh_times,
	const double omega,
	const double Rcut)
{
	T orbs_ccp(orbs.size());
	for( size_t i=0; i!=orbs.size(); ++i )
		orbs_ccp[i] = cal_orbs_ccp( orbs[i], ccp_type, rmesh_times, omega, Rcut );
	return orbs_ccp;
}

extern template
Numerical_Orbital_Lm Conv_Coulomb_Pot_K::cal_orbs_ccp<Numerical_Orbital_Lm>(
	const Numerical_Orbital_Lm & orbs,
	const Ccp_Type &ccp_type,
	const double rmesh_times,
	const double omega,
	const double Rcut);
	
#endif