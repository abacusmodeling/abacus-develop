#ifndef CONV_COULOMB_POT_K_H
#define CONV_COULOMB_POT_K_H

#include<vector>
#include<limits>

class Conv_Coulomb_Pot_K
{
public:

	enum class Ccp_Type{ Ccp, Hse };
	
	template<typename T> static T cal_orbs_ccp(
		const T & orbs,
		const Ccp_Type &ccp_type,
		const double rmesh_times = 1,
		const double omega = std::numeric_limits<double>::min(),
		const double Rcut = std::numeric_limits<double>::max());
	
private:

	static std::vector<double> cal_psi_ccp( const std::vector<double> & psif );

	static std::vector<double> cal_psi_hse( 
		const std::vector<double> & psif,
		const std::vector<double> & k_radial,
		const double omega);
};

#endif