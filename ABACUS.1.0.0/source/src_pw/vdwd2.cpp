//==========================================================
// AUTHOR : Peize Lin
// DATE : 2014-04-25
// UPDATE : 2019-04-26
//==========================================================

#include"vdwd2.h"
#include"src_global/global_function.h"
#include"src_global/constants.h"
#include<cmath>

Vdwd2::Vdwd2( const UnitCell_pseudo &unitcell):
	ucell(unitcell),
	init_set(false),
	energy_result(0)
{
	init_C6();
	init_R0();
}

void Vdwd2::C6_input(const std::string &file, const std::string &unit)
{
	if( file != "default" )
	{
		ifstream ifs(file);
		if(!ifs)
			WARNING_QUIT("Vdwd2::C6_input", "Can not find the file "+TO_STRING(file));
		std::string element;
		double value;
		while( ifs >> element >> value )
			C6[element]=value;
		ifs.close();
	}
	for(auto &c6 : C6)
	{
		if( unit == "Jnm6/mol")
			c6.second *= 1e6/(ELECTRONVOLT_SI*NA)/pow(BOHR_TO_A,6)/Ry_to_eV;
		else if( unit == "eVA6")
			c6.second /= pow(BOHR_TO_A,6)/Ry_to_eV;
//		else if( unit == "RyBohr6");
		else
			WARNING_QUIT("Input","vdwD2_C6_unit must be Jnm6/mol or eVA6");
	}
}

void Vdwd2::R0_input(const std::string &file, const std::string &unit)
{
	if( file != "default" )
	{
		ifstream ifs(file.c_str());
		if(!ifs)
			WARNING_QUIT("Vdwd2::R0_input", "Can not find the file "+TO_STRING(file));
		std::string element;
		double value;
		while( ifs >> element >> value )
			R0[element]=value;
		ifs.close();
	}
	for(auto &r0 : R0)
	{
		if( unit == "A")
			r0.second/= BOHR_TO_A;
		else if( unit == "Bohr") ;
		else
			WARNING_QUIT("Input","vdwD2_R0_unit must be A or Bohr");			
	}
}

void Vdwd2::initset()
{
	if(!init_set)
	{
		for(auto &c6 : C6)
			c6.second /= pow(ucell.lat0,6);
		for(auto &r0 : R0)
			r0.second /= ucell.lat0;
		if(model=="radius")
		{
			period.x = 2*ceil(radius/ucell.lat0/sqrt(ucell.a1.norm2())) +1;
			period.y = 2*ceil(radius/ucell.lat0/sqrt(ucell.a2.norm2())) +1;
			period.z = 2*ceil(radius/ucell.lat0/sqrt(ucell.a3.norm2())) +1;		
		}
		init_set=true;
	}
}

double Vdwd2::energy()
{
	initset();

	energy_result = 0;
	for( int it1=0; it1!=ucell.ntype; ++it1 )
	{
		for( int it2=0; it2!=ucell.ntype; ++it2 )
		{
			const double C6_product = sqrt( C6.at(ucell.atoms[it1].label) * C6.at(ucell.atoms[it2].label) );
			const double R0_sum = R0.at(ucell.atoms[it1].label) + R0.at(ucell.atoms[it2].label);
			if(!R0_sum)
				WARNING_QUIT("Input", "R0_sum can not be 0");		
			for( int ia1=0; ia1!=ucell.atoms[it1].na; ++ia1 )
			{
				for( int ia2=0; ia2!=ucell.atoms[it2].na; ++ia2 )
				{
					Vector3<int> ilat_loop;
					for( ilat_loop.x = -period.x/2; ilat_loop.x <= (period.x-1)/2; ++ilat_loop.x )
						for( ilat_loop.y = -period.y/2; ilat_loop.y <= (period.y-1)/2; ++ilat_loop.y )
							for( ilat_loop.z = -period.z/2; ilat_loop.z <= (period.z-1)/2; ++ilat_loop.z )
							{
								if( (!( ilat_loop.x || ilat_loop.y || ilat_loop.z )) && (it1==it2) && (ia1==ia2) )
									continue;
								const Vector3<double> tau1 = ucell.atoms[it1].tau[ia1];
								const Vector3<double> tau2 = ucell.atoms[it2].tau[ia2] + ilat_loop * ucell.latvec;
								const double r_sqr = (tau1 - tau2).norm2();
								const double r = sqrt(r_sqr);
								const double tmp_damp_recip = 1+ exp( -damping* (r/R0_sum-1) );
								energy_result -= C6_product/ pow(r_sqr,3)/ tmp_damp_recip/ 2;
							} // end for ilat_loop
				} // end for ia2
			} // end for ia1
		} // end for it2
	} // end for it1
	energy_result *= scaling;
	return energy_result;
}

const std::vector<Vector3<double>> &Vdwd2::force(matrix &stress_result, const bool stress_for_vdw)
{
	initset();

	force_result.clear();
	force_result.resize(ucell.nat);
	if(stress_for_vdw) stress_result.zero_out();
	
	for( int it1=0; it1!=ucell.ntype; ++it1 )
	{
		for( int it2=0; it2!=ucell.ntype; ++it2 )
		{
			const double C6_product = sqrt( C6.at(ucell.atoms[it1].label) * C6.at(ucell.atoms[it2].label) ) ;
			const double R0_sum = R0.at(ucell.atoms[it1].label) + R0.at(ucell.atoms[it2].label) ;
			if(!R0_sum)
				WARNING_QUIT("Input", "R0_sum can not be 0");
			for( int ia1=0; ia1!=ucell.atoms[it1].na; ++ia1 )
			{
				for( int ia2=0; ia2!=ucell.atoms[it2].na; ++ia2 )
				{
					Vector3<int> ilat_loop;
					for( ilat_loop.x = -period.x/2; ilat_loop.x <= (period.x-1)/2; ++ilat_loop.x )
						for( ilat_loop.y = -period.y/2; ilat_loop.y <= (period.y-1)/2; ++ilat_loop.y )
							for( ilat_loop.z = -period.z/2; ilat_loop.z <= (period.z-1)/2; ++ilat_loop.z )
							{
								if( (!( ilat_loop.x || ilat_loop.y || ilat_loop.z )) && (it1==it2) && (ia1==ia2) )
									continue;
								const Vector3<double> tau1 = ucell.atoms[it1].tau[ia1];
								const Vector3<double> tau2 = ucell.atoms[it2].tau[ia2] + ilat_loop * ucell.latvec;
								const double r_sqr = (tau1 - tau2).norm2();
								const double r = sqrt(r_sqr);
								const double tmp_exp = exp( -damping* (r/R0_sum-1) );
								const double tmp_factor = C6_product/ pow(r_sqr,3)/ r/ (1+tmp_exp)* ( -6/r + tmp_exp/(1+tmp_exp)*damping/R0_sum);
								force_result[ucell.itia2iat(it1,ia1)] += tmp_factor*(tau1-tau2);
								if(stress_for_vdw)//added by zhengdy 2018/10/28
								{
									double dr[3]={tau2.x - tau1.x, tau2.y - tau1.y, tau2.z - tau1.z};
									for(int ipol = 0;ipol<3;ipol++)
									{
										for(int jpol = 0;jpol<3;jpol++)
										{
											stress_result(ipol,jpol) += tmp_factor* dr[ipol] *dr[jpol]/2 ;
										}
//										stress_result(ipol,1) += tmp_factor* (tau1.y - tau2.y);
//										stress_result(ipol,2) += tmp_factor* (tau1.z - tau2.z);
									}
								} // end if stress
							} // end for ilat_loop
				} // end for ia2
			} // end for ia1
		} // end for it2
	} // end for it1
	for( int iat=0; iat!=ucell.nat; ++iat )
		force_result[iat] *= scaling/ucell.lat0;
	for(int ipol=0;ipol<3;ipol++)
	{
		for(int jpol=0;jpol<3;jpol++)
		{
			stress_result(ipol,jpol) *= scaling / ucell.omega;
		}
	}

	return force_result;
}
