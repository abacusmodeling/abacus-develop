//==========================================================
// AUTHOR : Peize Lin
// DATE : 2014-04-25
//==========================================================

#include"vdwd2.h"
#include"../src_global/vdwd2_parameters.h"
#include"../src_global/global_function.h"
#include"../src_global/element_name.h"
#include"../src_global/constants.h"


bool VdwD2::vdwD2=false;
double VdwD2::scaling=0.75;
double VdwD2::damping=20;
bool VdwD2::init_set=false;
vector<int> VdwD2::atom_kind;
string VdwD2::model="radius";
double VdwD2::radius=30;
int VdwD2::period[3]={3,3,3};
double VdwD2::energy_result=0;
vector< vector<double> > VdwD2::force_result;

VdwD2::VdwD2( const UnitCell_pseudo &unitcell):
	ucell(unitcell)
{
}

void VdwD2::C6_input(string &file, string &unit)
{
	if( file != "default" )
	{
		ifstream ifs(file.c_str());
		if(!ifs)
		{
			WARNING_QUIT("Input", "Can not find the file containing vdwD2_C6.");		
		}
		C6.clear();
		double value;
		while( ifs >> value )
			C6.push_back(value);
		ifs.close();
	}
	for(int i=0; i!=C6.size(); ++i)
	{
		if( unit == "Jnm6/mol")
			C6[i] *= 1e6/(ELECTRONVOLT_SI*NA)/pow(BOHR_TO_A,6)/Ry_to_eV;
		else if( unit == "eVA6")
			C6[i] /= pow(BOHR_TO_A,6)/Ry_to_eV;
//		else if( unit == "RyBohr6");
		else
			WARNING_QUIT("Input","vdwD2_C6_unit must be Jnm6/mol or eVA6");
	}
	return;
}

void VdwD2::R0_input(string &file, string &unit)
{
	if( file != "default" )
	{
		ifstream ifs(file.c_str());
		if(!ifs)
		{
			WARNING_QUIT("Input", "Can not find the file containing vdwD2_R0.");
		}
		R0.clear();
		double value;
		while( ifs >> value )
			R0.push_back(value);
		ifs.close();
	}
	for(int i=0; i!=R0.size(); ++i)
	{
		if( unit == "A")
			R0[i]/= BOHR_TO_A;
		else if( unit == "Bohr") ;
		else
			WARNING_QUIT("Input","vdwD2_R0_unit must be A or Bohr");			
	}
	return;
}

void VdwD2::atomkind (const UnitCell_pseudo &unitcell)
{
	atom_kind.resize(unitcell.ntype) ;	
	for ( vector<int>::size_type i=0 ; i!=unitcell.ntype ; ++i )
	{
		for (int j=0; j!=element_name.size(); ++j)
		{
			if ( unitcell.atoms[i].label == element_name[j] )
			{
				atom_kind[i] = j ;
				break;
			}
		}
	}
	return ;
}

void VdwD2::initset()
{
	if(!init_set)
	{
		atomkind(ucell);
		for( vector<double>::size_type i=0; i!=C6.size(); ++i)
			C6[i] /= pow(ucell.lat0,6);
		for( vector<double>::size_type i=0; i!=R0.size(); ++i)
			R0[i] /= ucell.lat0;
		if(model=="radius")
		{
			VdwD2::period[0]=2*static_cast<int>(radius/ucell.lat0/sqrt(ucell.a1.norm2())+0.999999999) +1;
			VdwD2::period[1]=2*static_cast<int>(radius/ucell.lat0/sqrt(ucell.a2.norm2())+0.999999999) +1;
			VdwD2::period[2]=2*static_cast<int>(radius/ucell.lat0/sqrt(ucell.a2.norm2())+0.999999999) +1;		
		}
		init_set=true;
	}
	return;
}

double VdwD2::energy()
{
	initset();

	int it1, it2, ia1, ia2, ilat_loop[3];
	double energy(0), C6_product, R0_sum, r_sqr, r, tmp_damp_recip; 
	Vector3<double> tau1, tau2;
	
	for( it1=0; it1!=ucell.ntype; ++it1)
	{
		for( ia1=0; ia1!=ucell.atoms[it1].na; ++ia1)
		{
			for( it2=0; it2!=ucell.ntype; ++it2)
			{
				C6_product = sqrt( C6[ atom_kind[it1] ] * C6[ atom_kind[it2] ] ) ;
				R0_sum = R0[ atom_kind[it1] ] + R0[ atom_kind[it2] ] ;
				if(!R0_sum)
				{
					WARNING_QUIT("Input", "R0_sum can not be 0");
				}
				for( ia2=0; ia2!=ucell.atoms[it2].na; ++ia2)
				{
					for( ilat_loop[0] = -period[0]/2; ilat_loop[0] <= (period[0]-1)/2; ++ilat_loop[0])
						for( ilat_loop[1] = -period[1]/2; ilat_loop[1] <= (period[1]-1)/2; ++ilat_loop[1])
							for( ilat_loop[2] = -period[2]/2; ilat_loop[2] <= (period[2]-1)/2; ++ilat_loop[2])
							{
								if( (!( ilat_loop[0] || ilat_loop[1] || ilat_loop[2] )) && (it1==it2) && (ia1==ia2) )
								{
									continue;
								}
								tau1 = ucell.atoms[it1].tau[ia1];
								tau2 = period_atom_tau(ucell.atoms[it2].tau[ia2], ilat_loop);
								r_sqr = (tau1 - tau2).norm2();
								r = sqrt(r_sqr);
								tmp_damp_recip = 1+ exp( -damping* (r/R0_sum-1) );
								energy -= C6_product/ pow(r_sqr,3)/ tmp_damp_recip/ 2;
							}
				}
			}
		}
	}
	energy*=scaling;
	energy_result=energy;
	return energy;
}

vector< vector<double> > VdwD2::force(matrix &stress_result, const bool stress_for_vdw)
{
	initset();

	force_result.resize(ucell.nat);
	for( int iat=0; iat!=ucell.nat; ++iat )
	{
		force_result[iat].assign(3,0);
	}
	if(stress_for_vdw) stress_result.zero_out();
	
	int it1, it2, ia1, ia2, ilat_loop[3], iat(0);
	double C6_product, R0_sum, r_sqr, r, tmp_exp, tmp_factor;
	Vector3<double> tau1, tau2;	
	for( it1=0; it1!=ucell.ntype; ++it1)
	{
		for( ia1=0; ia1!=ucell.atoms[it1].na; ++ia1)
		{
			for( it2=0; it2!=ucell.ntype; ++it2)
			{
				C6_product = sqrt( C6[ atom_kind[it1] ] * C6[ atom_kind[it2] ] ) ;
				R0_sum = R0[ atom_kind[it1] ] + R0[ atom_kind[it2] ] ;
				if(!R0_sum)
				{
					WARNING_QUIT("Input", "R0_sum can not be 0");
				}				
				for( ia2=0; ia2!=ucell.atoms[it2].na; ++ia2)
				{
					for( ilat_loop[0] = -period[0]/2; ilat_loop[0] <= (period[0]-1)/2; ++ilat_loop[0])
						for( ilat_loop[1] = -period[1]/2; ilat_loop[1] <= (period[1]-1)/2; ++ilat_loop[1])
							for( ilat_loop[2] = -period[2]/2; ilat_loop[2] <= (period[2]-1)/2; ++ilat_loop[2])
							{
								if( (!( ilat_loop[0] || ilat_loop[1] || ilat_loop[2] )) && (it1==it2) && (ia1==ia2) )
								{
									continue;
								}
								tau1 = ucell.atoms[it1].tau[ia1];
								tau2 = period_atom_tau(ucell.atoms[it2].tau[ia2], ilat_loop);
								r_sqr = (tau1 - tau2).norm2();
								r = sqrt(r_sqr);
								tmp_exp = exp( -damping* (r/R0_sum-1) );
								tmp_factor = C6_product/ pow(r_sqr,3)/ r/ (1+tmp_exp)* ( -6/r + tmp_exp/(1+tmp_exp)*damping/R0_sum);
								force_result[iat][0] += tmp_factor* (tau1.x - tau2.x);
								force_result[iat][1] += tmp_factor* (tau1.y - tau2.y);
								force_result[iat][2] += tmp_factor* (tau1.z - tau2.z);
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
								}//end if stress
							}
				}
			}
			++iat;
		}
	}
	for( iat=0; iat!=ucell.nat; ++iat)
	{
		force_result[iat][0] *= scaling/ucell.lat0;
		force_result[iat][1] *= scaling/ucell.lat0;
		force_result[iat][2] *= scaling/ucell.lat0;
	}
	for(int ipol=0;ipol<3;ipol++)
	{
		for(int jpol=0;jpol<3;jpol++)
		{
			stress_result(ipol,jpol) *= scaling / ucell.omega;
		}
	}
	return force_result;
}
