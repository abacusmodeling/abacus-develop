// This file contains wrapper for the mGGA functionals
// it includes 2 subroutines:
// 1. tau_xc
// 2. tau_xc_spin

#include "xc_functional.h"

//tau_xc and tau_xc_spin: interface for calling xc_mgga_exc_vxc from LIBXC
//XC_POLARIZED, XC_UNPOLARIZED: internal flags used in LIBXC, denote the polarized(nspin=1) or unpolarized(nspin=2) calculations, definition can be found in xc.h from LIBXC
#ifdef USE_LIBXC
void XC_Functional::tau_xc(const double &rho, const double &grho, const double &atau, double &sx, double &sc,
          double &v1x, double &v2x, double &v3x, double &v1c, double &v2c, double &v3c)
{

	double lapl_rho, vlapl_rho;
	int func_id;
    lapl_rho = grho;
	
//initialize X and C functionals
	xc_func_type x_func;
	xc_func_type c_func;
	const int xc_polarized = XC_UNPOLARIZED;

//exchange
    for(int id : XC_Functional::func_id)
    {
        switch( id )
        {  
            case 263:
                xc_func_init(&x_func, 263 ,xc_polarized);
                xc_mgga_exc_vxc(&x_func,1,&rho,&grho,&lapl_rho,&atau,&sx,&v1x,&v2x,&vlapl_rho,&v3x);
                xc_func_end(&x_func);
                sx = sx * rho;
                v2x = v2x * 2.0;
                break;
            case 267:
                xc_func_init(&c_func, 267 ,xc_polarized);
                xc_mgga_exc_vxc(&c_func,1,&rho,&grho,&lapl_rho,&atau,&sc,&v1c,&v2c,&vlapl_rho,&v3c);
                xc_func_end(&c_func);
                sc = sc * rho;
                v2c = v2c * 2.0;
                break;
            default:
                throw std::domain_error( "functional unfinished in "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
	    }
    }

	return;
}

void XC_Functional::tau_xc_spin(const Mgga_spin_in &mgga_spin_in, Mgga_spin_out &mgga_spin_out)
{
//initialize X and C functionals
	xc_func_type x_func;
	xc_func_type c_func;
	const int xc_polarized = XC_POLARIZED;

	std::vector<double> rho, grho2, tau; //density, gradient and kinetic energy density
	std::vector<double> v1x, v2x, v3x; //exchange potentials
	std::vector<double> v1c, v3c; //correlation potentials; v2c is in output
	std::vector<double> lapl_rho, vlapl_rho; //dummy variables, not used
	double sx, sc;

	rho.resize(2); grho2.resize(3); tau.resize(2);
	v1x.resize(2); v2x.resize(3); v3x.resize(2);
	v1c.resize(3); mgga_spin_out.v2c.resize(3); v3c.resize(2);
	lapl_rho.resize(6); vlapl_rho.resize(6);

	rho[0]=mgga_spin_in.rhoup; rho[1]=mgga_spin_in.rhodw;
	grho2[0]=mgga_spin_in.grhoup*mgga_spin_in.grhoup;
	grho2[1]=mgga_spin_in.grhoup*mgga_spin_in.grhodw;
	grho2[2]=mgga_spin_in.grhodw*mgga_spin_in.grhodw;
	tau[0]=mgga_spin_in.tauup; tau[1]=mgga_spin_in.taudw;

//exchange
    if (func_id[0]==263)
	{
		xc_func_init(&x_func, 263 ,xc_polarized);
		xc_mgga_exc_vxc(&x_func,1,rho.data(),grho2.data(),lapl_rho.data(),tau.data(),&sx,v1x.data(),v2x.data(),vlapl_rho.data(),v3x.data());
		xc_func_end(&x_func);
	}
	else
	{
		ModuleBase::WARNING_QUIT("tau_xc_spin","functional not implemented yet");
	}

//correlation
	if(func_id[1]==267)
	{
		xc_func_init(&c_func, 267 ,xc_polarized);
		xc_mgga_exc_vxc(&c_func,1,rho.data(),grho2.data(),lapl_rho.data(),tau.data(),&sc,v1c.data(),mgga_spin_out.v2c.data(),vlapl_rho.data(),v3c.data());
		xc_func_end(&c_func);
	}
	else
	{
		ModuleBase::WARNING_QUIT("tau_xc_spin","functional not implemented yet");
	}

	mgga_spin_out.ex = sx * (rho[0]+rho[1]);
	mgga_spin_out.v1xup = v1x[0];
	mgga_spin_out.v2xup = v2x[0] * 2.0;
	mgga_spin_out.v3xup = v3x[0];
	mgga_spin_out.v1xdw = v1x[1];
	mgga_spin_out.v2xdw = v2x[2] * 2.0;
	mgga_spin_out.v3xdw = v3x[1];

	mgga_spin_out.ec = sc * (rho[0]+rho[1]);
	mgga_spin_out.v1cup = v1c[0];
	mgga_spin_out.v3cup = v3c[0];
	mgga_spin_out.v1cdw = v1c[1];
	mgga_spin_out.v3cdw = v3c[1];

	mgga_spin_out.v2cup = mgga_spin_out.v2c[0]*mgga_spin_in.grhoup*2.0 
		+ mgga_spin_out.v2c[1]*mgga_spin_in.grhodw;
	mgga_spin_out.v2cdw = mgga_spin_out.v2c[2]*mgga_spin_in.grhodw*2.0 
		+ mgga_spin_out.v2c[1]*mgga_spin_in.grhoup;
	return;
}
#endif