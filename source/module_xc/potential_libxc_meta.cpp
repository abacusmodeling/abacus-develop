//==========================================================
// AUTHOR : Wenfei Li
// DATE :   2021-07-30
// UPDATE :
//==========================================================

#ifdef USE_LIBXC

#include "../src_pw/global.h"
#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "module_base/timer.h"
#include "src_parallel/parallel_reduce.h"
#include "./xc_functional.h"
#ifdef __LCAO
#include "../src_lcao/global_fp.h"
#endif

using namespace std;

//the interface to libxc xc_mgga_exc_vxc(xc_func,n,rho,grho,laplrho,tau,e,v1,v2,v3,v4)
//xc_func : LIBXC data type, contains information on xc functional
//n: size of array, nspin*nnr
//rho,grho,laplrho: electron density, its gradient and laplacian
//tau(kin_r): kinetic energy density
//e: energy density
//v1-v4: derivative of energy density w.r.t rho, gradient, laplacian and tau
//v1 and v2 are combined to give v; v4 goes into vofk

//XC_POLARIZED, XC_UNPOLARIZED: internal flags used in LIBXC, denote the polarized(nspin=1) or unpolarized(nspin=2) calculations, definition can be found in xc.h from LIBXC

// [etxc, vtxc, v, vofk] = Potential_Libxc::v_xc(...)
tuple<double,double,ModuleBase::matrix,ModuleBase::matrix> XC_Functional::v_xc_meta(
	const double * const * const rho_in,
	const double * const rho_core_in,
	const double * const * const kin_r_in)
{
    ModuleBase::TITLE("Potential_Libxc","v_xc");
    ModuleBase::timer::tick("Potential_Libxc","v_xc");

	//output of the subroutine
    double etxc = 0.0;
    double vtxc = 0.0;
	ModuleBase::matrix v(GlobalV::NSPIN,GlobalC::pw.nrxx);
	ModuleBase::matrix vofk(GlobalV::NSPIN,GlobalC::pw.nrxx);

	if(GlobalV::VXC_IN_H == 0 )
	{
    	ModuleBase::timer::tick("Potential_Libxc","v_xc_meta");
		return std::make_tuple( etxc, vtxc, move(v), move(vofk) );
	}
	//----------------------------------------------------------
	// xc_func_type is defined in Libxc package
	// to understand the usage of xc_func_type,
	// use can check on website, for example:
	// https://www.tddft.org/programs/libxc/manual/libxc-5.1.x/
	//----------------------------------------------------------

	//initialize X and C functionals
	xc_func_type x_func;
	xc_func_type c_func;
	const int xc_polarized = (GlobalV::NSPIN ? XC_UNPOLARIZED : XC_POLARIZED);

	//exchange
	if( func_id[0] == 263)
	{
		xc_func_init(&x_func, 263 ,xc_polarized);
	}
	else
	{
		throw domain_error("functional unfinished in "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
	}

	//correlation
	if( func_id[1] == 267)
	{
		xc_func_init(&c_func, 267 ,xc_polarized);
	}
    else
	{
		throw domain_error("functional unfinished in "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
	}
	
	//rho,grho,tau
	vector<vector<double>> rho;
	vector<vector<ModuleBase::Vector3<double>>> grho;
	vector<vector<double>> kin_r;

	//dExc/d rho,grho,tau
	vector<vector<double>> vrho;
	vector<vector<ModuleBase::Vector3<double>>> h;
	vector<vector<double>> kedtaur;

	//rho : from double** to vector<double>
	rho.resize(GlobalV::NSPIN);
	for( int is=0; is!=GlobalV::NSPIN; ++is )
	{
		rho[is].resize(GlobalC::pw.nrxx);
		for( int ir=0; ir!=GlobalC::pw.nrxx; ++ir )
		{
			rho[is][ir] = rho_in[is][ir] + (1.0/GlobalV::NSPIN)*rho_core_in[ir]; 
		}
	}

	//grho : calculate gradient	
	grho.resize(GlobalV::NSPIN);
	for( int is=0; is!=GlobalV::NSPIN; ++is )
	{
		grho[is].resize(GlobalC::pw.nrxx);
		
		vector<complex<double>> rhog(GlobalC::pw.ngmc);
		GlobalC::CHR.set_rhog(rho[is].data(), rhog.data());
		XC_Functional::grad_rho(rhog.data(), grho[is].data());
	}
		
	//kin_r : from double** to vector<double>
	kin_r.resize(GlobalV::NSPIN);
	if(GlobalV::NSPIN==1 || GlobalV::NSPIN==2)
	{
		for( int is=0; is!=GlobalV::NSPIN; ++is )
		{
			kin_r[is].resize(GlobalC::pw.nrxx);
			for( int ir=0; ir!=GlobalC::pw.nrxx; ++ir )
			{
				kin_r[is][ir] = kin_r_in[is][ir];
			}
		}	
	}

	if(GlobalV::NSPIN==1)
	{
		vrho.resize(GlobalV::NSPIN);
	    h.resize(GlobalV::NSPIN);
	    kedtaur.resize(GlobalV::NSPIN);

		for( int is=0; is!=GlobalV::NSPIN; ++is )
		{
			vrho[is].resize(GlobalC::pw.nrxx);
			h[is].resize(GlobalC::pw.nrxx);
			kedtaur[is].resize(GlobalC::pw.nrxx);

			double arho, grho2, atau, lapl;
			double ex, ec, v1x, v2x, v3x, v1c, v2c, v3c, vlapl;
			const double rho_th  = 1e-8;
			const double grho_th = 1e-12;
			const double tau_th  = 1e-8;


			for( int ir=0; ir!=GlobalC::pw.nrxx; ++ir )
			{
				arho  = abs(rho_in[is][ir]);
				grho2 = grho[0][ir]*grho[0][ir];
				atau  = kin_r[is][ir] / ModuleBase::e2;
				lapl  = grho2;//dummy argument, not used
			
				if(arho > rho_th && grho2 > grho_th && abs(atau) > tau_th)
				{
					xc_mgga_exc_vxc(&x_func,1,&arho,&grho2,&lapl,&atau,&ex,&v1x,&v2x,&vlapl,&v3x);
					xc_mgga_exc_vxc(&c_func,1,&arho,&grho2,&lapl,&atau,&ec,&v1c,&v2c,&vlapl,&v3c);
					ex = ex * rho[is][ir];
					ec = ec * rho[is][ir];
					v2x = v2x * ModuleBase::e2;
					v2c = v2c * ModuleBase::e2;

					vrho[is][ir] = (v1x + v1c) * ModuleBase::e2;
					h[is][ir] = (v2x + v2c) * ModuleBase::e2 * grho[is][ir];
					kedtaur[is][ir] = v3x + v3c;
				
					etxc += (ex + ec) * ModuleBase::e2;
					vtxc += (v1x + v1c) * ModuleBase::e2 * arho;	
				}
				else
				{
					vrho[is][ir] = 0.0;
					h[is][ir] = 0.0;
					kedtaur[is][ir] = 0.0;
				}
			}//loop over grid points
		}//loop over spin
	}//nspin=1
	else if(GlobalV::NSPIN==2)
	{
		vrho.resize(GlobalV::NSPIN);
	    h.resize(GlobalV::NSPIN);
	    kedtaur.resize(GlobalV::NSPIN);

		for( int is=0; is!=GlobalV::NSPIN; ++is )
		{
			vrho[is].resize(GlobalC::pw.nrxx);
			h[is].resize(GlobalC::pw.nrxx);
			kedtaur[is].resize(GlobalC::pw.nrxx);
		}

		double rh, ggrho2, atau;
		double rhoup, rhodw, tauup, taudw;
		double ex, v1xup, v1xdw, v2xup, v2xdw, v3xup, v3xdw;
		double ec, v1cup, v1cdw, v3cup, v3cdw;
		ModuleBase::Vector3<double> grhoup,grhodw,v2cup,v2cdw;

		XC_Functional::Mgga_spin_in mgga_spin_in;
		XC_Functional::Mgga_spin_out mgga_spin_out;

	
		const double rho_th  = 1e-8;
		const double grho_th = 1e-12;
		const double tau_th  = 1e-8;

		for( int ir=0; ir!=GlobalC::pw.nrxx; ++ir )
		{

			auto set_input_tau_spin = [&]()
			{
				mgga_spin_in.rhoup = rho_in[0][ir];
				mgga_spin_in.rhodw = rho_in[1][ir];
				rh = rho_in[0][ir] + rho_in[1][ir];
	
				mgga_spin_in.grhoup = grho[0][ir];
				mgga_spin_in.grhodw = grho[1][ir];
				ggrho2 = (grho[0][ir]*grho[0][ir] + grho[1][ir]*grho[1][ir]) * 4.0;
	
				mgga_spin_in.tauup = kin_r[0][ir] / ModuleBase::e2;
				mgga_spin_in.taudw = kin_r[1][ir] / ModuleBase::e2;
				atau = (kin_r[0][ir]+kin_r[1][ir])/ ModuleBase::e2;
			};

			auto set_output_tau_spin = [&]()
			{
				vrho[0][ir] = (mgga_spin_out.v1xup+mgga_spin_out.v1cup) * ModuleBase::e2;
				vrho[1][ir] = (mgga_spin_out.v1xdw+mgga_spin_out.v1cdw) * ModuleBase::e2;
	
				h[0][ir] = (mgga_spin_out.v2xup*grhoup+mgga_spin_out.v2cup)*ModuleBase::e2;
				h[1][ir] = (mgga_spin_out.v2xdw*grhodw+mgga_spin_out.v2cdw)*ModuleBase::e2;
				
				kedtaur[0][ir] = mgga_spin_out.v3xup+mgga_spin_out.v3cup;
				kedtaur[1][ir] = mgga_spin_out.v3xdw+mgga_spin_out.v3cdw;
	
				etxc += (mgga_spin_out.ex+mgga_spin_out.ec) * ModuleBase::e2;
				vtxc += (mgga_spin_out.v1xup+mgga_spin_out.v1xdw+mgga_spin_out.v1cup+mgga_spin_out.v1cdw) * ModuleBase::e2 * rh;
			};

			if (rh > rho_th && ggrho2 > grho_th && abs(atau) > tau_th)
			{
				set_input_tau_spin();
				XC_Functional::tau_xc_spin(mgga_spin_in, mgga_spin_out);
				set_output_tau_spin();
			}
			else
			{
				vrho[0][ir] = 0.0;
				vrho[1][ir] = 0.0;

				h[0][ir]=0.0;
				h[1][ir]=0.0;

				kedtaur[0][ir] = 0.0;
				kedtaur[1][ir] = 0.0;
			}
		}//end loop grid points
	}//nspin=2
	else
	{
		ModuleBase::WARNING_QUIT("potential_libxc_meta","meta-GGA for nspin=1,2 first");
	}

	vector<double> dh;
	dh.resize(GlobalC::pw.nrxx);

	for(int is=0;is<GlobalV::NSPIN;is++)
	{
		XC_Functional::grad_dot(ModuleBase::GlobalFunc::VECTOR_TO_PTR(h[is]),ModuleBase::GlobalFunc::VECTOR_TO_PTR(dh));	
		for( int ir=0; ir!=GlobalC::pw.nrxx; ++ir )
		{
			vrho[is][ir]-=dh[ir];
			vtxc-=dh[ir]*rho_in[is][ir];
		}
	}

	for( int is=0; is!=nspin0(); ++is )
	{
		for( int ir=0; ir!=GlobalC::pw.nrxx; ++ir )
		{
			v(is,ir) = vrho[is][ir];
			vofk(is,ir) = kedtaur[is][ir];
		}
	}
		
	//-------------------------------------------------
	// for MPI, reduce the exchange-correlation energy
	//-------------------------------------------------
	Parallel_Reduce::reduce_double_pool( etxc );
	Parallel_Reduce::reduce_double_pool( vtxc );

    etxc *= GlobalC::ucell.omega / GlobalC::pw.ncxyz;
    vtxc *= GlobalC::ucell.omega / GlobalC::pw.ncxyz;

    ModuleBase::timer::tick("Potential_Libxc","v_xc_meta");
	return std::make_tuple( etxc, vtxc, move(v), move(vofk) );

}

#endif	//ifdef USE_LIBXC

