//==========================================================
// AUTHOR : Wenfei Li
// DATE :   2021-07-30
// UPDATE :
//==========================================================

#ifdef USE_LIBXC

#include "potential_libxc.h"
#include "global.h"
#include "tools.h"
#include "xc_gga_pw.h"
#include "../module_base/global_function.h"
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

// [etxc, vtxc, v, vofk] = Potential_Libxc::v_xc(...)
tuple<double,double,matrix,matrix> Potential_Libxc::v_xc_meta(
	const double * const * const rho_in,
	const double * const rho_core_in,
	const double * const * const kin_r_in)
{
    TITLE("Potential_Libxc","v_xc");
    timer::tick("Potential_Libxc","v_xc");

    double etxc = 0.0;
    double vtxc = 0.0;
	matrix v(GlobalV::NSPIN,GlobalC::pw.nrxx);
	matrix vofk(GlobalV::NSPIN,GlobalC::pw.nrxx);

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
	if(GlobalC::xcf.igcx_now == 13)
	{
		xc_func_init(&x_func, 263 ,xc_polarized);
	}
	else
	{
		throw domain_error("iexch="+TO_STRING(GlobalC::xcf.iexch_now)+", igcx="+TO_STRING(GlobalC::xcf.igcx_now)
			+" unfinished in "+TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));
	}

	//correlation
	if(GlobalC::xcf.igcc_now == 9)
	{
		xc_func_init(&c_func, 267 ,xc_polarized);
	}
    else
	{
		throw domain_error("icorr="+TO_STRING(GlobalC::xcf.icorr_now)+", igcc="+TO_STRING(GlobalC::xcf.igcc_now)
			+" unfinished in "+TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));

	}
	
	//input of xc_mgga_exc_vxc 
	vector<vector<double>> rho;
	vector<vector<Vector3<double>>> grho;
	vector<vector<double>> sigma;
	vector<vector<double>> kin_r;

	//output of xc_mgga_exc_vxc
	vector<vector<double>> vrho;
	vector<vector<Vector3<double>>> h;
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
		GGA_PW::grad_rho(rhog.data(), grho[is].data());

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

	cout << "flag 3" << endl;
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

			lapl = 0.0;//dummy argument, not used

	cout << "flag 4" << endl;
			for( int ir=0; ir!=GlobalC::pw.nrxx; ++ir )
			{
				arho  = abs(rho[is][ir]);
				grho2 = grho[0][ir]*grho[0][ir];
				atau  = kin_r[is][ir] / e2;
			
				if(arho > rho_th && grho2 > grho_th && abs(atau) > tau_th)
				{
					xc_mgga_exc_vxc(&x_func,1,&arho,&grho2,&lapl,&atau,&ex,&v1x,&v2x,&vlapl,&v3x);
					xc_mgga_exc_vxc(&c_func,1,&arho,&grho2,&lapl,&atau,&ec,&v1c,&v2c,&vlapl,&v3c);
					ex = ex * rho[is][ir];
					ec = ec * rho[is][ir];
					v2x = v2x * e2;
					v2c = v2c * e2;

					vrho[is][ir] = (v1x + v1c) * e2;
					h[is][ir] = (v2x + v2c) * e2 * grho[is][ir];
					kedtaur[is][ir] = v3x + v3c;
				
					etxc += (ex + ec) * e2;
					vtxc += (v1x + v1c) * e2 * arho;	
				}
				else
				{
					vrho[is][ir] = 0.0;
					h[is][ir] = 0.0;
					cout << "should be 0" << h[is][ir];
					kedtaur[is][ir] = 0.0;
				}
			}//loop over grid points
	cout << "flag 5" << endl;
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

		double rh, ggrho2, atau, ex, ec;
		vector<double> arho, grho2, tau, lapl;
		vector<double> v1x, v2x, v3x, v1c, v2c, v3c, vlapl;
		Vector3<double> v2cup, v2cdw;
	
		arho.resize(2);grho2.resize(3);tau.resize(2);
		v1x.resize(2);v2x.resize(3);v3x.resize(2);
		v1c.resize(2);v2c.resize(3);v3c.resize(2);
		lapl.resize(6);vlapl.resize(6);
	
		const double rho_th  = 1e-8;
		const double grho_th = 1e-12;
		const double tau_th  = 1e-8;

		for( int ir=0; ir!=GlobalC::pw.nrxx; ++ir )
		{
			arho[0] = rho[0][ir]; arho[1] = rho[1][ir];
			rh = arho[0] + arho[1];

			grho2[0]=grho[0][ir]*grho[0][ir];
			grho2[1]=grho[0][ir]*grho[1][ir];
			grho2[2]=grho[1][ir]*grho[1][ir];
			ggrho2 = (grho2[0] + grho2[2]) * 4.0;

			tau[0] = kin_r[0][ir] / e2;tau[1] = kin_r[1][ir] / e2;
			atau = tau[0] + tau[1];

			if (rh > rho_th && ggrho2 > grho_th && abs(atau) > tau_th)
			{
				xc_mgga_exc_vxc(&x_func,1,arho.data(),grho2.data(),lapl.data(),tau.data(),&ex,v1x.data(),v2x.data(),vlapl.data(),v3x.data());
				xc_mgga_exc_vxc(&c_func,1,arho.data(),grho2.data(),lapl.data(),tau.data(),&ec,v1c.data(),v2c.data(),vlapl.data(),v3c.data());

				ex = ex*rh;
				v2x[0]=v2x[0]*e2;v2x[1]=v2x[1]*e2;
				
				ec = ec*rh;
				v2c[0]=v2c[0]*e2;v2c[1]=v2c[1]*e2;

				v2cup=v2c[0]*grho[0][ir]*e2 + v2c[1]*grho[1][ir];
				v2cdw=v2c[2]*grho[1][ir]*e2 + v2c[1]*grho[0][ir];

				vrho[0][ir] = (v1x[0]+v1c[0]) * e2;
				vrho[1][ir] = (v1x[1]+v1c[1]) * e2;

				h[0][ir] = (v2x[0]*grho[0][ir]+v2cup)*e2;
				h[1][ir] = (v2x[1]*grho[1][ir]+v2cdw)*e2;
				
				kedtaur[0][ir] = v3x[0]+v3c[0];
				kedtaur[1][ir] = v3x[1]+v3c[1];

				etxc += (ex+ec) * e2;
				vtxc += (v1x[0]+v1x[1]+v1c[0]+v1c[1]) * e2 * rh;
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
		WARNING_QUIT("potential_libxc_meta","meta-GGA for nspin=1,2 first");
	}

	vector<double> dh;
	dh.resize(GlobalC::pw.nrxx);
	cout << "flag 6" << endl;

	for(int is=0;is<GlobalV::NSPIN;is++)
	{
		GGA_PW::grad_dot(VECTOR_TO_PTR(h[is]),VECTOR_TO_PTR(dh));	
		for( int ir=0; ir!=GlobalC::pw.nrxx; ++ir )
		{
			vrho[is][ir]-=dh[ir];
			vtxc-=dh[ir]*rho_in[is][ir];
		}
	}
	cout << "flag 7" << endl;

	for( int is=0; is!=nspin0(); ++is )
	{
		for( int ir=0; ir!=GlobalC::pw.nrxx; ++ir )
		{
			v(is,ir) = vrho[is][ir];
			vofk(is,ir) = kedtaur[is][ir];
		}
	}
	cout << "flag 8" << endl;
		
	//-------------------------------------------------
	// for MPI, reduce the exchange-correlation energy
	//-------------------------------------------------
	Parallel_Reduce::reduce_double_pool( etxc );
	Parallel_Reduce::reduce_double_pool( vtxc );

    etxc *= GlobalC::ucell.omega / GlobalC::pw.ncxyz;
    vtxc *= GlobalC::ucell.omega / GlobalC::pw.ncxyz;

    timer::tick("Potential_Libxc","v_xc_meta");
	return std::make_tuple( etxc, vtxc, move(v), move(vofk) );

}

#endif	//ifdef USE_LIBXC

