//==========================================================
<<<<<<< HEAD
// AUTHOR : Wenfei Li
// DATE :   2021-07-30
// UPDATE :
=======
// AUTHOR : Peize Lin
// DATE :   2017-09-14
// UPDATE : 2021-02-28
>>>>>>> 82b11b24a6b4ee3d57b66c6058363c91564f079e
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

<<<<<<< HEAD
//the interface to libxc xc_mgga_exc_vxc(xc_func,n,rho,grho,laplrho,tau,e,v1,v2,v3,v4)
//xc_func : LIBXC data type, contains information on xc functional
//n: size of array, nspin*nnr
//rho,grho,laplrho: electron density, its gradient and laplacian
//tau(kin_r): kinetic energy density
//e: energy density
//v1-v4: derivative of energy density w.r.t rho, gradient, laplacian and tau
//v1 and v2 are combined to give v; v4 goes into vofk

// [etxc, vtxc, v, vofk] = Potential_Libxc::v_xc(...)
std::tuple<double,double,matrix,matrix> Potential_Libxc::v_xc_meta(
	const double * const * const rho_in,
	const double * const rho_core_in,
	const double * const * const kin_r_in)
=======
// [etxc, vtxc, v] = Potential_Libxc::v_xc(...)
std::tuple<double,double,matrix,matrix> Potential_Libxc::v_xc_meta(
	const double * const * const rho_in,
	const double * const rho_core_in,
	const double * const * const kin_r)
>>>>>>> 82b11b24a6b4ee3d57b66c6058363c91564f079e
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
	std::vector<xc_func_type> funcs = init_func();

	// the type of rho_sigma_gdr is automatically set to 'tuple'
	// [rho, sigma, gdr] = cal_input( funcs, rho_in );
	const auto rho_sigma_gdr = cal_input( funcs, rho_in, rho_core_in );
	const std::vector<double> &rho = std::get<0>(rho_sigma_gdr);
	const std::vector<double> &sigma = std::get<1>(rho_sigma_gdr);
<<<<<<< HEAD
	std::vector<double> kin_r;
	std::vector<double> lapl_rho; //dummy input, not used in SCAN
	
	//kin_r : from double** to vector<double>
	//lapl_rho : set to 0
	kin_r.resize(GlobalC::pw.nrxx*nspin0());
	lapl_rho.resize(GlobalC::pw.nrxx*nspin0());

	if(nspin0()==1 || GlobalV::NSPIN==2)
	{
		for( size_t is=0; is!=nspin0(); ++is )
			for( size_t ir=0; ir!=GlobalC::pw.nrxx; ++ir )
			{
				kin_r[ir*nspin0()+is] = kin_r_in[is][ir];
				lapl_rho[ir*nspin0()+is] = 0.0;
			}
	}
=======
>>>>>>> 82b11b24a6b4ee3d57b66c6058363c91564f079e

	for( xc_func_type &func : funcs )
	{
		// jiyy add for threshold
		constexpr double rho_threshold = 1E-6;
		constexpr double grho_threshold = 1E-10;
<<<<<<< HEAD
		constexpr double tau_threshold = 1E-6;
=======
>>>>>>> 82b11b24a6b4ee3d57b66c6058363c91564f079e
		xc_func_set_dens_threshold(&func, rho_threshold);
		// sgn for threshold mask
		const std::vector<double> sgn = [&]() -> std::vector<double>
		{
			std::vector<double> sgn( GlobalC::pw.nrxx * nspin0(), 1.0);
			if(nspin0()==2 && func.info->family != XC_FAMILY_LDA && func.info->kind==XC_CORRELATION)
			{
				for( size_t ir=0; ir!=GlobalC::pw.nrxx; ++ir )
				{
<<<<<<< HEAD
					double atau;
					atau = abs(kin_r[ir*2])/2.0;
					if ( rho[ir*2]<rho_threshold || sqrt(abs(sigma[ir*3]))<grho_threshold || atau<tau_threshold) 
						sgn[ir*2] = 0.0;
					atau = abs(kin_r[ir*2+1])/2.0;
					if ( rho[ir*2+1]<rho_threshold || sqrt(abs(sigma[ir*3+2]))<grho_threshold || atau<tau_threshold) 
=======
					if ( rho[ir*2]<rho_threshold || sqrt(abs(sigma[ir*3]))<grho_threshold ) 
						sgn[ir*2] = 0.0;
					if ( rho[ir*2+1]<rho_threshold || sqrt(abs(sigma[ir*3+2]))<grho_threshold ) 
>>>>>>> 82b11b24a6b4ee3d57b66c6058363c91564f079e
						sgn[ir*2+1] = 0.0;
				}
			}
			return sgn;
		}();

<<<<<<< HEAD
		std::vector<double> exc    ( GlobalC::pw.nrxx                       );
		std::vector<double> vrho   ( GlobalC::pw.nrxx * nspin0()            );
		std::vector<double> vsigma ( GlobalC::pw.nrxx * ((1==nspin0())?1:3) );
		std::vector<double> vlapl  ( GlobalC::pw.nrxx * nspin0()            ); //dummy output, not used in SCAN
	    std::vector<double> kedtaur( GlobalC::pw.nrxx * nspin0()			);
	
=======
		std::vector<double> exc   ( GlobalC::pw.nrxx                       );
		std::vector<double> vrho  ( GlobalC::pw.nrxx * nspin0()            );
		std::vector<double> vsigma( GlobalC::pw.nrxx * ((1==nspin0())?1:3) );
		
>>>>>>> 82b11b24a6b4ee3d57b66c6058363c91564f079e
		// cal etxc from rho, exc
		auto process_exc = [&]()
		{
			for( size_t is=0; is!=nspin0(); ++is )
				for( size_t ir=0; ir!=GlobalC::pw.nrxx; ++ir )
					etxc += e2 * exc[ir] * rho[ir*nspin0()+is] * sgn[ir*nspin0()+is];
		};

<<<<<<< HEAD
		auto process_kedtau = [&]()
		{
			if(nspin0()==1 || GlobalV::NSPIN==2)
			{
				for( size_t is=0; is!=nspin0(); ++is )
				{
					for( size_t ir=0; ir!=GlobalC::pw.nrxx; ++ir )
					{
						vofk(is,ir) = kedtaur[ir*nspin0()+is] * sgn[ir*nspin0()+is];
					}
				}
			}
			else // may need updates for SOC
			{
				WARNING_QUIT("Potential_Libxc::v_xc_meta","meta-GGA: implement for nspin=1,2 first");
			}
		};

=======
>>>>>>> 82b11b24a6b4ee3d57b66c6058363c91564f079e
		// cal vtx, v from rho_in, vrho
		auto process_vrho = [&]()
		{
			if(nspin0()==1 || GlobalV::NSPIN==2)
			{
				for( size_t is=0; is!=nspin0(); ++is )
				{
					for( size_t ir=0; ir!=GlobalC::pw.nrxx; ++ir )
					{
						const double v_tmp = e2 * vrho[ir*nspin0()+is] * sgn[ir*nspin0()+is];
						v(is,ir) += v_tmp;
						vtxc += v_tmp * rho_in[is][ir];
					}
				}
			}
			else // may need updates for SOC
			{
				constexpr double vanishing_charge = 1.0e-12;
				for( size_t ir=0; ir!=GlobalC::pw.nrxx; ++ir )
				{
					std::vector<double> v_tmp(4);
					v_tmp[0] = e2 * (0.5 * (vrho[ir*2] + vrho[ir*2+1]));
					const double vs = 0.5 * (vrho[ir*2] - vrho[ir*2+1]);
					const double amag = sqrt( pow(rho_in[1][ir],2) 
						+ pow(rho_in[2][ir],2) 
						+ pow(rho_in[3][ir],2) );

					if(amag>vanishing_charge)
					{
						for(int ipol=1; ipol<4; ++ipol)
						{
							v_tmp[ipol] = e2 * vs * rho_in[ipol][ir] / amag;
						}
					}
					for(int ipol=0; ipol<4; ++ipol)
					{
						v(ipol, ir) += v_tmp[ipol];
						vtxc += v_tmp[ipol] * rho_in[ipol][ir];
					}
				}
			}
			
		};

		// cal vtxc, v from rho_in, rho, gdr, vsigma
		auto process_vsigma = [&]()
		{
			const std::vector<std::vector<Vector3<double>>> &gdr = std::get<2>(rho_sigma_gdr);
			
			std::vector<std::vector<Vector3<double>>> h( nspin0(), std::vector<Vector3<double>>(GlobalC::pw.nrxx) );
			if( 1==nspin0() )
			{
				for( size_t ir=0; ir!=GlobalC::pw.nrxx; ++ir )
				{
					h[0][ir] = e2 * gdr[0][ir] * vsigma[ir] * 2.0 * sgn[ir];
				}
			}
			else
			{
				for( size_t ir=0; ir!=GlobalC::pw.nrxx; ++ir )
				{
					h[0][ir] = e2 * (gdr[0][ir] * vsigma[ir*3  ] * 2.0 * sgn[ir*2  ]
						           + gdr[1][ir] * vsigma[ir*3+1]       * sgn[ir*2]   * sgn[ir*2+1]);
					h[1][ir] = e2 * (gdr[1][ir] * vsigma[ir*3+2] * 2.0 * sgn[ir*2+1]
					               + gdr[0][ir] * vsigma[ir*3+1]       * sgn[ir*2]   * sgn[ir*2+1]);
				}
			}

			// define two dimensional array dh [ nspin, GlobalC::pw.nrxx ]
			std::vector<std::vector<double>> dh(nspin0(), std::vector<double>(GlobalC::pw.nrxx));
			for( size_t is=0; is!=nspin0(); ++is )
				GGA_PW::grad_dot( VECTOR_TO_PTR(h[is]), VECTOR_TO_PTR(dh[is]) );

			for( size_t is=0; is!=nspin0(); ++is )
				for( size_t ir=0; ir!=GlobalC::pw.nrxx; ++ir )
					vtxc -= dh[is][ir] * rho[ir*nspin0()+is];

			if(nspin0()==1 || GlobalV::NSPIN==2)
			{
				for( size_t is=0; is!=nspin0(); ++is )
					for( size_t ir=0; ir!=GlobalC::pw.nrxx; ++ir )
						v(is,ir) -= dh[is][ir];
			}
			else // may need updates for SOC
			{
				constexpr double vanishing_charge = 1.0e-12;
				for( size_t ir=0; ir!=GlobalC::pw.nrxx; ++ir )
				{
					v(0,ir) -= 0.5 * (dh[0][ir] + dh[1][ir]);
					const double amag = sqrt( pow(rho_in[1][ir],2) + pow(rho_in[2][ir],2) + pow(rho_in[3][ir],2) );
					const double neg = (GlobalC::ucell.magnet.lsign_
					                 && rho_in[1][ir]*GlobalC::ucell.magnet.ux_[0]+rho_in[2][ir]*GlobalC::ucell.magnet.ux_[1]+rho_in[3][ir]*GlobalC::ucell.magnet.ux_[2]<=0)
									 ? -1 : 1;
					if(amag > vanishing_charge)
					{
						for(int i=1;i<4;i++)
							v(i,ir) -= neg * 0.5 * (dh[0][ir]-dh[1][ir]) * rho_in[i][ir] / amag;
					}
				}
			}
		};

		//------------------------------------------------------------------
		// "func.info->family" includes LDA, GGA, Hybrid functional (HYB).
		//------------------------------------------------------------------
		switch( func.info->family )
		{
			case XC_FAMILY_LDA:
			case XC_FAMILY_GGA:
			case XC_FAMILY_HYB_GGA:
				WARNING_QUIT("Potential_Libxc::v_xc_meta","func.info->family should be meta-GGA");
<<<<<<< HEAD
				break;
			case XC_FAMILY_MGGA:
				// call Libxc function: xc_mgga_exc_vxc
				xc_mgga_exc_vxc( &func, GlobalC::pw.nrxx, rho.data(), sigma.data(), lapl_rho.data(),kin_r.data(),
					exc.data(), vrho.data(), vsigma.data(), vlapl.data(), kedtaur.data());
				process_exc();
				process_vrho();
				process_vsigma();
				process_kedtau();
				break;
=======
			case XC_FAMILY_MGGA:
				throw domain_error("meta-GGA unfinished in "+TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));
				break;
				// call Libxc function: xc_mgga_exc_vxc
				//xc_mgga_exc_vxc( &func, GlobalC::pw.nrxx, rho.data(), sigma.data(), sigma.data(),
				//	exc.data(), vrho.data(), vsigma.data() );
				//process_exc();
				//process_vrho();
				//process_vsigma();
>>>>>>> 82b11b24a6b4ee3d57b66c6058363c91564f079e
			default:
				throw domain_error("func.info->family ="+TO_STRING(func.info->family)
					+" unfinished in "+TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));
				break;
		}
		// Libxc function: Deallocate memory
		xc_func_end(&func);
	}

	//-------------------------------------------------
	// for MPI, reduce the exchange-correlation energy
	//-------------------------------------------------
	Parallel_Reduce::reduce_double_pool( etxc );
	Parallel_Reduce::reduce_double_pool( vtxc );

    etxc *= GlobalC::ucell.omega / GlobalC::pw.ncxyz;
    vtxc *= GlobalC::ucell.omega / GlobalC::pw.ncxyz;

    timer::tick("Potential_Libxc","v_xc_meta");
	return std::make_tuple( etxc, vtxc, std::move(v), vofk );
}

#endif	//ifdef USE_LIBXC

