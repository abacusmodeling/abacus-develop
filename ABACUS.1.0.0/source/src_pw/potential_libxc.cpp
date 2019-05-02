#include "potential_libxc.h"
#include "src_pw/global.h"
#include "src_lcao/global_fp.h"
#include "src_pw/tools.h"
#include "src_global/global_function.h"
#include "src_pw/gga_pw.h"

void Potential_Libxc::v_xc(
	const double * const * const rho_in,
    double &etxc,
    double &vtxc,
    matrix &v)
{
    TITLE("potential","v_xc");
    timer::tick("Potential_Libxc","v_xc");

    etxc = 0.0;
    vtxc = 0.0;
	v.zero_out();

	vector<XC(func_type)> funcs = init_func();

	const auto input_tmp = cal_input( funcs, rho_in );
	const vector<double> &rho = std::get<0>(input_tmp);
	const vector<double> &sigma = std::get<1>(input_tmp);

	for( XC(func_type) &func : funcs )
	{
		vector<double> exc   ( pw.nrxx                    );
		vector<double> vrho  ( pw.nrxx * NSPIN            );
		vector<double> vsigma( pw.nrxx * ((1==NSPIN)?1:3) );

		auto process_exc = [&]()
		{
			for( size_t is=0; is!=NSPIN; ++is )
				for( size_t ir=0; ir!=pw.nrxx; ++ir )
					etxc += e2 * exc[ir] * rho[ir*NSPIN+is];
		};
		auto process_vrho = [&]()
		{
			for( size_t is=0; is!=NSPIN; ++is )
				for( size_t ir=0; ir!=pw.nrxx; ++ir )
				{
					v(is,ir) += e2 * vrho[ir*NSPIN+is];
					vtxc += e2 * vrho[ir*NSPIN+is] * rho_in[is][ir];
				}
		};
		auto process_vsigma = [&]()
		{
			const std::vector<std::vector<Vector3<double>>> &gdr = std::get<2>(input_tmp);
			
			vector<vector<Vector3<double>>> h( NSPIN, vector<Vector3<double>>(pw.nrxx) );
			if( 1==NSPIN )
				for( size_t ir=0; ir!=pw.nrxx; ++ir )
					h[0][ir] = e2 * gdr[0][ir] * vsigma[ir];
			else
				for( size_t ir=0; ir!=pw.nrxx; ++ir )
				{
					h[0][ir] = e2 * gdr[0][ir] * vsigma[ir*3  ] + gdr[1][ir] * vsigma[ir*3+1];
					h[1][ir] = e2 * gdr[1][ir] * vsigma[ir*3+2] + gdr[0][ir] * vsigma[ir*3+1];
				}
			for( size_t is=0; is!=NSPIN; ++is )
			{
				vector<double> dh(pw.nrxx);
				GGA_PW::grad_dot( VECTOR_TO_PTR(h[is]), VECTOR_TO_PTR(dh) );
				for( size_t ir=0; ir!=pw.nrxx; ++ir )
				{
					v(is,ir) -= e2 * dh[ir];
					vtxc -= e2 * dh[ir] * rho_in[is][ir];
				}
			}
		};

		switch( func.info->family )
		{
			case XC_FAMILY_LDA:
				XC(lda_exc_vxc)( &func, pw.nrxx,VECTOR_TO_PTR(rho), VECTOR_TO_PTR(exc), VECTOR_TO_PTR(vrho) );
				process_exc();
				process_vrho();
				break;
			case XC_FAMILY_GGA:
			case XC_FAMILY_HYB_GGA:
				XC(gga_exc_vxc)( &func, pw.nrxx, VECTOR_TO_PTR(rho), VECTOR_TO_PTR(sigma), VECTOR_TO_PTR(exc), VECTOR_TO_PTR(vrho), VECTOR_TO_PTR(vsigma) );
				process_exc();
				process_vrho();
				process_vsigma();
				break;
			default:
				throw domain_error("func.info->family ="+TO_STRING(func.info->family)+" unfinished "+TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));
				break;
		}

		XC(func_end)(&func);
	}

	Parallel_Reduce::reduce_double_pool( etxc );
	Parallel_Reduce::reduce_double_pool( vtxc );

    etxc *= ucell.omega / pw.ncxyz;
    vtxc *= ucell.omega / pw.ncxyz;

    timer::tick("Potential_Libxc","v_xc");
    return;
}



std::vector<XC(func_type)> Potential_Libxc::init_func()
{
	auto test_parameter_PBE0 = [&]()
	{
		TITLE("test_parameter_PBE0");
		XC(func_type) func;
		XC(func_init)( &func, XC_HYB_GGA_XC_PBEH, XC_POLARIZED );
		
		cout<<XC(hyb_exx_coef)(&func)<<endl;						// 0.25
	};
	auto test_parameter_HSE = [&]()
	{
		TITLE("test_parameter_HSE");
		XC(func_type) func;
		XC(func_init)( &func, XC_HYB_GGA_XC_HSE06, XC_POLARIZED );
		
		double omega=-1,alpha=-1,beta=-1;
		cout<<XC(hyb_exx_coef)(&func)<<endl;						// 0
		XC(hyb_cam_coef)(&func, &omega, &alpha, &beta);
		cout<<omega<<"\t"<<alpha<<"\t"<<beta<<endl;					// 0.11, 0, 0.25
		
		double parameter_hse[3] = {0.123,0.456,0.789};				// beta, omega_HF, omega_PBE 
		xc_func_set_ext_params(&func, parameter_hse);
		
		omega=-2,alpha=-2,beta=-2;
		cout<<XC(hyb_exx_coef)(&func)<<endl;						// 0
		XC(hyb_cam_coef)(&func, &omega, &alpha, &beta);
		cout<<omega<<"\t"<<alpha<<"\t"<<beta<<endl;					// 0.456, 0, 0.123
	};	
	
	
	std::vector<XC(func_type)> funcs;

	const int xc_polarized = (1==NSPIN) ? XC_UNPOLARIZED : XC_POLARIZED;

	auto add_func = [&]( const int function )
	{
		funcs.push_back({});
		XC(func_init)( &funcs.back(), function, xc_polarized );
	};
	
	if( 6==xcf.iexch_now && 8==xcf.igcx_now && 4==xcf.icorr_now && 4==xcf.igcc_now )
	{	
		add_func( XC_HYB_GGA_XC_PBEH );		
		double parameter_hse[3] = { exx_global.info.hybrid_alpha, exx_global.info.hse_omega, exx_global.info.hse_omega };
		xc_func_set_ext_params(&funcs.back(), parameter_hse);
		return funcs;	
	}
	else if( 9==xcf.iexch_now && 12==xcf.igcx_now && 4==xcf.icorr_now && 4==xcf.igcc_now )
	{
		add_func( XC_HYB_GGA_XC_HSE06 );	
		double parameter_hse[3] = { exx_global.info.hybrid_alpha, exx_global.info.hse_omega, exx_global.info.hse_omega };
		xc_func_set_ext_params(&funcs.back(), parameter_hse);
		return funcs;	
	}
	
	if( 1==xcf.iexch_now && 0==xcf.igcx_now )
		add_func( XC_LDA_X );
	else if( 1==xcf.iexch_now && 3==xcf.igcx_now )
		add_func( XC_GGA_X_PBE );
	else
		throw domain_error("iexch="+TO_STRING(xcf.iexch_now)+", igcx="+TO_STRING(xcf.igcx_now)+" unfinished");

	if( 1==xcf.icorr_now && 0==xcf.igcc_now )
		add_func( XC_LDA_C_PZ );
	else if( 4==xcf.icorr_now && 4==xcf.igcc_now )
		add_func( XC_GGA_C_PBE );
	else
		throw domain_error("icorr="+TO_STRING(xcf.icorr_now)+", igcc="+TO_STRING(xcf.igcc_now)+" unfinished");

	return funcs;
}



std::tuple< std::vector<double>, std::vector<double>, std::vector<std::vector<Vector3<double>>> > 
Potential_Libxc::cal_input( 
	const std::vector<XC(func_type)> &funcs, 
	const double * const * const rho_in )
{
	// ..., ↑_{i},↓_{i}, ↑_{i+1},↓_{i+1}, ...
	std::vector<double> rho;
	// ..., ↑↑_{i},↑↓_{i},↓↓_{i}, ↑↑_{i+1},↑↓_{i+1},↓↓_{i+1}, ...
	std::vector<double> sigma;
	// [...,↑_{i},↑_{i+1},...], [...,↓_{i},↓_{i+1},...]
	std::vector<std::vector<Vector3<double>>> gdr;

	bool finished_rho   = false;
	bool finished_gdr   = false;
	bool finished_sigma = false;
	
	auto cal_rho = [&]()
	{
		if(!finished_rho)
		{
			rho.resize( pw.nrxx * NSPIN );
			for( size_t is=0; is!=NSPIN; ++is )
				for( size_t ir=0; ir!=pw.nrxx; ++ir )
				{
					rho[ir*NSPIN+is] = rho_in[is][ir] + 1.0/NSPIN*chr.rho_core[ir];
				}
		}
		finished_rho = true;
	};		
		
	auto cal_gdr = [&]()
	{
		if(!finished_gdr)
		{
			for( size_t is=0; is!=NSPIN; ++is )
				chr.set_rhog(chr.rho[is], chr.rhog[is]);
			chr.set_rhog(chr.rho_core, chr.rhog_core);

			gdr.resize( NSPIN );
			for( size_t is=0; is!=NSPIN; ++is )
			{
				vector<complex<double>>rhog(pw.ngmc);
				for(int ig=0; ig<pw.ngmc; ig++)
					rhog[ig] = chr.rhog[is][ig] + 1.0/NSPIN * chr.rhog_core[ig];

				gdr[is].resize(pw.nrxx);
				GGA_PW::grad_rho( VECTOR_TO_PTR(rhog), VECTOR_TO_PTR(gdr[is]) );
			}
		}
		finished_gdr = true;
	};
	
	auto cal_sigma = [&]()
	{
		if(!finished_sigma)
		{
			sigma.resize( pw.nrxx * ((1==NSPIN)?1:3) );

			if( 1==NSPIN )
				for( size_t ir=0; ir!=pw.nrxx; ++ir )
					sigma[ir] = gdr[0][ir]*gdr[0][ir];
			else
				for( size_t ir=0; ir!=pw.nrxx; ++ir )
				{
					sigma[ir*3]   = gdr[0][ir]*gdr[0][ir];
					sigma[ir*3+1] = gdr[0][ir]*gdr[1][ir];
					sigma[ir*3+2] = gdr[1][ir]*gdr[1][ir];
				}
		}		
		finished_sigma = true;
	};

	
	for( const XC(func_type) &func : funcs )
	{
		switch( func.info->family )
		{
			case XC_FAMILY_LDA:
				cal_rho();
				break;
			case XC_FAMILY_GGA:
			case XC_FAMILY_HYB_GGA:
				cal_rho();
				cal_gdr();
				cal_sigma();
				break;
			default:
				throw domain_error("func.info->family ="+TO_STRING(func.info->family)+" unfinished "+TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));
				break;
		}
	}
	
	return std::make_tuple( rho, sigma, gdr );
}
