#include "xc_functional.h"
#include "../src_pw/global.h"
#include "../module_base/global_function.h"
#include <stdexcept>
#include <xc.h>

XC_Functional::XC_Functional(){}

XC_Functional::~XC_Functional(){}

std::vector<int> XC_Functional::func_id(1);
int XC_Functional::func_type = 0;
bool XC_Functional::use_libxc = true;

int XC_Functional::get_func_type()
{
    return func_type;
}

// The setting values of functional id according to the index in LIBXC
// for detail, refer to https://www.tddft.org/programs/libxc/functionals/
void XC_Functional::set_xc_type(const std::string xc_func_in)
{
    func_id.clear();
    std::string xc_func = xc_func_in;
    std::transform(xc_func.begin(), xc_func.end(), xc_func.begin(), (::toupper));
	if( xc_func == "LDA" || xc_func == "PZ" || xc_func == "SLAPZNOGXNOGC") //SLA+PZ
	{
        // I should use XC_LDA_X and XC_LDA_C_PZ here,
        // but since we are not compiling with libxc as default,
        // I chose to set it manually instead.
        // Same for the rest.
        func_id.push_back(XC_LDA_X);
        func_id.push_back(XC_LDA_C_PZ);
        func_type = 1;
        use_libxc = false;
	}
    else if (xc_func == "PWLDA")
    {
        func_id.push_back(XC_LDA_X);
        func_id.push_back(XC_LDA_C_PW);
        func_type = 1;
        use_libxc = false;
    }
	else if ( xc_func == "PBE" || xc_func == "SLAPWPBXPBC") //PBX+PBC
	{
        func_id.push_back(XC_GGA_X_PBE);
        func_id.push_back(XC_GGA_C_PBE);
        func_type = 2;
        use_libxc = false;
	} 
	else if( xc_func == "revPBE" ) //rPBX+PBC
	{
		func_id.push_back(XC_GGA_X_RPBE);
        func_id.push_back(XC_GGA_C_PBE);
        func_type = 2;
        use_libxc = false;
	}
	else if ( xc_func == "WC") //WC+PBC
	{
        func_id.push_back(XC_GGA_X_WC);
        func_id.push_back(XC_GGA_C_PBE);
        func_type = 2;
        use_libxc = false;
	}	
	else if ( xc_func == "BLYP") //B88+LYP
	{
        func_id.push_back(XC_GGA_X_B88);
        func_id.push_back(XC_GGA_C_LYP);
        func_type = 2;
        use_libxc = false;
	}
	else if ( xc_func == "BP") //B88+P86
	{
        func_id.push_back(XC_GGA_X_B88);
        func_id.push_back(XC_GGA_C_P86);
        func_type = 2;
        use_libxc = false;
	} 
	else if ( xc_func == "PW91") //PW91_X+PW91_C
	{
        func_id.push_back(XC_GGA_X_PW91);
        func_id.push_back(XC_GGA_C_PW91);
        func_type = 2;
        use_libxc = false;
	} 
	else if ( xc_func == "HCTH") //HCTH_X+HCTH_C
	{
        func_id.push_back(XC_GGA_X_HCTH_A);
        func_id.push_back(XC_GGA_C_HCTH_A);
        func_type = 2;
        use_libxc = false;
	}
	else if ( xc_func == "OLYP") //OPTX+LYP
	{
        func_id.push_back(XC_GGA_X_OPTX);
        func_id.push_back(XC_GGA_C_LYP);
        func_type = 2;
        use_libxc = false;
	}
	else if ( xc_func == "SCAN")
	{
        func_id.push_back(XC_MGGA_X_SCAN);
        func_id.push_back(XC_MGGA_C_SCAN);
        func_type = 3;
	}
   	else if( xc_func == "PBE0")
	{
        func_id.push_back(XC_HYB_GGA_XC_PBEH);
        func_type = 4;
        use_libxc = false;
	}
    else if( xc_func == "HF" || xc_func == "OPT_ORB" ||  xc_func == "NONE")
    {
        // not doing anything
        if(xc_func == "HF") func_type = 4;
    }
    else if( xc_func == "HSE")
    {
        func_id.push_back(XC_HYB_GGA_XC_HSE06);
        func_type = 4;
    }
    else
    {
        std::cout << "functional name not recognized!" << std::endl;
    }

	if (func_id[0] == XC_GGA_X_OPTX)
	{
		std::cerr << "\n OPTX untested please test,";
	}

    if(func_type ==3 && GlobalV::BASIS_TYPE != "pw")
    {
        ModuleBase::WARNING_QUIT("set_xc_type","mGGA not realized for LCAO yet");
    }

#ifndef USE_LIBXC
    use_libxc = false;
#endif

}

void XC_Functional::xc_libxc(const double &rho, double &exc, double &vxc)
{
#ifdef USE_LIBXC
    double e,v;
    exc = vxc = 0.00;

	std::vector<xc_func_type> funcs = init_func(XC_UNPOLARIZED);

    for(xc_func_type &func : funcs)
    {
        if( func.info->family == XC_FAMILY_LDA)
        {
            // call Libxc function: xc_lda_exc_vxc
            xc_lda_exc_vxc( &func, 1, &rho, &e, &v);
        }
        exc += e;
        vxc += v;
    }
	return;
#endif 
}


void XC_Functional::xc(const double &rho, double &exc, double &vxc)
{

	double third = 1.0 / 3.0;
	double pi34 = 0.6203504908994e0 ; // pi34=(3/4pi)^(1/3)
	double rs;
    double e,v;
    
    exc = vxc = 0.00;
	
	rs = pi34 / std::pow(rho, third);

    for(int id : func_id)
    {
        switch( id )
        {
            // Exchange functionals containing slater exchange
            case XC_LDA_X: case XC_GGA_X_PBE: case XC_GGA_X_RPBE: 
            case XC_GGA_X_WC: case XC_GGA_X_B88: case XC_GGA_X_PW91:
            //  SLA,PBX,rPBX,PBXsol,WC,B88,PW91_X
                XC_Functional::slater(rs, e, v);break;

            // Exchange functionals containing attenuated slater exchange
            case XC_HYB_GGA_XC_PBEH:
            //  PBE0
                XC_Functional::slater(rs, e, v);
                e *= 0.75; v*= 0.75;
                break;

            // Correlation functionals containing PW correlation
            case XC_GGA_C_PBE: case XC_GGA_C_PW91: case XC_LDA_C_PW:
            //   PBC,PW91,PWLDA
                XC_Functional::pw(rs, 0, e, v);break;

            // Correlation functionals containing PZ correlation
            case XC_LDA_C_PZ: case XC_GGA_C_P86:
            //  PZ,P86
                XC_Functional::pz(rs, 0, e, v);break;

            // Correlation functionals containing LYP correlation
            case XC_GGA_C_LYP:
            //  BLYP
                XC_Functional::lyp(rs, e, v);break;

            // Functionals that are realized only using LIBXC
            case XC_HYB_GGA_XC_HSE06: case XC_MGGA_X_SCAN: case XC_MGGA_C_SCAN:
            //  HSE,SCAN_X,SCAN_C
                throw std::domain_error("functional unfinished in "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
            
            default:
                e = v = 0.0;
        }
        exc += e;
        vxc += v;
    }
	return;
}

void XC_Functional::xc_spin_libxc(const double &rhoup, const double &rhodw,
		double &exc, double &vxcup, double &vxcdw)
{
#ifdef USE_LIBXC
    double e, vup, vdw;
    double *rho_ud, *vxc_ud;
    exc = vxcup = vxcdw = 0.0;

    rho_ud = new double[2];
    vxc_ud = new double[2];
    rho_ud[0] = rhoup;
    rho_ud[1] = rhodw;

    std::vector<xc_func_type> funcs = init_func(XC_POLARIZED);

    for(xc_func_type &func : funcs)
    {
        if( func.info->family == XC_FAMILY_LDA)
        {
            // call Libxc function: xc_lda_exc_vxc
            xc_lda_exc_vxc( &func, 1, rho_ud, &e, vxc_ud);
        }
        exc += e;
        vxcup += vxc_ud[0];
        vxcdw += vxc_ud[1];
    }    


    delete[] rho_ud;
    delete[] vxc_ud;
#endif 
}

void XC_Functional::xc_spin(const double &rho, const double &zeta,
		double &exc, double &vxcup, double &vxcdw)
{
	static const double small = 1.e-10;
    double e, vup, vdw;
    exc = vxcup = vxcdw = 0.0;

	static const double third = 1.0 / 3.0;
	static const double pi34 = 0.62035049089940; 
	const double rs = pi34 / pow(rho, third);//wigner_sitz_radius;

    for(int id : func_id)
    {
        switch( id )
        {
            // Exchange functionals containing slater exchange
            case XC_LDA_X: case XC_GGA_X_PBE: case XC_GGA_X_RPBE: 
            case XC_GGA_X_WC: case XC_GGA_X_B88: case XC_GGA_X_PW91:
            //  SLA,PBX,rPBX,PBXsol,WC,B88,PW91_X
                XC_Functional::slater_spin(rho, zeta, e, vup, vdw);	break;
            
            // Exchange functionals containing attenuated slater exchange
            case XC_HYB_GGA_XC_PBEH:
            //  PBE0
                XC_Functional::slater_spin(rho, zeta, e, vup, vdw);
                e *= 0.75; vup *= 0.75; vdw *=0.75;
                break;

            // Correlation functionals containing PZ correlation
            case XC_LDA_C_PZ: case XC_GGA_C_P86:
            //  PZ,P86
                XC_Functional::pz_spin(rs, zeta, e, vup, vdw); break;

            // Correlation functionals containing PW correlationtests/integrate/101_PW_OU_pseudopot
            case XC_GGA_C_PBE:
            //   PBC,PBCsol
                XC_Functional::pw_spin(rs, zeta, e, vup, vdw); break;

            // Correlation functionals that already contains LDA components
            // so sets to 0 here
            case XC_GGA_X_HCTH_A: case XC_GGA_C_HCTH_A: case XC_GGA_X_OPTX:
            //HCTH_X  HTCH_C      OPTX
                e = vup = vdw =0.0;

            // Cases that are only realized in LIBXC
            default:
                throw std::domain_error("functional unfinished in "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));	break;

        }
        exc += e;
        vxcup += vup;
        vxcdw += vdw;
	}
	return;
}

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

void XC_Functional::gcxc(const double &rho, const double &grho, double &sxc,
          double &v1xc, double &v2xc)
{
    //-----------------------------------------------------------------------
    //     gradient corrections for exchange and correlation - Hartree a.u.
    //     exchange  :  Becke88
    //                  GGA (Generalized Gradient Approximation), PW91
    //                  PBE
    //                  revPBE
    //     correlation: Perdew86
    //                  GGA (PW91)
    //                  Lee-Yang-Parr
    //                  PBE
    //
    //     input:  rho, grho=|\nabla rho|^2
    //     definition:  E_x = \int E_x(rho,grho) dr
    //     output: sx = E_x(rho,grho)
    //             v1x= D(E_x)/D(rho)
    //             v2x= D(E_x)/D( D rho/D r_alpha ) / |\nabla rho|
    //             sc, v1c, v2c as above for correlation
    //
    // use funct
    // USE kinds
    // implicit none
    // real rho, grho, sx, sc, v1x, v2x, v1c, v2c;
    double small = 1.e-10;
    double s,v1,v2;
    sxc = v1xc = v2xc = 0.0;

    if (rho <= small)
    {
        return;
    }

    for(int id : func_id)
    {
        switch( id )
        {
            case 106: //B88
                XC_Functional::becke88(rho, grho, s, v1, v2);break;
            case 109: //PW91_X
                XC_Functional::ggax(rho, grho, s, v1, v2);break;
            case 101: //PBX
                XC_Functional::pbex(rho, grho, 0, s, v1, v2);break;
            case 117: //revised PBX
                XC_Functional::pbex(rho, grho, 1, s, v1, v2);break;
            case 34: //HCTH_X
                XC_Functional::hcth(rho, grho, s, v1, v2);break; //XC together
            case 97: //HCTH_C
                s = 0.0; v1 = 0.0; v2 = 0.0;break;
            case 110: //OPTX
                XC_Functional::optx(rho, grho, s, v1, v2);break;
            case 406: //PBE0
                XC_Functional::pbex(rho, grho, 0, s, v1, v2);
                s *= 0.75; v1 *= 0.75; v2 *= 0.75;break;
            case 116: //PBXsol
                XC_Functional::pbex(rho, grho, 2, s, v1, v2);break;
            case 118: //Wu-Cohen
                XC_Functional::wcx (rho, grho, s, v1, v2);break;
            case 132: //P86
                XC_Functional::perdew86(rho, grho, s, v1, v2);break;
            case 134: //PW91_C
                XC_Functional::ggac(rho, grho, s, v1, v2);break;
            case 130: //PBC
                XC_Functional::pbec(rho, grho, 0, s, v1, v2);break;
            case 133: //PBCsol
                XC_Functional::pbec(rho, grho, 1, s, v1, v2);break;
            case 131: //BLYP
                XC_Functional::glyp(rho, grho, s, v1, v2); break;
            default: //SCAN_X,SCAN_C,HSE, and so on
                throw std::domain_error("functional unfinished in "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
        }
        sxc += s;
        v1xc += v1;
        v2xc += v2;
    }

    return;
}

void XC_Functional::gcxc_libxc(const double &rho, const double &grho, double &sxc,
          double &v1xc, double &v2xc)
{
#ifdef USE_LIBXC
    double s,v1,v2;
    sxc = v1xc = v2xc = 0.0;

	std::vector<xc_func_type> funcs = init_func(XC_UNPOLARIZED);

    for(xc_func_type &func : funcs)
    {
        if( func.info->family == XC_FAMILY_GGA || func.info->family == XC_FAMILY_HYB_GGA)
        {
            // call Libxc function: xc_gga_exc_vxc
            xc_gga_exc_vxc( &func, 1, &rho, &grho, &s, &v1, &v2);
        }
        sxc += s;
        v1xc += v1;
        v2xc += v2;
    }
	return;
#endif
}

void XC_Functional::gcxc_spin_libxc(double rhoup, double rhodw, 
        ModuleBase::Vector3<double> gdr1, ModuleBase::Vector3<double> gdr2,
        double &sxc, double &v1xcup, double &v1xcdw, double &v2xcup, double &v2xcdw, double &v2xcud)
{
#ifdef USE_LIBXC
	std::vector<xc_func_type> funcs = init_func(XC_POLARIZED);
    double *rho, *grho, *v1xc, *v2xc, *sgn, s;
    sxc = v1xcup = v1xcdw = 0.0;
    v2xcup = v2xcdw = v2xcud = 0.0;
    rho = new double[2];
    grho= new double[3];
    v1xc= new double[2];
    v2xc= new double[3];
    sgn = new double[2];
    
    rho[0] = rhoup;
    rho[1] = rhodw;
    grho[0] = gdr1.norm2();
    grho[1] = gdr1 * gdr2;
    grho[2] = gdr2.norm2();

    const double rho_threshold = 1E-6;
    const double grho_threshold = 1E-10;

    for(xc_func_type &func : funcs)
    {
        if( func.info->family == XC_FAMILY_GGA || func.info->family == XC_FAMILY_HYB_GGA)
        {
            sgn[0] = sgn[1] = 1.0;
            // call Libxc function: xc_gga_exc_vxc
            xc_gga_exc_vxc( &func, 1, rho, grho, &s, v1xc, v2xc);
            if(func.info->kind==XC_CORRELATION)
            {
                if ( rho[0]<rho_threshold || sqrt(abs(grho[0]))<grho_threshold )
                    sgn[0] = 0.0;
                if ( rho[1]<rho_threshold || sqrt(abs(grho[2]))<grho_threshold )
                    sgn[1] = 0.0;
            }
            sxc += s * (rho[0] * sgn[0] + rho[1] * sgn[1]);
            v1xcup += v1xc[0] * sgn[0];
            v1xcdw += v1xc[1] * sgn[1];
            v2xcup += 2.0 * v2xc[0] * sgn[0];
            v2xcud += v2xc[1] * sgn[0] * sgn[1];
            v2xcdw += 2.0 * v2xc[2] * sgn[1];
        }
    }
    delete[] grho;
    delete[] rho;
    delete[] v1xc;
    delete[] v2xc;
    delete[] sgn;

    return;
#endif
} 

//-----------------------------------------------------------------------
void XC_Functional::gcx_spin(double rhoup, double rhodw, double grhoup2, double grhodw2,
              double &sx, double &v1xup, double &v1xdw, double &v2xup, double &v2xdw)
{
    //--------------------------------------------------------------------
    //     gradient corrections for exchange - Hartree a.u.
    //     Implemented:  Becke88, GGA (PW91), PBE, revPBE
    //
    // use funct
    // USE kinds
    // implicit none

    //     dummy arguments

    // real rhoup, rhodw, grhoup2, grhodw2, sx, v1xup, v1xdw,
    //	v2xup, v2xdw;
    // up and down charge
    // up and down gradient of the charge
    // exchange and correlation energies
    // derivatives of exchange wr. rho
    // derivatives of exchange wr. grho

    // parameter :
    double small = 1.e-10;
    double sxup, sxdw;
    int iflag;

    // exchange
    double rho = rhoup + rhodw;

    sx = 0.00;
    v1xup = 0.00; v2xup = 0.00;
    v1xdw = 0.00; v2xdw = 0.00;

    if (rho <= small)
    {
        return;
    }

    // not the correct way to do things, will change later
    // should put exchange and correlation together
    // like the others
    //for(int id : func_id)
    //{
        int id = func_id[0];
        switch( id )
        {
            case 106: //B88
                if (rhoup > small && sqrt(fabs(grhoup2)) > small)
                {
                    XC_Functional::becke88_spin(rhoup, grhoup2, sxup, v1xup, v2xup);
                }
                if (rhodw > small && sqrt(fabs(grhodw2)) > small)
                {
                    XC_Functional::becke88_spin(rhodw, grhodw2, sxdw, v1xdw, v2xdw);
                }
                break;
            case 101: //PBX
                if (rhoup > small && sqrt(fabs(grhoup2)) > small)
                {
                    XC_Functional::pbex(2.0 * rhoup, 4.0 * grhoup2, 0, sxup, v1xup, v2xup);
                }
                if (rhodw > small && sqrt(fabs(grhodw2)) > small)
                {
                    XC_Functional::pbex(2.0 * rhodw, 4.0 * grhodw2, 0, sxdw, v1xdw, v2xdw);
                }
                break;
            case 117: //revised PBX
                if (rhoup > small && sqrt(fabs(grhoup2)) > small)
                {
                    XC_Functional::pbex(2.0 * rhoup, 4.0 * grhoup2, 1, sxup, v1xup, v2xup);
                }
                if (rhodw > small && sqrt(fabs(grhodw2)) > small)
                {
                    XC_Functional::pbex(2.0 * rhodw, 4.0 * grhodw2, 1, sxdw, v1xdw, v2xdw);
                }
                break;
            case 406: //PBE0
                if (rhoup > small && sqrt(fabs(grhoup2)) > small)
                {
                    XC_Functional::pbex(2.0 * rhoup, 4.0 * grhoup2, 0, sxup, v1xup, v2xup);
                    sxup *= 0.75; v1xup *= 0.75; v2xup *= 0.75;
                }
                if (rhodw > small && sqrt(fabs(grhodw2)) > small)
                {
                    XC_Functional::pbex(2.0 * rhodw, 4.0 * grhodw2, 0, sxdw, v1xdw, v2xdw);
        	    	sxdw *= 0.75; v1xdw *= 0.75; v2xdw *= 0.75;            
                }
                break;
            case 116: //PBXsol
                if (rhoup > small && sqrt(fabs(grhoup2)) > small)
                {
                    XC_Functional::pbex(2.0 * rhoup, 4.0 * grhoup2, 2, sxup, v1xup, v2xup);
                }
                if (rhodw > small && sqrt(fabs(grhodw2)) > small)
                {
                    XC_Functional::pbex(2.0 * rhodw, 4.0 * grhodw2, 2, sxdw, v1xdw, v2xdw);
                }
                break;
            default:
                sxup = 0.0; sxdw = 0.0;
                v1xup = 0.0; v2xup = 0.0;
                v1xdw = 0.0; v2xdw = 0.0;
        }
        sx = 0.50 * (sxup + sxdw);
        v2xup = 2.0 * v2xup;
        v2xdw = 2.0 * v2xdw;       
    //}

    return;
} //end subroutine gcx_spin

//
//-----------------------------------------------------------------------
void XC_Functional::gcc_spin(double rho, double &zeta, double grho, double &sc,
              double &v1cup, double &v1cdw, double &v2c)
{
    //-------------------------------------------------------------------
    //     gradient corrections for correlations - Hartree a.u.
    //     Implemented:  Perdew86, GGA (PW91), PBE

    // use funct
    // USE kinds
    // implicit none

    //     dummy arguments

    // real(kind=DP) :: rho, zeta, grho, sc, v1cup, v1cdw, v2c
    // the total charge
    // the magnetization
    // the gradient of the charge squared
    // exchange and correlation energies
    // derivatives of correlation wr. rho
    // derivatives of correlation wr. grho

    // parameter :
    double small = 1.0e-10;
	double epsr = 1.0e-6;

    double x;

    sc = 0.00;
    v1cup = 0.00; v1cdw = 0.00;
    v2c = 0.00;
    if (abs(zeta) - 1.0 > small || rho <= small || sqrt(abs(grho)) <= small)
    {
        return;
    }
    else
    {
        // ... ( - 1.0 + epsr )  <  zeta  <  ( 1.0 - epsr )
        // zeta = SIGN( MIN( ABS( zeta ), ( 1.D0 - epsr ) ) , zeta )
        x = std::min(abs(zeta), (1.0 - epsr));
		if(zeta>0)
		{
			zeta = x;
		}
		else
		{
			zeta = -x;
		}
    } //endif

    //for(int id : func_id)
    //{
        int id = func_id[1];
        switch( id )
        {          
            case 132: //P86
                XC_Functional::perdew86_spin(rho, zeta, grho, sc, v1cup, v1cdw, v2c);break;
            case 134: //PW91_C
                XC_Functional::ggac_spin(rho, zeta, grho, sc, v1cup, v1cdw, v2c);break;
            case 130: //PBC
                XC_Functional::pbec_spin(rho, zeta, grho, 1, sc, v1cup, v1cdw, v2c);break;
            case 133: //PBCsol
                XC_Functional::pbec_spin(rho, zeta, grho, 2, sc, v1cup, v1cdw, v2c);break;
        }
    //}
    return;
} //end subroutine gcc_spin

#ifdef USE_LIBXC
std::vector<xc_func_type> XC_Functional::init_func(const int xc_polarized)
{
	// 'funcs' is the return value
	std::vector<xc_func_type> funcs;

	//-------------------------------------------
	// define a function named 'add_func', which 
	// initialize a functional according to its ID
	//-------------------------------------------
	auto add_func = [&]( const int func_id )
	{
		funcs.push_back({});
		// 'xc_func_init' is defined in Libxc
		xc_func_init( &funcs.back(), func_id, xc_polarized );
	};

	for(int id : func_id)
	{
		if( id == 406 ) // PBE0
		{
			add_func( XC_HYB_GGA_XC_PBEH );		
			double parameter_hse[3] = { GlobalC::exx_global.info.hybrid_alpha, 
				GlobalC::exx_global.info.hse_omega, 
				GlobalC::exx_global.info.hse_omega };
			xc_func_set_ext_params(&funcs.back(), parameter_hse);	
		}
		else if( id == 428 ) // HSE06 hybrid functional
		{
			add_func( XC_HYB_GGA_XC_HSE06 );	
			double parameter_hse[3] = { GlobalC::exx_global.info.hybrid_alpha, 
				GlobalC::exx_global.info.hse_omega, 
				GlobalC::exx_global.info.hse_omega };
			xc_func_set_ext_params(&funcs.back(), parameter_hse);
		}
		else
		{
			add_func( id );
		}
	}
	return funcs;
}
#endif