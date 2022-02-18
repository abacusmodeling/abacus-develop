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
	else if ( xc_func == "PBEsol") //PBXsol+PBCsol
	{
        func_id.push_back(XC_GGA_X_PBE_SOL);
        func_id.push_back(XC_GGA_C_PBE_SOL);
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
		GlobalV::DFT_META = 1;
	}
   	else if( xc_func == "PBE0")
	{
        func_id.push_back(XC_HYB_GGA_XC_PBEH);
        func_type = 4;
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

	if (func_id[0] == 110)
	{
		std::cerr << "\n OPTX untested please test,";
	}
}

void XC_Functional::xc_libxc(const double &rho, double &exc, double &vxc)
{

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
            case XC_LDA_X: case XC_GGA_X_PBE: case XC_GGA_X_RPBE: case XC_GGA_X_PBE_SOL: 
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
            case XC_GGA_C_PBE: case XC_GGA_C_PBE_SOL: case XC_GGA_C_PW91: case XC_LDA_C_PW:
            //   PBC,PBCsol
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
            case XC_LDA_X: case XC_GGA_X_PBE: case XC_GGA_X_RPBE: case XC_GGA_X_PBE_SOL: 
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
            case XC_GGA_C_PBE: case XC_GGA_C_PBE_SOL:
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

//Slater exchange with alpha=2/3
void XC_Functional::slater(const double &rs, double &ex, double &vx)
{
	// f = -9/8*(3/2pi)^(2/3)
	const double f = -0.687247939924714e0;
	const double alpha = 2.00 / 3.00;
	ex = f * alpha / rs;
	vx = 4.0 / 3.0 * f * alpha / rs;
	return;
}

// Slater exchange with alpha=2/3, spin-polarized case
void XC_Functional::slater_spin( const double &rho, const double &zeta, 
		double &ex, double &vxup, double &vxdw)
{
    const double f = - 1.107838149573033610;
    const double alpha = 2.00 / 3.00;
    // f = -9/8*(3/pi)^(1/3)
    const double third = 1.0 / 3.0;
    const double p43 = 4.0 / 3.0;

    double rho13 = pow(((1.0 + zeta) * rho) , third);
    double exup = f * alpha * rho13;
    vxup = p43 * f * alpha * rho13;
    rho13 = pow(((1.0 - zeta) * rho) , third);
    double exdw = f * alpha * rho13;
    vxdw = p43 * f * alpha * rho13;
    ex = 0.50 * ((1.0 + zeta) * exup + (1.0 - zeta) * exdw);

    return;
}


//Slater exchange with alpha=1, corresponding to -1.374/r_s Ry
//used to recover old results
void XC_Functional::slater1(const double &rs, double &ex, double &vx)
{
	const double f = -0.687247939924714e0;
	const double alpha = 1.0;
	ex = f * alpha / rs;
	vx = 4.0 / 3.0 * f * alpha / rs;
	return;
}


// Slater exchange with alpha=2/3, spin-polarized case
void XC_Functional::slater1_spin( const double &rho, const double &zeta, double &ex, double &vxup, double &vxdw)
{
    const double f = - 1.107838149573033610;
	const double alpha = 1.00;
	const double third = 1.0 / 3.0;
	const double p43 = 4.0 / 3.0;
    // f = -9/8*(3/pi)^(1/3)

    double rho13 = pow(((1.0 + zeta) * rho) , third);
    double exup = f * alpha * rho13;
    vxup = p43 * f * alpha * rho13;
    rho13 = pow(((1.0 - zeta) * rho) , third);
    double exdw = f * alpha * rho13;
    vxdw = p43 * f * alpha * rho13;
    ex = 0.50 * ((1.0 + zeta) * exup + (1.0 - zeta) * exdw);

    return;
} // end subroutine slater1_spin


// Slater exchange with alpha=2/3 and Relativistic exchange
void XC_Functional::slater_rxc(const double &rs, double &ex, double &vx)
{
	static const double pi = 3.14159265358979;
    const double trd = 1.0 / 3.0;
    //const double ftrd = 4.0 / 3.0;
    //const double tftm = pow(2.0, ftrd) - 2.0;
    const double a0 = pow((4.0 / (9.0 * pi)), trd);
    // X-alpha parameter:
    const double alp = 2 * trd;

    double vxp = -3 * alp / (2 * pi * a0 * rs);
    double exp = 3 * vxp / 4;
    const double beta = 0.014 / rs;
    const double sb = sqrt(1 + beta * beta);
    const double alb = log(beta + sb);
    vxp = vxp * (-0.5 + 1.5 * alb / (beta * sb));
    double x = (beta * sb - alb) / (beta * beta);
    exp = exp * (1.0 - 1.5 * x * x);
    vx = vxp;
    ex = exp;
	return;	
}

// Slater exchange with alpha=2/3, relativistic exchange case
void XC_Functional::slater_rxc_spin( const double &rho, const double &z, 
		double &ex, double &vxup, double &vxdw)
{
    if (rho <= 0.0)
    {
        ex = vxup = vxdw = 0.0;
        return;
    }

	const double pi = 3.141592653589790;

	const double trd = 1.0 / 3.0;
	const double ftrd = 4.0 / 3.0;
	double tftm = pow(2.0, ftrd) - 2;
	double a0 = pow((4 / (9 * pi)), trd);

	double alp = 2 * trd;

	double fz = (pow((1 + z), ftrd) + pow((1 - z), ftrd) - 2) / tftm;
	double fzp = ftrd * (pow((1 + z), trd) - pow((1 - z), trd)) / tftm;

    double rs = pow((3 / (4 * pi * rho)), trd);
    double vxp = -3 * alp / (2 * pi * a0 * rs);
    double exp = 3 * vxp / 4;

    double beta = 0.014 / rs;
    double sb = sqrt(1 + beta * beta);
    double alb = log(beta + sb);
    vxp = vxp * (-0.5 + 1.5 * alb / (beta * sb));
    exp = exp * (1.0- 1.5*( (beta*sb-alb) / (beta*beta) )
			        * 1.5*( (beta*sb-alb) / (beta*beta) ));

    double x = (beta * sb - alb) / (beta * beta);

    exp = exp * (1.0 - 1.5 * x * x);

    double vxf = pow(2.0, trd) * vxp;

    double exf = pow(2.0, trd) * exp;

    vxup  = vxp + fz * (vxf - vxp) + (1 - z) * fzp * (exf - exp);

    vxdw  = vxp + fz * (vxf - vxp) - (1 + z) * fzp * (exf - exp);

    ex    = exp + fz * (exf - exp);

    return;
}

//LDA parameterization form Monte Carlo data
//iflag=0: J.P. Perdew and A. Zunger, PRB 23, 5048 (1981)
//iflag=1: G. Ortiz and P. Ballone, PRB 50, 1391 (1994)
void XC_Functional::pz(const double &rs, const int &iflag, double &ec, double &vc)
{
//	ModuleBase::TITLE("XC_Functional","pz");
	const double a[2] = {0.0311, 0.031091};
    const double b[2] = { -0.048, -0.046644 };
    const double c[2] = { 0.0020, 0.00419 };
    const double d[2] = { -0.0116, -0.00983 };
    const double gc[2] = { -0.1423, -0.103756};
    const double b1[2] = { 1.0529, 0.56371 };
    const double b2[2] = { 0.3334, 0.27358 };

    if (rs < 1.0)
    {
        // high density formula
        const double lnrs = log(rs);
        ec = a [iflag] * lnrs + b [iflag] + c [iflag] * rs * lnrs + d [iflag] * rs;
        vc = a [iflag] * lnrs + (b [iflag] - a [iflag] / 3.0) + 2.0 /
             3.0 * c [iflag] * rs * lnrs + (2.0 * d [iflag] - c [iflag])
             / 3.0 * rs;
    }
    else
    {
        // interpolation formula
        const double rs12 = sqrt(rs);
        const double ox = 1.0 + b1 [iflag] * rs12 + b2 [iflag] * rs;
        const double dox = 1.0 + 7.0 / 6.0 * b1 [iflag] * rs12 + 4.0 / 3.0 *
              b2 [iflag] * rs;
        ec = gc [iflag] / ox;
        vc = ec * dox / ox;
    }
    return;
}

// J.P. Perdew and Y. Wang, PRB 45, 13244 (1992)
void XC_Functional::pz_spin( const double &rs, const double &zeta, 
	double &ec, double &vcup, double &vcdw)
{
//	ModuleBase::TITLE("XC_Functional","pz_spin");
	double ecu, vcu, ecp, vcp ;
    static const double p43 = 4.00 / 3.0;
    static const double third = 1.0 / 3.0;

    // unpolarized part (Perdew-Zunger formula)
    XC_Functional::pz(rs, 0, ecu, vcu);

    // polarization contribution
	XC_Functional::pz_polarized(rs, ecp, vcp);

    double fz = (pow((1.00 + zeta) , p43) + pow((1.0 - zeta) , p43) - 2.0) /
         (pow(2.0 , p43) - 2.0);
    double dfz = p43 * (pow((1.00 + zeta) , third) - pow((1.0 - zeta) , third))
          / (pow(2.0, p43) - 2.0);

	//correlation energy
    ec = ecu + fz * (ecp - ecu);
	
	// spin up correlation potential
    vcup = vcu + fz * (vcp - vcu) + (ecp - ecu) * dfz * (1.0 - zeta);
	// spin down correlation potential
    vcdw = vcu + fz * (vcp - vcu) + (ecp - ecu) * dfz * (- 1.0 - zeta);

    return;
}

// J.P. Perdew and A. Zunger, PRB 23, 5048 (1981)
void XC_Functional::pz_polarized( const double &rs, double &ec, double &vc)
{
//	ModuleBase::TITLE("XC_Functional","pz_polarized");
	static const double a = 0.015550;
	static const double b = -0.02690;
	static const double c = 0.00070;
	static const double d = -0.00480;
	static const double gc = -0.08430;
	static const double b1 = 1.39810;
	static const double b2 = 0.26110;
	
	if (rs < 1.00)
	{
		double lnrs = log(rs);
		ec = a * lnrs + b + c * rs * lnrs + d * rs;
		vc = a * lnrs + (b - a / 3.0) + 2.0 / 3.0 * c * rs * lnrs +
			(2.0 * d - c) / 3.0 * rs;
	}
	else
	{
		// interpolation formula
		double rs12 = sqrt(rs);
		double ox = 1.0 + b1 * rs12 + b2 * rs;
		double dox = 1.0 + 7.0 / 6.0 * b1 * rs12 + 4.0 / 3.0 * b2 * rs;
		ec = gc / ox;
		vc = ec * dox / ox;
	}
	return;
}

// S.H. Vosko, L. Wilk, and M. Nusair, Can. J. Phys. 58, 1200 (1980)
void XC_Functional::vwn(const double &rs, double &ec, double &vc)
{
    const double a = 0.0310907;
    const double b = 3.72744;
    const double c = 12.9352;
    const double x0 = -0.10498;
    double q, f1, f2, f3, rs12, fx, qx, tx, tt;

    q = sqrt(4.0 * c - b * b);
    f1 = 2.0 * b / q;
    f2 = b * x0 / (x0 * x0 + b * x0 + c);
    f3 = 2.0 * (2.0 * x0 + b) / q;
    rs12 = sqrt(rs);
    fx = rs + b * rs12 + c;
    qx = atan(q / (2.0 * rs12 + b));

    double x;
    x = (rs12 - x0);
    ec = a * (log(rs / fx) + f1 * qx - f2 * (log((x * x) /
              fx) + f3 * qx));
    tx = 2.0 * rs12 + b;
    tt = tx * tx + q * q;

    vc = ec - rs12 * a / 6.0 * (2.0 / rs12 - tx / fx - 4.0 * b /
	tt - f2 * (2.0 / (rs12 - x0) - tx / fx - 4.0 * (2.0 * x0 + b)/ tt));

    return;
}

// C. Lee, W. Yang, and R.G. Parr, PRB 37, 785 (1988)
// LDA part only
void XC_Functional::lyp(const double &rs, double &ec, double &vc)
{
    const double a = 0.04918e0;
    const double b = 0.1320 * 2.87123400018819108e0;
    // pi43 = (4pi/3)^(1/3)
    const double pi43 = 1.61199195401647e0;
    const double c = 0.2533e0 * pi43;
    const double d = 0.349e0 * pi43;
    double ecrs, ox;

    ecrs = b * exp(- c * rs);
    ox = 1.0 / (1.0 + d * rs);
    ec = - a * ox * (1.0 + ecrs);
    vc = ec - rs / 3.0 * a * ox * (d * ox + ecrs * (d * ox + c));

    return;
}

void XC_Functional::pw(const double &rs, const int &iflag, double &ec, double &vc)
{
    const double a = 0.0310910;
    const double b1 = 7.59570;
    const double b2 = 3.58760;
    const double c0 = a;
    const double c1 = 0.0466440;
    const double c2 = 0.006640;
    const double c3 = 0.010430;
    const double d0 = 0.43350;
    const double d1 = 1.44080;
    double lnrs, rs12, rs32, rs2, om, dom, olog;
    double a1[2] = { 0.213700, 0.0264810 };
	double b3[2] = { 1.63820, -0.466470 };
	double b4[2] = { 0.492940, 0.133540 };

	// high- and low-density formulae implemented but not used in PW case
	// (reason: inconsistencies in PBE/PW91 functionals)
	
   	if (rs < 1 && iflag == 1)
    {
        // high density formula
        lnrs = log(rs);
        ec = c0 * lnrs - c1 + c2 * rs * lnrs - c3 * rs;
        vc = c0 * lnrs - (c1 + c0 / 3.0) + 2.0 / 3.0 * c2 * rs *
             lnrs - (2.0 * c3 + c2) / 3.0 * rs;
    }
    else if (rs > 100.0 && iflag == 1)
    {
        // low density formula
        ec = - d0 / rs + d1 / pow(rs, 1.50);
        vc = - 4.0 / 3.0 * d0 / rs + 1.50 * d1 / pow(rs, 1.50);
    }
    else
    {
        // interpolation formula
        rs12 = sqrt(rs);
        rs32 = rs * rs12;
        rs2 = rs * rs;
        om = 2.0 * a * (b1 * rs12 + b2 * rs + b3 [iflag] * rs32 + b4 [
                            iflag] * rs2);
        dom = 2.0 * a * (0.50 * b1 * rs12 + b2 * rs + 1.50 * b3 [
                             iflag] * rs32 + 2.0 * b4 [iflag] * rs2);
        olog = log(1.0 + 1.0 / om);
        ec = - 2.0 * a * (1.0 + a1 [iflag] * rs) * olog;
        vc = - 2.0 * a * (1.0 + 2.0 / 3.0 * a1 [iflag] * rs)
             * olog - 2.0 / 3.0 * a * (1.0 + a1 [iflag] * rs) * dom /
             (om * (om + 1.0));
    }
	return;
}

// J.P. Perdew and Y. Wang, PRB 45, 13244 (1992)
void XC_Functional::pw_spin( const double &rs, const double &zeta, 
		double &ec, double &vcup, double &vcdw)
{
	const double a = 0.0310910;
    const double a1 = 0.213700;
    const double b1 = 7.59570;
    const double b2 = 3.58760;
    const double b3 = 1.63820;
    const double b4 = 0.492940;
    //const double c0 = a;
    //const double c1 = 0.0466440;
    //const double c2 = 0.006640;
    //const double c3 = 0.010430;
    //const double d0 = 0.43350;
    //const double d1 = 1.44080;
	// polarized parameters
    const double ap = 0.0155450;
    const double a1p = 0.205480;
    const double b1p = 14.11890;
    const double b2p = 6.19770;
    const double b3p = 3.36620;
    const double b4p = 0.625170;
    //const double c0p = ap;
    //const double c1p = 0.0255990;
    //const double c2p = 0.003190;
    //const double c3p = 0.003840;
    //const double d0p = 0.32870;
    //const double d1p = 1.76970;
	// antiferro
    const double aa = 0.0168870;
    const double a1a = 0.111250;
    const double b1a = 10.3570;
    const double b2a = 3.62310;
    const double b3a = 0.880260;
    const double b4a = 0.496710;
    //const double c0a = aa;
    //const double c1a = 0.0354750;
    //const double c2a = 0.001880;
    //const double c3a = 0.005210;
    //const double d0a = 0.22400;
    //const double d1a = 0.39690;

	const double fz0 = 1.7099210;
	double rs12, rs32, rs2, zeta2, zeta3, zeta4, fz, dfz;
	double om, dom, olog, epwc, vpwc;
    double omp, domp, ologp, epwcp, vpwcp;
    double oma, doma, ologa, alpha, vpwca;

	//     if(rs.lt.0.5d0) then
    // high density formula (not implemented)
    //
    //     else if(rs.gt.100.d0) then
    // low density formula  (not implemented)
    //
    //     else
    // interpolation formula
    zeta2 = zeta * zeta;
    zeta3 = zeta2 * zeta;
    zeta4 = zeta3 * zeta;
    rs12 = sqrt(rs);
    rs32 = rs * rs12;
    rs2 = rs * rs;
    // unpolarised
    om = 2.0 * a * (b1 * rs12 + b2 * rs + b3 * rs32 + b4 * rs2);
    dom = 2.0 * a * (0.50 * b1 * rs12 + b2 * rs + 1.50 * b3 * rs32
                     + 2.0 * b4 * rs2);
    olog = log(1.0 + 1.00 / om);
    epwc = - 2.0 * a * (1.0 + a1 * rs) * olog;
    vpwc = - 2.0 * a * (1.0 + 2.0 / 3.0 * a1 * rs) * olog - 2.0 /
           3.0 * a * (1.0 + a1 * rs) * dom / (om * (om + 1.0));
    // polarized
    omp = 2.0 * ap * (b1p * rs12 + b2p * rs + b3p * rs32 + b4p * rs2);
    domp = 2.0 * ap * (0.50 * b1p * rs12 + b2p * rs + 1.50 * b3p *
                       rs32 + 2.0 * b4p * rs2);
    ologp = log(1.0 + 1.00 / omp);
    epwcp = - 2.0 * ap * (1.0 + a1p * rs) * ologp;
    vpwcp = - 2.0 * ap * (1.0 + 2.0 / 3.0 * a1p * rs) * ologp -
            2.0 / 3.0 * ap * (1.0 + a1p * rs) * domp / (omp * (omp + 1.0));
    // antiferro
    oma = 2.0 * aa * (b1a * rs12 + b2a * rs + b3a * rs32 + b4a * rs2);
    doma = 2.0 * aa * (0.50 * b1a * rs12 + b2a * rs + 1.50 * b3a *
                       rs32 + 2.0 * b4a * rs2);
    ologa = log(1.0 + 1.00 / oma);
    alpha = 2.0 * aa * (1.0 + a1a * rs) * ologa;
    vpwca = + 2.0 * aa * (1.0 + 2.0 / 3.0 * a1a * rs) * ologa +
            2.0 / 3.0 * aa * (1.0 + a1a * rs) * doma / (oma * (oma + 1.0));

	fz = (pow((1.0 + zeta) , (4.0 / 3.0)) + pow((1.0 - zeta) , (4.0 /
            3.0)) - 2.0) / (pow(2.0, (4.0 / 3.0)) - 2.0);
    dfz = (pow((1.0 + zeta) , (1.0 / 3.0)) - pow((1.0 - zeta) , (1.0 /
            3.0))) * 4.0 / (3.0 * (pow(2.0 , (4.0 / 3.0)) - 2.0));

    ec = epwc + alpha * fz * (1.0 - zeta4) / fz0 + (epwcp - epwc)
         * fz * zeta4;

    vcup = vpwc + vpwca * fz * (1.0 - zeta4) / fz0 
	+ (vpwcp - vpwc) * fz * zeta4 
	+ (alpha / fz0 * (dfz * (1.0 - zeta4) - 4.0 * fz * zeta3) 
	+ (epwcp - epwc) * (dfz * zeta4 + 4.0 * fz * zeta3)) * (1.0 - zeta);

    vcdw = vpwc + vpwca * fz * (1.0 - zeta4) / fz0 
	+ (vpwcp - vpwc) * fz * zeta4 
	- (alpha / fz0 * (dfz * (1.0 - zeta4) - 4.0 * fz * zeta3) 
	+ (epwcp - epwc) * (dfz * zeta4 + 4.0 * fz * zeta3)) * (1.0 + zeta);
	return;
}

void XC_Functional::wigner( const double &rs, double &ec, double &vc)
{
    const double pi34 = 0.6203504908994e0;
    // pi34=(3/4pi)^(1/3), rho13=rho^(1/3)

	double rho13 = pi34 / rs;
    double x = 1.0 + 12.570 * rho13;
    vc = - rho13 * ((0.9436560 + 8.89630 * rho13) / (x * x));
    ec = - 0.7380 * rho13 * (0.9590 / (1.0 + 12.570 * rho13));
    return;
}

// L. Hedin and  B.I. Lundqvist,  J. Phys. C 4, 2064 (1971)
void XC_Functional::hl( const double &rs, double &ec, double &vc)
{
    double a, x;
    a = log(1.00 + 21.0 / rs);
    x = rs / 21.00;
    ec = a + (pow(x, 3) * a - x * x) + x / 2.0 - 1.00 / 3.00;
    ec = - 0.02250 * ec;
    vc = - 0.02250 * a;
    return;
}

// O. Gunnarsson and B. I. Lundqvist, PRB 13, 4274 (1976)
void XC_Functional::gl( const double &rs, double &ec, double &vc)
{
    const double c = 0.0333;
    const double r = 11.4;
    // c=0.0203, r=15.9 for the paramagnetic case
    double x = rs / r;
    vc = - c * log(1.0 + 1.0 / x);
    ec = - c * ((1.0 + pow(x, 3)) * log(1.0 + 1.0 / x) - 1.00 /
                3.00 + x * (0.50 - x));
	return;	
}


void XC_Functional::pbex(const double &rho, const double &grho, const int &iflag, 
double &sx, double &v1x, double &v2x)
{
    // PBE exchange (without Slater exchange):
    // iflag=0  J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)
    // iflag=1  "revised' PBE: Y. Zhang et al., PRL 80, 890 (1998)

    // input: charge and squared gradient
    // output: energy
    // output: potential
    // (3*pi2*|rho|)^(1/3)
    // |grho|
    // |grho|/(2*kf*|rho|)
    // s^2
    // n*ds/dn
    // n*ds/d(gn)
    // exchange energy LDA part
    // exchange energy gradient part
    
	// numerical coefficients (NB: c2=(3 pi^2)^(1/3) )
    const double third = 1.0 / 3.0;
    const double pi = 3.14159265358979323846;
    const double c1 = 0.750 / pi;
    const double c2 = 3.0936677262801360;
    const double c5 = 4.0 * third;
    // parameters of the functional
    double k[3] = { 0.8040, 1.24500, 0.8040 };
    const double mu[3] = {0.2195149727645171, 0.2195149727645171, 0.12345679012345679} ;//modified by zhengdy, to ensure the same parameters with another dft code.

    const double agrho = sqrt(grho);
    const double kf = c2 * pow(rho, third);
    const double dsg = 0.50 / kf;
    const double s1 = agrho * dsg / rho;
    const double s2 = s1 * s1;
    const double ds = - c5 * s1;

    // Energy
    const double f1 = s2 * mu[iflag] / k [iflag];
    const double f2 = 1.0 + f1;
    const double f3 = k [iflag] / f2;
    const double fx = k [iflag] - f3;
    const double exunif = - c1 * kf;
    sx = exunif * fx;

    // Potential
    const double dxunif = exunif * third;
    const double dfx1 = f2 * f2;
    const double dfx = 2.0 * mu[iflag] * s1 / dfx1;
    
	v1x = sx + dxunif * fx + exunif * dfx * ds;
    v2x = exunif * dfx * dsg / agrho;
    sx = sx * rho;

	return;
}

void XC_Functional::pwcorr(const double r, const double c[], double &g, double &dg)
{
    // USE kinds
    // implicit none
    // real(kind=DP) :: r, g, dg, c(6)
    double r12, r32, r2, rb, drb, sb;

    r12 = sqrt(r);
    r32 = r * r12;
    r2 = r * r;
    rb = c[3] * r12 + c[4] * r + c[5] * r32 + c[6] * r2;	//c[i] i=0--5;
    sb = 1.00 + 1.00 / (2.00 * c[1] * rb);
    g = -2.00 * c[1] * (1.00 + c[2] * r) * log(sb);
    drb = c[3] / (2.00 * r12) + c[4] + 1.50 * c[5] * r12 + 2.00 * c[6] * r;
    dg = (1.00 + c[2] * r) * drb / (rb * rb * sb) - 2.00 * c[1] * c[2] * log(sb);

    return;
} //end subroutine pwcorr

void XC_Functional::hcth(const double rho, const double grho, double &sx, double &v1x, double &v2x)
{
    //     ===============================================================
    //     HCTH/120, JCP 109, p. 6264 (1998)
    //     Parameters set-up after N.L. Doltsisnis & M. Sprik (1999)
    //     Present release: Mauro Boero, Tsukuba, 11/05/2004
    //--------------------------------------------------------------------------
    //     rhoa = rhob = 0.5 * rho
    //     grho is the SQUARE of the gradient of rho// --> gr=sqrt(grho)
    //     sx  : total exchange correlation energy at point r
    //     v1x : d(sx)/drho  (eq. dfdra = dfdrb in original)
    //     v2x : 1/gr*d(sx)/d(gr) (eq. 0.5 * dfdza = 0.5 * dfdzb in original)
    //--------------------------------------------------------------------------
    // USE kinds
    // implicit none
    // real(kind=DP) :: rho, grho, sx, v1x, v2x

    // parameter :
    double o3 = 1.00 / 3.00;
    double o34 = 4.00 / 3.00;
    double fr83 = 8.0 / 3.0;
    //double cg0[6], cg1[6], caa[6], cab[6], cx[6];
    double cg0[7], cg1[7], caa[7], cab[7], cx[7];	//mohan modify 2007-10-13
    double r3q2, r3pi, gr, rho_o3, rho_o34, xa, xa2, ra, rab,
    dra_drho, drab_drho, g, dg, era1, dera1_dra, erab0, derab0_drab,
    ex, dex_drho, uaa, uab, ux, ffaa, ffab,  dffaa_drho, dffab_drho,
    denaa, denab, denx, f83rho, bygr, gaa, gab, gx, taa, tab, txx,
    dgaa_drho, dgab_drho, dgx_drho, dgaa_dgr, dgab_dgr, dgx_dgr;

    r3q2 = std::pow(2.0, (-o3));
    r3pi = std::pow((3.0 / ModuleBase::PI), o3);
    //.....coefficients for pwf correlation......................................
    cg0[1] = 0.0310910;
    cg0[2] = 0.2137000;
    cg0[3] = 7.5957000;
    cg0[4] = 3.5876000;
    cg0[5] = 1.6382000;
    cg0[6] = 0.4929400;

    cg1[1] = 0.0155450;
    cg1[2] = 0.2054800;
    cg1[3] = 14.1189000;
    cg1[4] = 6.1977000;
    cg1[5] = 3.3662000;
    cg1[6] = 0.6251700;
    //......hcth-19-4.....................................
    caa[1] =  0.489508e+00;
    caa[2] = -0.260699e+00;
    caa[3] =  0.432917e+00;
    caa[4] = -0.199247e+01;
    caa[5] =  0.248531e+01;
    caa[6] =  0.200000e+00;

    cab[1] =  0.514730e+00;
    cab[2] =  0.692982e+01;
    cab[3] = -0.247073e+02;
    cab[4] =  0.231098e+02;
    cab[5] = -0.113234e+02;
    cab[6] =  0.006000e+00;

    cx[1] =  0.109163e+01;
    cx[2] = -0.747215e+00;
    cx[3] =  0.507833e+01;
    cx[4] = -0.410746e+01;
    cx[5] =  0.117173e+01;
    cx[6] =   0.004000e+00;
    //...........................................................................
    gr = sqrt(grho);
    rho_o3 = pow(rho, (o3));
    rho_o34 = pow(rho, (o34));
    xa = 1.259921050 * gr / rho_o34;
    xa2 = xa * xa;
    ra = 0.7815926420 / rho_o3;
    rab = r3q2 * ra;
    dra_drho = -0.2605308810 / rho_o34;
    drab_drho = r3q2 * dra_drho;
    XC_Functional::pwcorr(ra, cg1, g, dg);
    era1 = g;
    dera1_dra = dg;
    XC_Functional::pwcorr(rab, cg0, g, dg);
    erab0 = g;
    derab0_drab = dg;
    ex = -0.750 * r3pi * rho_o34;
    dex_drho = -r3pi * rho_o3;
    uaa = caa[6] * xa2;
    uaa = uaa / (1.00 + uaa);
    uab = cab[6] * xa2;
    uab = uab / (1.00 + uab);
    ux = cx[6] * xa2;
    ux = ux / (1.00 + ux);
    ffaa = rho * era1;
    ffab = rho * erab0 - ffaa;
    dffaa_drho = era1 + rho * dera1_dra * dra_drho;
    dffab_drho = erab0 + rho * derab0_drab * drab_drho - dffaa_drho;
    // mb-> i-loop removed
    denaa = 1.0 / (1.00 + caa[6] * xa2);
    denab = 1.0 / (1.00 + cab[6] * xa2);
    denx = 1.0 / (1.00 + cx[6] * xa2);
    f83rho = fr83 / rho;
    bygr = 2.00 / gr;
    gaa = caa[1] + uaa * (caa[2] + uaa * (caa[3] + uaa * (caa[4] + uaa * caa[5])));
    gab = cab[1] + uab * (cab[2] + uab * (cab[3] + uab * (cab[4] + uab * cab[5])));
    gx = cx[1] + ux * (cx[2] + ux * (cx[3] + ux * (cx[4] + ux * cx[5])));
    taa = denaa * uaa * (caa[2] + uaa * (2.0 * caa[3] + uaa
                                         * (3.0 * caa[4] + uaa * 4.0 * caa[5])));
    tab = denab * uab * (cab[2] + uab * (2.0 * cab[3] + uab
                                         * (3.0 * cab[4] + uab * 4.0 * cab[5])));
    txx = denx * ux * (cx[2] + ux * (2.0 * cx[3] + ux
                                     * (3.0 * cx[4] + ux * 4.0 * cx[5])));
    dgaa_drho = -f83rho * taa;
    dgab_drho = -f83rho * tab;
    dgx_drho = -f83rho * txx;
    dgaa_dgr = bygr * taa;
    dgab_dgr = bygr * tab;
    dgx_dgr = bygr * txx;
    // mb
    sx = ex * gx + ffaa * gaa + ffab * gab;
    v1x = dex_drho * gx + ex * dgx_drho
          + dffaa_drho * gaa + ffaa * dgaa_drho
          + dffab_drho * gab + ffab * dgab_drho;
    v2x = (ex * dgx_dgr + ffaa * dgaa_dgr + ffab * dgab_dgr) / gr;
    return;
} //end subroutine hcth

void XC_Functional::pbec(const double &rho, const double &grho, const int &iflag, double &sc, double &v1c, double &v2c)
{
	// PBE correlation (without LDA part)
	// iflag=0: J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996).
	// iflag=1: J.P.Perdew et al., PRL 100, 136406 (2008).
	const double ga = 0.0310906908696548950;
	const double be[2] = {0.06672455060314922, 0.046};
	
	const double third = 1.0 / 3.0;
	const double pi34 = 0.62035049089940;
	const double xkf = 1.9191582926775130;
	const double xks = 1.1283791670955130;
	
	// pi34=(3/4pi)^(1/3), xkf=(9 pi/4)^(1/3), xks= sqrt(4/pi)
	double ec, vc;
		
	const double rs = pi34 / pow(rho, third);
	
	XC_Functional::pw(rs, 0, ec, vc);
	
	const double kf = xkf / rs;
	const double ks = xks * sqrt(kf);
	const double t = sqrt(grho) / (2.0 * ks * rho);
	const double expe = exp(- ec / ga);
	const double af = be[iflag] / ga * (1.0 / (expe - 1.0));
	const double bf = expe * (vc - ec);
	const double y = af * t * t;
	const double xy = (1.0 + y) / (1.0 + y + y * y);
		
	const double x = 1.0 + y + y * y;
	const double qy = y * y * (2.0 + y) / (x * x);
	const double s1 = 1.0 + be[iflag] / ga * t * t * xy;
	const double h0 = ga * log(s1);
	const double dh0 = be[iflag] * t * t / s1 * (- 7.0 / 3.0 * xy - qy * (af * bf /be[iflag] - 7.0 / 3.0));
	const double ddh0 = be[iflag] / (2.0 * ks * ks * rho) * (xy - qy) / s1;

	sc = rho * h0;
	v1c = h0 + dh0;
	v2c = ddh0;
	
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
}

void XC_Functional::optx(const double rho, const double grho, double &sx, double &v1x, double &v2x)
{
    //     OPTX, Handy et al. JCP 116, p. 5411 (2002) and refs. therein
    //     Present release: Mauro Boero, Tsukuba, 10/9/2002
    //--------------------------------------------------------------------------
    //     rhoa = rhob = 0.5 * rho in LDA implementation
    //     grho is the SQUARE of the gradient of rho// --> gr=sqrt(grho)
    //     sx  : total exchange correlation energy at point r
    //     v1x : d(sx)/drho
    //     v2x : 1/gr*d(sx)/d(gr)
    //--------------------------------------------------------------------------
    // use kinds, only: DP
    // implicit none
    // real(kind=DP) :: rho, grho, sx, v1x, v2x

    // parameter :
    double small = 1.e-30;
    double smal2 = 1.e-10;
    //.......coefficients and exponents....................
    // parameter :
    double o43 = 4.00 / 3.00,
                 two13 = 1.2599210498948730,
                         two53 = 3.1748021039363990,
                                 gam = 0.0060,
                                       a1cx = 0.97845711702844210,
                                              a2 = 1.431690;
    double gr, rho43, xa, gamx2, uden, uu;
    //.......OPTX in compact form..........................

    if (rho <= small)
    {
        sx = 0.00;
        v1x = 0.00;
        v2x = 0.00;
    }
    else
    {
        gr = (grho > smal2) ? grho : smal2;	//max()
        rho43 = pow(rho, o43);
        xa = two13 * sqrt(gr) / rho43;
        gamx2 = gam * xa * xa;
        uden = 1.e+00 / (1.e+00 + gamx2);
        uu = a2 * gamx2 * gamx2 * uden * uden;
        uden = rho43 * uu * uden;
        sx = -rho43 * (a1cx + uu) / two13;
        v1x = o43 * (sx + two53 * uden) / rho;
        v2x = -two53 * uden / gr;
    } //endif

    return;
} // end subroutine optx

void XC_Functional::perdew86(const double rho, const double grho, double &sc, double &v1c, double &v2c)
{
    // Perdew gradient correction on correlation: PRB 33, 8822 (1986)

    // USE kinds
    // implicit none
    // real(kind=DP) :: rho, grho, sc, v1c, v2c
    double p1, p2, p3, p4, pc1, pc2, pci;
    // parameter :
    p1 = 0.0232660;
    p2 = 7.389e-6;
    p3 = 8.7230;
    p4 = 0.4720;
    // parameter :
    pc1 = 0.0016670;
    pc2 = 0.0025680;
    pci = pc1 + pc2;
    double third, pi34;
    // parameter :
    third = 1.0 / 3.0;
    pi34 = 0.62035049089940;  // pi34=(3/4pi)^(1/3)
    double rho13, rho43, rs, rs2, rs3, cna, cnb, cn, drs;
    double dcna, dcnb, dcn, phi, ephi;

    rho13 = pow(rho, third);
    rho43 = pow(rho13, 4);
    rs = pi34 / rho13;
    rs2 = rs * rs;
    rs3 = rs * rs2;
    cna = pc2 + p1 * rs + p2 * rs2;
    cnb = 1.0 + p3 * rs + p4 * rs2 + 1.e4 * p2 * rs3;
    cn = pc1 + cna / cnb;
    drs = - third * pi34 / rho43;
    dcna = (p1 + 2.0 * p2 * rs) * drs;
    dcnb = (p3 + 2.0 * p4 * rs + 3.e4 * p2 * rs2) * drs;
    dcn = dcna / cnb - cna / (cnb * cnb) * dcnb;
    phi = 0.1920 * pci / cn * sqrt(grho) * pow(rho, (- 7.0 / 6.0));
    // SdG: in the original paper 1.745*0.11=0.19195 is used
    ephi = exp(- phi);
    sc = grho / rho43 * cn * ephi;
    v1c = sc * ((1.0 + phi) * dcn / cn - ((4.0 / 3.0) - (7.0 /
                                          6.0) * phi) / rho);
    v2c = cn * ephi / rho43 * (2.0 - phi);

    return;
} //end subroutine perdew86

void XC_Functional::ggax(const double &rho, const double &grho, double &sx, double &v1x, double &v2x)
{
    //-----------------------------------------------------------------------
    // Perdew-Wang GGA (PW91), exchange part:
    // J.P. Perdew et al.,PRB 46, 6671 (1992)
    const double f1 = 0.196450;
    const double f2 = 7.79560;
    const double f3 = 0.27430;
    const double f4 = 0.15080;
    const double f5 = 0.0040;
    const double fp1 = -0.0192920212964260;
    const double fp2 = 0.1616204596739950;
    // fp1 = -3/(16 pi)*(3 pi^2)^(-1/3)
    double rhom43, s, s2, s3, s4, exps, as, sa2b8, shm1, bs, das,
    dbs, dls;

    rhom43 = pow(rho, (- 4.0 / 3.0));
    s = fp2 * sqrt(grho) * rhom43;
    s2 = s * s;
    s3 = s2 * s;
    s4 = s2 * s2;
    exps = f4 * exp(- 100.0 * s2);
    as = f3 - exps - f5 * s2;
    sa2b8 = sqrt(1.00 + f2 * f2 * s2);
    shm1 = log(f2 * s + sa2b8);
    bs = 1.0 + f1 * s * shm1 + f5 * s4;
    das = (200.0 * exps - 2.0 * f5) * s;
    dbs = f1 * (shm1 + f2 * s / sa2b8) + 4.0 * f5 * s3;
    dls = (das / as - dbs / bs);
    
	sx = fp1 * grho * rhom43 * as / bs;
	v1x = - 4.0 / 3.0 * sx / rho * (1.0 + s * dls);
    v2x = fp1 * rhom43 * as / bs * (2.0 + s * dls);

    return;
} //end subroutine ggax



void XC_Functional::ggac(const double &rho,const double &grho, double &sc, double &v1c, double &v2c)
{
    // Perdew-Wang GGA (PW91) correlation part
    double al, pa, pb, pc, pd, cx, cxc0, cc0;
    // parameter :
    al = 0.090;
    pa = 0.0232660;
    pb = 7.389e-6;
    pc = 8.7230;
    pd = 0.4720;
    // parameter :
    cx = -0.0016670;
    cxc0 = 0.0025680;
    cc0 = - cx + cxc0;
    double third, pi34, nu, be, xkf, xks;
    // parameter :
    third = 1.0 / 3.0;
    pi34 = 0.62035049089940;
    // parameter :
    nu = 15.7559203494831440;
    be = nu * cc0;
    // parameter :
    xkf = 1.9191582926775130;
    xks = 1.1283791670955130;
    // pi34=(3/4pi)^(1/3),  nu=(16/pi)*(3 pi^2)^(1/3)
    // xkf=(9 pi/4)^(1/3), xks= sqrt(4/pi)
    double kf, ks, rs, rs2, rs3, ec, vc, t, expe, af, bf, y, xy,
    qy, s1;
    double h0, dh0, ddh0, ee, cn, dcn, cna, dcna, cnb, dcnb, h1,
    dh1, ddh1;
	
	rs = pi34 / pow(rho, third);
    rs2 = rs * rs;
    rs3 = rs * rs2;
	
	// call function: pw
	XC_Functional::pw(rs, 0, ec, vc);
    kf = xkf / rs;
    ks = xks * sqrt(kf);
    t = sqrt(grho) / (2.0 * ks * rho);
    expe = exp(- 2.0 * al * ec / (be * be));
    af = 2.0 * al / be * (1.0 / (expe - 1.0));
    bf = expe * (vc - ec);
    y = af * t * t;
    xy = (1.0 + y) / (1.0 + y + y * y);

    double x;
    x = 1.0 + y + y * y;
    qy = y * y * (2.0 + y) / (x * x);
    s1 = 1.0 + 2.0 * al / be * t * t * xy;
    h0 = be * be / (2.0 * al) * log(s1);
    dh0 = be * t * t / s1 * (- 7.0 / 3.0 * xy - qy * (af * bf /
                             be - 7.0 / 3.0));
    ddh0 = be / (2.0 * ks * ks * rho) * (xy - qy) / s1;

    x = ks / kf * t;
    ee = - 100.0 * (x * x);
    cna = cxc0 + pa * rs + pb * rs2;
    dcna = pa * rs + 2.0 * pb * rs2;
    cnb = 1.0 + pc * rs + pd * rs2 + 1.e4 * pb * rs3;
    dcnb = pc * rs + 2.0 * pd * rs2 + 3.e4 * pb * rs3;
    cn = cna / cnb - cx;
    dcn = dcna / cnb - cna * dcnb / (cnb * cnb);
    h1 = nu * (cn - cc0 - 3.0 / 7.0 * cx) * t * t * exp(ee);
    dh1 = - third * (h1 * (7.0 + 8.0 * ee) + nu * t * t * exp(ee)
                     * dcn);
    ddh1 = 2.0 * h1 * (1.0 + ee) * rho / grho;
    sc = rho * (h0 + h1);
    v1c = h0 + h1 + dh0 + dh1;
    v2c = ddh0 + ddh1;

    return;
} //end subroutine ggac

void XC_Functional::wcx(const double &rho,const double &grho, double &sx, double &v1x, double &v2x)
{
  double kf, agrho, s1, s2, es2, ds, dsg, exunif, fx;
  // (3*pi2*|rho|)^(1/3)
  // |grho|
  // |grho|/(2*kf*|rho|)
  // s^2
  // n*ds/dn
  // n*ds/d(gn)
  // exchange energy LDA part
  // exchange energy gradient part
  double dxunif, dfx, f1, f2, f3, dfx1, x1, x2, x3, dxds1, dxds2, dxds3;
  // numerical coefficients (NB: c2=(3 pi^2)^(1/3) )
  double third, c1, c2, c5, teneightyone;	// c6
  
  third = 1.0/3.0;
  c1 = 0.75/ModuleBase::PI;
  c2 = 3.093667726280136;
  c5 = 4.0 * third;
  teneightyone = 0.123456790123;
  // parameters of the functional
  double k, mu, cwc;
  k = 0.804;
  mu = 0.2195149727645171;
  cwc = 0.00793746933516;
  //
  agrho = sqrt (grho);
  kf = c2 * pow(rho,third);
  dsg = 0.5 / kf;
  s1 = agrho * dsg / rho;
  s2 = s1 * s1;
  es2 = exp(-s2);
  ds = - c5 * s1;
  //
  //   Energy
  //
  // x = 10/81 s^2 + (mu - 10/81) s^2 e^-s^2 + ln (1 + c s^4)
  x1 = teneightyone * s2;
  x2 = (mu - teneightyone) * s2 * es2;
  x3 = log(1.0 + cwc * s2 * s2);
  f1 = (x1 + x2 + x3) / k;
  f2 = 1.0 + f1;
  f3 = k / f2;
  fx = k - f3;
  exunif = - c1 * kf;
  sx = exunif * fx;
  //
  //   Potential
  //
  dxunif = exunif * third;
  dfx1 = f2 * f2;
  dxds1 = teneightyone;
  dxds2 = (mu - teneightyone) * es2 * (1.0 - s2);
  dxds3 = 2.0 * cwc * s2 / (1.0 + cwc * s2 *s2);
  dfx = 2.0 * s1 * (dxds1 + dxds2 + dxds3) / dfx1;
  v1x = sx + dxunif * fx + exunif * dfx * ds;
  v2x = exunif * dfx * dsg / agrho;

  sx = sx * rho;
  return;
}

void XC_Functional::becke88(const double &rho, const double &grho, double &sx, double &v1x, double &v2x)
{
    //-----------------------------------------------------------------------
    // Becke exchange: A.D. Becke, PRA 38, 3098 (1988)
    // only gradient-corrected part, no Slater term included
    //
    double beta, third, two13;
    beta = 0.00420;
    third = 1.0 / 3.0;
    two13 = 1.2599210498948730;
    double rho13, rho43, xs, xs2, sa2b8, shm1, dd, dd2, ee;

    rho13 = pow(rho, third);
    rho43 = pow(rho13, 4);
    xs = two13 * sqrt(grho) / rho43;
    xs2 = xs * xs;
    sa2b8 = sqrt(1.00 + xs2);
    shm1 = log(xs + sa2b8);
    dd = 1.00 + 6.00 * beta * xs * shm1;
    dd2 = dd * dd;
    ee = 6.00 * beta * xs2 / sa2b8 - 1.0;
    sx = two13 * grho / rho43 * (- beta / dd);
    v1x = - (4.0 / 3.0) / two13 * xs2 * beta * rho13 * ee / dd2;
    v2x = two13 * beta * (ee - dd) / (rho43 * dd2);

    return;
} // end subroutine becke88


void XC_Functional::glyp(const double &rho, const double &grho, double &sc, double &v1c, double &v2c)
{
    //-----------------------------------------------------------------------
    // Lee Yang Parr: gradient correction part

    double a, b, c, d;
    a = 0.049180;
    b = 0.1320;
    c = 0.25330;
    d = 0.3490;
    double rhom13, rhom43, rhom53, om, xl, ff, dom, dxl;

    rhom13 = pow(rho, (- 1.0 / 3.0));
    om = exp(- c * rhom13) / (1.0 + d * rhom13);
    xl = 1.0 + (7.0 / 3.0) * (c * rhom13 + d * rhom13 / (1.0 + d *
                              rhom13));
    ff = a * b * grho / 24.0;
    rhom53 = pow(rhom13, 5);
    sc = ff * rhom53 * om * xl;
    dom = - om * (c + d + c * d * rhom13) / (1.0 + d * rhom13);

    double x;
    x = 1.0 + d * rhom13;
    dxl = (7.0 / 3.0) * (c + d + 2.0 * c * d * rhom13 + c * d * d *
                         rhom13 * rhom13) / (x * x);
    rhom43 = pow(rhom13, 4);
    v1c = - ff * rhom43 / 3.0 * (5.0 * rhom43 * om * xl + rhom53 *
                                 dom * xl + rhom53 * om * dxl);
    v2c = 2.0 * sc / grho;

    return;
} // end subroutine glyp

void XC_Functional::gcxc_spin_libxc(double rhoup, double rhodw, 
        ModuleBase::Vector3<double> gdr1, ModuleBase::Vector3<double> gdr2,
        double &sxc, double &v1xcup, double &v1xcdw, double &v2xcup, double &v2xcdw, double &v2xcud)
{
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

void XC_Functional::becke88_spin(double rho, double grho, double &sx, double &v1x, double &v2x)
{
    //-------------------------------------------------------------------
    // Becke exchange: A.D. Becke, PRA 38, 3098 (1988) - Spin polarized case

    // USE kinds
    // implicit none
    // real(kind=DP) :: rho, grho, sx, v1x, v2x
    // input: charge
    // input: gradient
    // output: the up and down energies
    // output: first part of the potential
    // output: the second part of the potential

    double beta, third;
    // parameter :
    beta = 0.00420;
    third = 1.0 / 3.0;
    double rho13, rho43, xs, xs2, sa2b8, shm1, dd, dd2, ee;

    rho13 = pow(rho, third);
    rho43 = pow(rho13, 4);
    xs = sqrt(grho) / rho43;
    xs2 = xs * xs;
    sa2b8 = sqrt(1.00 + xs2);
    shm1 = log(xs + sa2b8);
    dd = 1.00 + 6.00 * beta * xs * shm1;
    dd2 = dd * dd;
    ee = 6.00 * beta * xs2 / sa2b8 - 1.0;
    sx = grho / rho43 * (- beta / dd);
    v1x = - (4.0 / 3.0) * xs2 * beta * rho13 * ee / dd2;
    v2x = beta * (ee - dd) / (rho43 * dd2);

    return;
} //end subroutine becke88_spin

//
//-----------------------------------------------------------------------
void XC_Functional::perdew86_spin(double rho, double zeta, double grho, double &sc,
                   double &v1cup, double &v1cdw, double &v2c)
{
    //-------------------------------------------------------------------
    // Perdew gradient correction on correlation: PRB 33, 8822 (1986)
    // spin-polarized case

    // USE kinds
    // implicit none
    // real(kind=DP) :: rho, zeta, grho, sc, v1cup, v1cdw, v2c
    double p1, p2, p3, p4, pc1, pc2, pci;
    // parameter :
    p1 = 0.0232660;
    p2 = 7.389e-6;
    p3 = 8.7230;
    p4 = 0.4720;
    // parameter :
    pc1 = 0.0016670;
    pc2 = 0.0025680;
    pci = pc1 + pc2;
    double third, pi34;
    // parameter :
    third = 1.0 / 3.0;
    pi34 = 0.62035049089940;
    // pi34=(3/4pi)^(1/3)

    double rho13, rho43, rs, rs2, rs3, cna, cnb, cn, drs;
    double dcna, dcnb, dcn, phi, ephi, dd, ddd;

    rho13 = pow(rho, third);
    rho43 = pow(rho13, 4);
    rs = pi34 / rho13;
    rs2 = rs * rs;
    rs3 = rs * rs2;
    cna = pc2 + p1 * rs + p2 * rs2;
    cnb = 1.0 + p3 * rs + p4 * rs2 + 1.e4 * p2 * rs3;
    cn = pc1 + cna / cnb;
    drs = - third * pi34 / rho43;
    dcna = (p1 + 2.0 * p2 * rs) * drs;
    dcnb = (p3 + 2.0 * p4 * rs + 3.e4 * p2 * rs2) * drs;
    dcn = dcna / cnb - cna / (cnb * cnb) * dcnb;
    phi = 0.1920 * pci / cn * sqrt(grho) * pow(rho, (- 7.0 / 6.0));
    //SdG: in the original paper 1.745*0.11=0.19195 is used
    dd = pow((2.0) , third) * sqrt(pow(((1.0 + zeta) * 0.50) , (5.0 /
                                       3.0)) + pow(((1.0 - zeta) * 0.50) , (5.0 / 3.0)));
    ddd = pow((2.0) , (- 4.0 / 3.0)) * 5.0 * (pow(((1.0 + zeta)
            * 0.50) , (2.0 / 3.0)) - pow(((1.0 - zeta) * 0.50) , (2.0 /
                                         3.0))) / (3.0 * dd);
    ephi = exp(- phi);
    sc = grho / rho43 * cn * ephi / dd;
    v1cup = sc * ((1.0 + phi) * dcn / cn - ((4.0 / 3.0) -
                                            (7.0 / 6.0) * phi) / rho) - sc * ddd / dd * (1.0 - zeta)
            / rho;
    v1cdw = sc * ((1.0 + phi) * dcn / cn - ((4.0 / 3.0) -
                                            (7.0 / 6.0) * phi) / rho) + sc * ddd / dd * (1.0 + zeta)
            / rho;
    v2c = cn * ephi / rho43 * (2.0 - phi) / dd;

    return;
} //end subroutine perdew86_spin

//
//-----------------------------------------------------------------------
void XC_Functional::ggac_spin(double rho, double zeta, double grho, double &sc,
               double &v1cup, double &v1cdw, double &v2c)
{
    //-------------------------------------------------------------------
    // Perdew-Wang GGA (PW91) correlation part - spin-polarized

    // USE kinds
    // implicit none
    // real(kind=DP) :: rho, zeta, grho, sc, v1cup, v1cdw, v2c
    double al, pa, pb, pc, pd, cx, cxc0, cc0;
    // parameter :
    al = 0.090;
    pa = 0.0232660;
    pb = 7.389e-6;
    pc = 8.7230;
    pd = 0.4720;
    // parameter :
    cx = - 0.0016670;
    cxc0 = 0.0025680;
    cc0 = - cx + cxc0;

    double third, pi34, nu, be, xkf, xks;
    // parameter :
    third = 1.0 / 3.0;
    pi34 = 0.6203504908994e0l;
    // parameter :
    nu = 15.7559203494831440;
    be = nu * cc0;
    // parameter :
    xkf = 1.9191582926775130;
    xks = 1.1283791670955130;
    // pi34=(3/4pi)^(1/3),  nu=(16/pi)*(3 pi^2)^(1/3)
    // xkf=(9 pi/4)^(1/3), xks= sqrt(4/pi)
    double kf, ks, rs, rs2, rs3, ec, vcup, vcdw, t, expe, af, y,
    xy, qy, s1, h0, ddh0, ee, cn, dcn, cna, dcna, cnb, dcnb, h1, dh1,
    ddh1, fz, fz2, fz3, fz4, dfz, bfup, dh0up, dh0zup, //bfdw, dh0dw,
    dh0zdw, dh1zup, dh1zdw;

    rs = pi34 / pow(rho, third);
    rs2 = rs * rs;
    rs3 = rs * rs2;
	XC_Functional::pw_spin(rs, zeta, ec, vcup, vcdw);
    kf = xkf / rs;
    ks = xks * sqrt(kf);
    fz = 0.50 * (pow((1.0 + zeta) , (2.0 / 3.0)) + pow((1.0 - zeta) , (
                     2.0 / 3.0)));
    fz2 = fz * fz;
    fz3 = fz2 * fz;
    fz4 = fz3 * fz;
    dfz = (pow((1.0 + zeta) , (- 1.0 / 3.0)) - pow((1.0 - zeta) , (-
            1.0 / 3.0))) / 3.0;
    t = sqrt(grho) / (2.0 * fz * ks * rho);
    expe = exp(- 2.0 * al * ec / (fz3 * be * be));
    af = 2.0 * al / be * (1.0 / (expe - 1.0));
    bfup = expe * (vcup - ec) / fz3;
    //bfdw = expe * (vcdw - ec) / fz3;
    y = af * t * t;
    xy = (1.0 + y) / (1.0 + y + y * y);
    qy = y * y * (2.0 + y) / (1.0 + y + y * y) ;	//**2;
    qy *= qy;
    s1 = 1.0 + 2.0 * al / be * t * t * xy;
    h0 = fz3 * be * be / (2.0 * al) * log(s1);
    dh0up = be * t * t * fz3 / s1 * (- 7.0 / 3.0 * xy - qy *
                                     (af * bfup / be - 7.0 / 3.0));
    //dh0dw = be * t * t * fz3 / s1 * (- 7.0 / 3.0 * xy - qy *
    //                                 (af * bfdw / be - 7.0 / 3.0));
    dh0zup = (3.0 * h0 / fz - be * t * t * fz2 / s1 * (2.0 * xy -
              qy * (3.0 * af * expe * ec / fz3 / be + 2.0))) * dfz * (1.0 -
                      zeta);
    dh0zdw = - (3.0 * h0 / fz - be * t * t * fz3 / s1 * (2.0 * xy -
                qy * (3.0 * af * expe * ec / fz3 / be + 2.0))) * dfz * (1.0 +
                        zeta);
    ddh0 = be * fz / (2.0 * ks * ks * rho) * (xy - qy) / s1;
    ee = - 100.0 * fz4 * (ks / kf * t);	// **2;
    ee *= ee;
    cna = cxc0 + pa * rs + pb * rs2;
    dcna = pa * rs + 2.0 * pb * rs2;
    cnb = 1.0 + pc * rs + pd * rs2 + 1.e4 * pb * rs3;
    dcnb = pc * rs + 2.0 * pd * rs2 + 3.e4 * pb * rs3;
    cn = cna / cnb - cx;
    dcn = dcna / cnb - cna * dcnb / (cnb * cnb);
    h1 = nu * (cn - cc0 - 3.0 / 7.0 * cx) * fz3 * t * t * exp(ee);
    dh1 = - third * (h1 * (7.0 + 8.0 * ee) + fz3 * nu * t * t * exp
                     (ee) * dcn);
    ddh1 = 2.0 * h1 * (1.0 + ee) * rho / grho;
    dh1zup = (1.0 - zeta) * dfz * h1 * (1.0 + 2.0 * ee / fz);
    dh1zdw = - (1.0 + zeta) * dfz * h1 * (1.0 + 2.0 * ee / fz);
    sc = rho * (h0 + h1);
    v1cup = h0 + h1 + dh0up + dh1 + dh0zup + dh1zup;
    v1cdw = h0 + h1 + dh0up + dh1 + dh0zdw + dh1zdw;
    v2c = ddh0 + ddh1;
    return;
} // end subroutine ggac_spin

//
//---------------------------------------------------------------
void XC_Functional::pbec_spin(double rho, double zeta, double grho, const int &iflag, double &sc,
               double &v1cup, double &v1cdw, double &v2c)
{
    //-----------------------------------------------------------

    // PBE correlation (without LDA part) - spin-polarized
    // J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996).

    // USE kinds
    // implicit none
    // real(kind=DP) :: rho, zeta, grho, sc, v1cup, v1cdw, v2c
    double ga, be[3];//mohan add
    // parameter :
    ga = 0.0310910;
    be[1] = 0.06672455060314922;//zhengdy add 2019-09-12, ensure same parameter with another dft code.
	be[2] = 0.0460000;//mohan add 2012-05-28
    double third, pi34, xkf, xks;
    // parameter :
    third = 1.0 / 3.0;
    pi34 = 0.62035049089940;
    // parameter :
    xkf = 1.9191582926775130;
    xks = 1.1283791670955130;
    // pi34=(3/4pi)^(1/3), xkf=(9 pi/4)^(1/3), xks= sqrt(4/pi)
    double kf, ks, rs, ec, vcup, vcdw, t, expe, af, y, xy, qy,
    s1, h0, ddh0;
    double fz, fz2, fz3, dfz, bfup, bfdw, dh0up, dh0dw, dh0zup, dh0zdw; // fz4

    rs = pi34 / pow(rho, third);
	XC_Functional::pw_spin(rs, zeta, ec, vcup, vcdw); //mohan fix bug 2012-05-28
    kf = xkf / rs;
    ks = xks * sqrt(kf);
    fz = 0.50 * (pow((1.0 + zeta) , (2.0 / 3.0)) + pow((1.0 - zeta) , (
                     2.0 / 3.0)));
    fz2 = fz * fz;
    fz3 = fz2 * fz;
    //fz4 = fz3 * fz;
    dfz = (pow((1.0 + zeta) , (- 1.0 / 3.0)) - pow((1.0 - zeta) , (-
            1.0 / 3.0))) / 3.0;
    t = sqrt(grho) / (2.0 * fz * ks * rho);
    expe = exp(- ec / (fz3 * ga));
    af = be[iflag] / ga * (1.0 / (expe - 1.0));
    bfup = expe * (vcup - ec) / fz3;
    bfdw = expe * (vcdw - ec) / fz3;
    y = af * t * t;
    xy = (1.0 + y) / (1.0 + y + y * y);
    qy = y * y * (2.0 + y) / pow( (1.0 + y + y * y),2 );
//    qy *= qy; //bug
    s1 = 1.0 + be[iflag] / ga * t * t * xy;
    h0 = fz3 * ga * log(s1);
    dh0up = be[iflag] * t * t * fz3 / s1 * (- 7.0 / 3.0 * xy - qy *
                                     (af * bfup / be[iflag] - 7.0 / 3.0));
    dh0dw = be[iflag] * t * t * fz3 / s1 * (- 7.0 / 3.0 * xy - qy *
                                     (af * bfdw / be[iflag] - 7.0 / 3.0));
    dh0zup = (3.0 * h0 / fz - be[iflag] * t * t * fz2 / s1 * (2.0 * xy -
              qy * (3.0 * af * expe * ec / fz3 / be[iflag] + 2.0))) * dfz * (1.0 -
                      zeta);
    dh0zdw = - (3.0 * h0 / fz - be[iflag] * t * t * fz2 / s1 * (2.0 * xy -
                qy * (3.0 * af * expe * ec / fz3 / be[iflag] + 2.0))) * dfz * (1.0 +
                        zeta);
    ddh0 = be[iflag] * fz / (2.0 * ks * ks * rho) * (xy - qy) / s1;
    sc = rho * h0;
    v1cup = h0 + dh0up + dh0zup;
    v1cdw = h0 + dh0dw + dh0zdw;
    v2c = ddh0;
    return;
} // end subroutine pbec_spin