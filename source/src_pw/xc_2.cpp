#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <stdexcept>

using namespace std;

#include "../src_pw/global.h"
#include "../src_pw/xc_type.h"
#include "../src_pw/xc_functional.h"
#include "myfunc.h"

xcfunc xc2;

//-----------------------------------------------------------------------
void gcx_spin(double rhoup, double rhodw, double grhoup2, double grhodw2,
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

	//cout << " xcf.igcx_now=" << xcf.igcx_now << endl;

    if (rho <= small || xcf.igcx_now == 0)
    {
        sx = 0.00;
        v1xup = 0.00;
        v2xup = 0.00;
        v1xdw = 0.00;
        v2xdw = 0.00;
    }
    else if (xcf.igcx_now == 1)
    {
        if (rhoup > small && sqrt(fabs(grhoup2)) > small)
        {
            becke88_spin(rhoup, grhoup2, sxup, v1xup, v2xup);
        }
        else
        {
            sxup = 0.0;
            v1xup = 0.0;
            v2xup = 0.0;
        } //  endif

        if (rhodw > small && sqrt(fabs(grhodw2)) > small)
        {
            becke88_spin(rhodw, grhodw2, sxdw, v1xdw, v2xdw);
        }
        else
        {
            sxdw = 0.0;
            v1xdw = 0.0;
            v2xdw = 0.0;
        } //endif

        sx = sxup + sxdw;
    }
    else if (xcf.igcx_now == 2)
    {
        if (rhoup > small && sqrt(abs(grhoup2)) > small)
        {
            //XC_Functional::ggax(2.0 * rhoup, 4.0 * grhoup2, sxup, v1xup, v2xup);
			throw runtime_error("sxup, v1xup, v2xup uninitialized in "+TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));
        }
        else
        {
            sxup = 0.0;
            v1xup = 0.0;
            v2xup = 0.0;
        } //endif

        if (rhodw > small && sqrt(abs(grhodw2)) > small)
        {
            //XC_Functional::ggax(2.0 * rhodw, 4.0 * grhodw2, sxdw, v1xdw, v2xdw);
			throw runtime_error("sxdw, v1xdw, v2xdw uninitialized in "+TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));
        }
        else
        {
            sxdw = 0.0;
            v1xdw = 0.0;
            v2xdw = 0.0;
        } //endif

        sx = 0.50 * (sxup + sxdw);

        v2xup = 2.0 * v2xup;

        v2xdw = 2.0 * v2xdw;
    }
	// igcx=3: PBE
	// igcx=4: revised PBE
	// igcx=8: PBE0
	// igcx=10: PBEsol
    else if (xcf.igcx_now == 3 || xcf.igcx_now == 4 || xcf.igcx_now==8 || xcf.igcx_now == 10)
    {
        if (xcf.igcx_now == 4)//mohan fix bug 2012-05-28
        {
            iflag = 1; //minus 1 in C++
        }
		else if(xcf.igcx_now == 10)
		{
			iflag = 2;
		}
        else
        {
            iflag = 0;
        } 

        if (rhoup > small && sqrt(fabs(grhoup2)) > small)
        {
            XC_Functional::pbex(2.0 * rhoup, 4.0 * grhoup2, iflag, sxup, v1xup, v2xup);
        }
        else
        {
            sxup = 0.0;
            v1xup = 0.0;
            v2xup = 0.0;
        } 
		
		if (xcf.igcx_now == 8)						// Peize Lin add 2017-10-23
		{
			sxup *= 0.75;
			v1xup *= 0.75;
			v2xup *= 0.75;
		}

        if (rhodw > small && sqrt(fabs(grhodw2)) > small)
        {
            XC_Functional::pbex(2.0 * rhodw, 4.0 * grhodw2, iflag, sxdw, v1xdw, v2xdw);
        }
        else
        {
            sxdw = 0.0;
            v1xdw = 0.0;
            v2xdw = 0.0;
        } 
		
		if (xcf.igcx_now == 8)						// Peize Lin add 2017-10-23
		{
			sxdw *= 0.75;
			v1xdw *= 0.75;
			v2xdw *= 0.75;
		}		

        sx = 0.5 * (sxup + sxdw);
        v2xup = 2.0 * v2xup;
        v2xdw = 2.0 * v2xdw;
    }
	// igcx=12: HSE
	else if (xcf.igcx_now == 12)
	{
		throw domain_error( "HSE unfinished in "+TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));
	}
    else
    {
        cout << "\n gcx_spin, not implemented, igcx_now";
    } //endif

    return;
} //end subroutine gcx_spin

//
//-----------------------------------------------------------------------
void gcc_spin(double rho, double &zeta, double grho, double &sc,
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

    if (abs(zeta) - 1.0 > small)
    {
        sc = 0.00;
        v1cup = 0.00;
        v1cdw = 0.00;
        v2c = 0.00;
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


//	cout << "xcf.igcc_now=" << xcf.igcc_now << endl;

    if (xcf.igcc_now == 0 || rho <= small || sqrt(abs(grho)) <= small)
    {
        sc = 0.00;
        v1cup = 0.00;
        v1cdw = 0.00;
        v2c = 0.00;
    }
    else if (xcf.igcc_now == 1)
    {
        perdew86_spin(rho, zeta, grho, sc, v1cup, v1cdw, v2c);
    }
    else if (xcf.igcc_now == 2)
    {
        ggac_spin(rho, zeta, grho, sc, v1cup, v1cdw, v2c);
    }
    else if (xcf.igcc_now == 3 || xcf.igcc_now > 4)
    {
        cout << "\n lsda_functionals, not implemented, igcc_now = "
             << xcf.igcc_now;
    }
    else if (xcf.igcc_now == 4)
    {
        pbec_spin(rho, zeta, grho, 1, sc, v1cup, v1cdw, v2c);
    }
	else if (xcf.igcc_now == 8)//mohan add 2012-05-28
	{
        pbec_spin(rho, zeta, grho, 2, sc, v1cup, v1cdw, v2c);
	}
    else
    {
        sc = 0.00;
        v1cup = 0.00;
        v1cdw = 0.00;
        v2c = 0.00;
    } //endif

    return;
} //end subroutine gcc_spin

//
//-----------------------------------------------------------------------
void pz_polarized(double rs, double &ec, double &vc)
{
    //-------------------------------------------------------------------
    //     J.P. Perdew and A. Zunger, PRB 23, 5048 (1981)
    //     spin-polarized energy and potential

    // USE kinds
    // implicit none
    // real(kind=DP) :: rs, ec, vc
    double a, b, c, d, gc, b1, b2;
    // parameter
    a = 0.015550;
    b = - 0.02690;
    c = 0.00070;
    d = - 0.00480;
    gc = - 0.08430;
    b1 = 1.39810;
    b2 = 0.26110;
    double lnrs, rs12, ox, dox;
//	bool xc_rel;
    // PARAMETER :
    //double xcprefact = 0.022575584;
    //double pi34 = 0.62035049089940;
//	double betha, etha, csi, prefact;

    if (rs < 1.00)
    {
        // high density formula
        lnrs = log(rs);
        ec = a * lnrs + b + c * rs * lnrs + d * rs;
        vc = a * lnrs + (b - a / 3.0) + 2.0 / 3.0 * c * rs * lnrs +
             (2.0 * d - c) / 3.0 * rs;
    }
    else
    {
        // interpolation formula
        rs12 = sqrt(rs);
        ox = 1.0 + b1 * rs12 + b2 * rs;
        dox = 1.0 + 7.0 / 6.0 * b1 * rs12 + 4.0 / 3.0 * b2 * rs;
        ec = gc / ox;
        vc = ec * dox / ox;
    } //endif

    /*
    if ( lxc_rel ){
    	betha = prefact * pi34 / rs;
    	etha = DSQRT( 1 + betha**2 );
    	csi = betha + etha;
    	prefact = 1.00 - (3.00/2.00) * ( (betha*etha - log(csi))/betha**2 )**2;
    	ec = ec * prefact;
    	vc = vc * prefact;
    } //ENDIF
    */
    return;
} //end subroutine pz_polarized

void becke88_spin(double rho, double grho, double &sx, double &v1x, double &v2x)
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
void perdew86_spin(double rho, double zeta, double grho, double &sc,
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
void ggac_spin(double rho, double zeta, double grho, double &sc,
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
void pbec_spin(double rho, double zeta, double grho, const int &iflag, double &sc,
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

double min3(double x1, double x2, double x3)
{
    double x;

    if (x1 < x2)
    {
        if (x1 < x3)
            x = x1;
        else
            x = x3;
    }
    else
    {
        if (x2 < x3)
            x = x2;
        else
            x = x3;
    }

    return x;
}
