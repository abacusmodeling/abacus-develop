#include "../src_pw/tools.h"
#include "../src_pw/global.h"
#include "../src_pw/functional.h"
#include "myfunc.h"

//-----------------------------------------------------------------------
void becke88(double rho, double grho, double &sx, double &v1x, double &v2x)
{
    //-----------------------------------------------------------------------
    // Becke exchange: A.D. Becke, PRA 38, 3098 (1988)
    // only gradient-corrected part, no Slater term included
    //
    // USE kinds
    // implicit none
    // real(kind=DP) :: rho, grho, sx, v1x, v2x
    double beta, third, two13;
    // parameter :
    beta = 0.00420;
    // parameter :
    third = 1.0 / 3.0;
    two13 = 1.2599210498948730;
    // two13 = 2^(1/3)
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

//
//-----------------------------------------------------------------------
void ggax(double rho, double grho, double &sx, double &v1x, double &v2x)
{
    //-----------------------------------------------------------------------
    // Perdew-Wang GGA (PW91), exchange part:
    // J.P. Perdew et al.,PRB 46, 6671 (1992)

    // USE kinds
    // implicit none
    // real(kind=DP) :: rho, grho, sx, v1x, v2x
    double f1, f2, f3, f4, f5;
    // parameter :
    f1 = 0.196450;
    f2 = 7.79560;
    f3 = 0.27430;
    f4 = 0.15080;
    f5 = 0.0040;
    double fp1, fp2;
    // parameter :
    fp1 = -0.0192920212964260;
    fp2 = 0.1616204596739950;
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

//
//-----------------------------------------------------------------------
void perdew86(double rho, double grho, double &sc, double &v1c, double &v2c)
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

//
//-----------------------------------------------------------------------
void glyp(double rho, double grho, double &sc, double &v1c, double &v2c)
{
    //-----------------------------------------------------------------------
    // Lee Yang Parr: gradient correction part

    // USE kinds
    // implicit none
    // real(kind=DP) :: rho, grho, sc, v1c, v2c
    double a, b, c, d;
    // parameter :
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


//     ==================================================================
void hcth(double rho, double grho, double &sx, double &v1x, double &v2x)
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
    r3pi = std::pow((3.0 / PI), o3);
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
    pwcorr(ra, cg1, g, dg);
    era1 = g;
    dera1_dra = dg;
    pwcorr(rab, cg0, g, dg);
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

//-------------------------------------------------------------------=
void pwcorr(double r, double c[], double &g, double &dg)
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

//-----------------------------------------------------------------------------
//     ==================================================================
void optx(double rho, double grho, double &sx, double &v1x, double &v2x)
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
