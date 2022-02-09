//==========================================================
// AUTHOR : mohan
// DATE : 2009-04-04
//==========================================================
#ifndef XC_FUNCTIONAL_H
#define XC_FUNCTIONAL_H

#ifdef USE_LIBXC
#include <xc.h>
#endif	// ifdef USE_LIBXC

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/vector3.h"
class XC_Functional
{
	public:

	XC_Functional();
	~XC_Functional();

	// LDA
	static void xc(const double &rho, double &ex, double &ec, double &vx, double &vc);

	// LSDA
	static void xc_spin(const double &rho, const double &zeta,
			double &ex, double &ec,
			double &vxup, double &vxdw,
			double &vcup, double &vcdw);

	// GGA
	static void gcxc(const double &rho, const double &grho, double &sx, double &sc,
        	double &v1x, double &v2x, double &v1c, double &v2c);

	// spin polarized GGA
	static void gcx_spin(double rhoup, double rhodw, double grhoup2, double grhodw2,
            double &sx, double &v1xup, double &v1xdw, double &v2xup,
            double &v2xdw);
	static void gcc_spin(double rho, double &zeta, double grho, double &sc,
            double &v1cup, double &v1cdw, double &v2c);

#ifdef USE_LIBXC
	struct Mgga_spin_in
	{
		double rhoup, rhodw;//electron densities
		ModuleBase::Vector3<double> grhoup, grhodw;//gradient of electron densities
		double tauup, taudw;//kinetic energy densities
	};

	struct Mgga_spin_out
	{
		double ex, ec;//xc energy densities
		double v1xup, v1xdw;//vx: lda part
		double v2xup, v2xdw;//vx: gga part
		double v3xup, v3xdw;//vx: mgga part
		double v1cup, v1cdw;//vc: lda part
		ModuleBase::Vector3<double> v2cup, v2cdw;
		std::vector<double> v2c;//vc: gga part, two different formats	
		double v3cup, v3cdw;//vc: mgga part
	};	
	
	// mGGA
	static void tau_xc(const double &rho, const double &grho, const double &atau, double &sx, double &sc,
          double &v1x, double &v2x, double &v3x, double &v1c, double &v2c, double &v3c);
	static void tau_xc_spin(const Mgga_spin_in &mgga_spin_in, Mgga_spin_out &mgga_spin_out);
#endif

	// PBEx, PBEc
	static void pbex(const double &rho, const double &grho, const int &iflag,
		double &sx, double &v1x, double &v2x);

	static void pbec(const double &rho, const double &grho, const int &flag,
		double &sc, double &v1c, double &v2c);



	private:

	// For LDA exchange energy
	static void slater(const double &rs, double &ex, double &vx);
	static void slater1(const double &rs, double &ex, double &vx);
	static void slater_rxc(const double &rs, double &ex, double &vx);

	// For LSDA exchange energy
	static void slater_spin( const double &rho, const double &zeta,
		double &ex, double &vxup, double &vxdw);
	static void slater1_spin( const double &rho, const double &zeta,
		double &ex, double &vxup, double &vxdw);
	static void slater_rxc_spin( const double &rho, const double &z,
		double &ex, double &vxup, double &vxdw);

	// For LDA correlation energy
	static void pz(const double &rs, const int &iflag, double &ec, double &vc);
	static void vwn(const double &rs, double &ec, double &vc);
	static void pw(const double &rs, const int &iflag, double &ec, double &vc);
	static void lyp(const double &rs, double &ec, double &vc);
	static void wigner(const double &rs, double &ec, double &vc);
	static void hl(const double &rs, double &ec, double &vc);
	static void gl(const double &rs, double &ec, double &vc);

	// For LSDA correlation energy
	// pz and pz_polarized should be put together
	static void pz_spin( const double &rs, const double &zeta,
    	double &ec, double &vcup, double &vcdw);
	static void pz_polarized( const double &rs, double &ec, double &vc);
	static void pw_spin( const double &rs, const double &zeta,
        double &ec, double &vcup, double &vcdw);

	// some GGA functionals
	static void hcth(const double rho, const double grho, double &sx, double &v1x, double &v2x);
	static void pwcorr(const double r, const double c[], double &g, double &dg);
	static void optx(const double rho, const double grho, double &sx, double &v1x, double &v2x);
	
	// For PW86 correlation functional
	static void perdew86(const double rho, const double grho, double &sc, double &v1c, double &v2c);

	// For PW91
	static void ggax(const double &rho, const double &grho, double &sx, double &v1x, double &v2x);
	static void ggac(const double &rho,const double &grho, double &sc, double &v1c, double &v2c);

	// For Wu-Cohen
	static void wcx(const double &rho,const double &grho, double &sx, double &v1x, double &v2x);

	// For BLYP
	static void becke88(const double &rho, const double &grho, double &sx, double &v1x, double &v2x);
	static void glyp(const double &rho, const double &grho, double &sc, double &v1c, double &v2c);

	// spin-polarized GGA
	static void becke88_spin(double rho, double grho, double &sx, double &v1x,
			double &v2x);
	static void perdew86_spin(double rho, double zeta, double grho, double &sc,
			double &v1cup, double &v1cdw, double &v2c);
	static void ggac_spin(double rho, double zeta, double grho, double &sc,
			double &v1cup, double &v1cdw, double &v2c);
	static void pbec_spin(double rho, double zeta, double grho, const int &flag, double &sc,
			double &v1cup, double &v1cdw, double &v2c);
};

#endif //XC_FUNCTION_H
