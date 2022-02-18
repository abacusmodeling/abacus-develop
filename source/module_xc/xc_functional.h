//==========================================================
// AUTHOR : mohan
// DATE : 2009-04-04
//==========================================================
#ifndef XC_FUNCTIONAL_H
#define XC_FUNCTIONAL_H

#ifdef USE_LIBXC
#include <xc.h>
#else
#include "xc_funcs.h"
#endif	// ifdef USE_LIBXC

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/vector3.h"
#include "../module_base/matrix.h"
#include "exx_global.h"
class XC_Functional
{
	public:

	friend class Run_lcao;
	friend class Run_pw;
	friend class Ions;
	friend class ELEC_scf;
	friend class LOOP_elec;
	friend class Run_MD_PW;

	XC_Functional();
	~XC_Functional();

	// compute the exchange-correlation energy 
	// [etxc, vtxc, v] = v_xc(...)
    static std::tuple<double,double,ModuleBase::matrix> v_xc(
		const int &nrxx, // number of real-space grid
		const int &ncxyz, // total number of charge grid
		const double &omega, // volume of cell
		const double*const*const rho_in, 
		const double*const rho_core); // core charge density

	static std::tuple<double,double,ModuleBase::matrix,ModuleBase::matrix> v_xc_meta(
		const double * const * const rho_in,
		const double * const rho_core_in,
		const double * const * const kin_r_in);

	// GGA
	static void gcxc(const double &rho, const double &grho,
			double &sxc, double &v1xc, double &v2xc);

	static void gcxc_libxc(const double &rho, const double &grho,
			double &sxc, double &v1xc, double &v2xc);

	// spin polarized GGA
	static void gcx_spin(double rhoup, double rhodw, double grhoup2, double grhodw2,
            double &sx, double &v1xup, double &v1xdw, double &v2xup,
            double &v2xdw);
	static void gcc_spin(double rho, double &zeta, double grho, double &sc,
            double &v1cup, double &v1cdw, double &v2c);

	static void gcxc_spin_libxc(double rhoup, double rhodw, 
		ModuleBase::Vector3<double> gdr1, ModuleBase::Vector3<double> gdr2,
        double &sxc, double &v1xcup, double &v1xcdw, double &v2xcup, double &v2xcdw, double &v2xcud);

	static void gradcorr(double &etxc, double &vtxc, ModuleBase::matrix &v);
	static void grad_rho( const std::complex<double> *rhog, ModuleBase::Vector3<double> *gdr );
	static void grad_wfc( const std::complex<double> *rhog, const int ik, std::complex<double> **grad, const int npw );
	static void grad_dot( const ModuleBase::Vector3<double> *h, double *dh);
	static void noncolin_rho(double *rhoout1,double *rhoout2,double *seg);

	// mGGA
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
	static void tau_xc(const double &rho, const double &grho, const double &atau, double &sx, double &sc,
          double &v1x, double &v2x, double &v3x, double &v1c, double &v2c, double &v3c);
	static void tau_xc_spin(const Mgga_spin_in &mgga_spin_in, Mgga_spin_out &mgga_spin_out);
#endif 

	static int get_func_type();

	private:

	static std::vector<int> func_id; // id of exchange functional
	static int func_type; //0:none, 1:lda, 2:gga, 3:mgga, 4:hybrid
	static bool use_libxc;

	static void set_xc_type(const std::string xc_func_in);
#ifdef USE_LIBXC
	static std::vector<xc_func_type> init_func(const int xc_polarized);
#endif
	// LDA
	static void xc(const double &rho, double &exc, double &vxc);
	static void xc_libxc(const double &rho, double &exc, double &vxc);


	// LSDA
	static void xc_spin(const double &rho, const double &zeta,
			double &exc, double &vxcup, double &vxcdw);
	static void xc_spin_libxc(const double &rhoup, const double &rhodw,
			double &exc, double &vxcup, double &vxcdw);

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

	// PBEx, PBEc
	static void pbex(const double &rho, const double &grho, const int &iflag,
		double &sx, double &v1x, double &v2x);

	static void pbec(const double &rho, const double &grho, const int &flag,
		double &sc, double &v1c, double &v2c);

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

	//----------------------------
	// decide the value of spin
	//----------------------------
	static int nspin0() // may need updates from SOC
	{
		if     (GlobalV::NSPIN==1 || (GlobalV::NSPIN==4 && (!GlobalV::DOMAG && !GlobalV::DOMAG_Z)))		return 1;
		else if(GlobalV::NSPIN==2 || (GlobalV::NSPIN==4 && ( GlobalV::DOMAG ||  GlobalV::DOMAG_Z)))		return 2;
		else throw std::runtime_error(ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
	}
};

#endif //XC_FUNCTION_H
