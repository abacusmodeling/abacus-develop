#include "xc_functional.h"
#include "../src_parallel/parallel_reduce.h"
#include "../module_base/timer.h"

// [etxc, vtxc, v] = XC_Functional::v_xc(...)
std::tuple<double,double,ModuleBase::matrix> XC_Functional::v_xc
(
	const int &nrxx, // number of real-space grid
	const int &ncxyz, // total number of charge grid
	const double &omega, // volume of cell
    const double*const*const rho_in,
	const double*const rho_core) // core charge density
{
    ModuleBase::TITLE("H_XC_pw","v_xc");
    ModuleBase::timer::tick("H_XC_pw","v_xc");

	#ifndef USE_LIBXC
	if(GlobalV::DFT_META)
	{
		ModuleBase::WARNING_QUIT("Potential::v_of_rho","to use metaGGA, please link LIBXC");
	}
	#endif
    //Exchange-Correlation potential Vxc(r) from n(r)
    double etxc = 0.0;
    double vtxc = 0.0;
	ModuleBase::matrix v(GlobalV::NSPIN, nrxx);

    // the square of the e charge
    // in Rydeberg unit, so * 2.0.
    double e2 = 2.0;

    double rhox = 0.0;
    double arhox = 0.0;
    double zeta = 0.0;
    double exc = 0.0;
    double vxc[2];

    int ir, is;

    double vanishing_charge = 1.0e-10;

    if (GlobalV::NSPIN == 1 || ( GlobalV::NSPIN ==4 && !GlobalV::DOMAG && !GlobalV::DOMAG_Z))
    {
        // spin-unpolarized case
        for (int ir = 0;ir < nrxx;ir++)
        {
			// total electron charge density
            rhox = rho_in[0][ir] + rho_core[ir];
            arhox = abs(rhox);
            if (arhox > vanishing_charge)
            {
                if(use_libxc)
                {
                    XC_Functional::xc_libxc(arhox, exc, vxc[0]);
                }
                else
                {
                    XC_Functional::xc(arhox, exc, vxc[0]);
                }
                v(0,ir) = e2 * vxc[0];
				// consider the total charge density
                etxc += e2 * exc * rhox;
				// only consider rho_in
                vtxc += v(0, ir) * rho_in[0][ir];
            } // endif
        } //enddo
    }
    else if(GlobalV::NSPIN ==2)
    {
        // spin-polarized case
        for (ir = 0;ir < nrxx;ir++)
        {
            rhox = rho_in[0][ir] + rho_in[1][ir] + rho_core[ir]; //HLX(05-29-06): bug fixed
            arhox = abs(rhox);

            if (arhox > vanishing_charge)
            {
                zeta = (rho_in[0][ir] - rho_in[1][ir]) / arhox; //HLX(05-29-06): bug fixed

                if (abs(zeta)  > 1.0)
                {
                    zeta = (zeta > 0.0) ? 1.0 : (-1.0);
                }

                // call
                if(use_libxc)
                {
                    double rhoup = arhox * (1.0+zeta) / 2.0;
                    double rhodw = arhox * (1.0-zeta) / 2.0;
                    XC_Functional::xc_spin_libxc(rhoup, rhodw, exc, vxc[0], vxc[1]);
                }
                else
                {
                    double rhoup = arhox * (1.0+zeta) / 2.0;
                    double rhodw = arhox * (1.0-zeta) / 2.0;
                    XC_Functional::xc_spin(arhox, zeta, exc, vxc[0], vxc[1]);
                }

                for (is = 0;is < GlobalV::NSPIN;is++)
                {
                    v(is, ir) = e2 * vxc[is];
                }

                etxc += e2 * exc * rhox;
                vtxc += v(0, ir) * rho_in[0][ir] + v(1, ir) * rho_in[1][ir];
            }
        }

    }
    else if(GlobalV::NSPIN == 4)//noncollinear case added by zhengdy
    {
        for( ir = 0;ir<nrxx; ir++)
        {
            double amag = sqrt( pow(rho_in[1][ir],2) + pow(rho_in[2][ir],2) + pow(rho_in[3][ir],2) );

            rhox = rho_in[0][ir] + rho_core[ir];

            arhox = abs( rhox );

            if ( arhox > vanishing_charge )
            {
                zeta = amag / arhox;

                if ( abs( zeta ) > 1.0 )
                {
                    zeta = (zeta > 0.0) ? 1.0 : (-1.0);
                }//end if

                if(use_libxc)
                {
                    double rhoup = arhox * (1.0+zeta) / 2.0;
                    double rhodw = arhox * (1.0-zeta) / 2.0;
                    XC_Functional::xc_spin_libxc(rhoup, rhodw, exc, vxc[0], vxc[1]);
                }
                else
                {
                    double rhoup = arhox * (1.0+zeta) / 2.0;
                    double rhodw = arhox * (1.0-zeta) / 2.0;
                    XC_Functional::xc_spin(arhox, zeta, exc, vxc[0], vxc[1]);
                }

                etxc += e2 * exc * rhox;

                v(0, ir) = e2*( 0.5 * ( vxc[0] + vxc[1]) );
                vtxc += v(0,ir) * rho_in[0][ir];

                double vs = 0.5 * ( vxc[0] - vxc[1] );
                if ( amag > vanishing_charge )
                {
                    for(int ipol = 1;ipol< 4;ipol++)
                    {
                        v(ipol, ir) = e2 * vs * rho_in[ipol][ir] / amag;
                        vtxc += v(ipol,ir) * rho_in[ipol][ir];
                    }//end do
                }//end if
            }//end if
        }//end do
    }//end if
    // energy terms, local-density contributions

    // add gradient corrections (if any)
    // mohan modify 2009-12-15
    gradcorr(etxc, vtxc, v);

    // parallel code : collect vtxc,etxc
    // mohan add 2008-06-01
    Parallel_Reduce::reduce_double_pool( etxc );
    Parallel_Reduce::reduce_double_pool( vtxc );

    etxc *= omega / ncxyz;
    vtxc *= omega / ncxyz;

    ModuleBase::timer::tick("H_XC_pw","v_xc");
    return std::make_tuple(etxc, vtxc, std::move(v));
}
