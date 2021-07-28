#include "FORCE_gamma.h"
#include "../src_pw/global.h"

void Force_LCAO_gamma::cal_fvl_dphi(
	matrix& dm2d, 
	const bool isforce, 
	const bool isstress, 
	matrix& fvl_dphi, 
	matrix& svl_dphi)
{   
    TITLE("Force_LCAO_gamma","cal_fvl_dphi");
    timer::tick("Force_LCAO_gamma","cal_fvl_dphi");
	
	if(isforce)
	{
		ZEROS (GlobalC::LM.DHloc_fixed_x, GlobalC::ParaO.nloc);
		ZEROS (GlobalC::LM.DHloc_fixed_y, GlobalC::ParaO.nloc);
		ZEROS (GlobalC::LM.DHloc_fixed_z, GlobalC::ParaO.nloc);
	}
    if(isstress)
    {
        ZEROS (GlobalC::LM.DHloc_fixed_11, GlobalC::ParaO.nloc);
        ZEROS (GlobalC::LM.DHloc_fixed_12, GlobalC::ParaO.nloc);
        ZEROS (GlobalC::LM.DHloc_fixed_13, GlobalC::ParaO.nloc);
        ZEROS (GlobalC::LM.DHloc_fixed_22, GlobalC::ParaO.nloc);
        ZEROS (GlobalC::LM.DHloc_fixed_23, GlobalC::ParaO.nloc);
        ZEROS (GlobalC::LM.DHloc_fixed_33, GlobalC::ParaO.nloc);
    }


    double* tmpDHx = new double[GlobalC::ParaO.nloc];
    double* tmpDHy = new double[GlobalC::ParaO.nloc];
    double* tmpDHz = new double[GlobalC::ParaO.nloc];
    ZEROS( tmpDHx, GlobalC::ParaO.nloc );
    ZEROS( tmpDHy, GlobalC::ParaO.nloc );
    ZEROS( tmpDHz, GlobalC::ParaO.nloc );
    for(int i=0; i<GlobalC::ParaO.nloc; ++i)
    {
        tmpDHx[i] = GlobalC::LM.DHloc_fixed_x[i];
        tmpDHy[i] = GlobalC::LM.DHloc_fixed_y[i];
        tmpDHz[i] = GlobalC::LM.DHloc_fixed_z[i];
        //cout << "  GlobalC::LM.DHloc_fixed_x=" <<  GlobalC::LM.DHloc_fixed_x[i] << endl;
        //cout << "  GlobalC::LM.DHloc_fixed_y=" <<  GlobalC::LM.DHloc_fixed_y[i] << endl;
        //cout << "  GlobalC::LM.DHloc_fixed_z=" <<  GlobalC::LM.DHloc_fixed_z[i] << endl;
    }

    //calculate dVL
    //calculate <dphi | VL | phi>

	// mohan add 2021, needs reconstruction!!!
    int istep = 1;
    GlobalC::pot.init_pot(istep, GlobalC::pw.strucFac);

    for(int is=0; is<GlobalV::NSPIN; ++is)
    {
        GlobalV::CURRENT_SPIN = is;

        ZEROS (GlobalC::LM.DHloc_fixed_x, GlobalC::ParaO.nloc);
        ZEROS (GlobalC::LM.DHloc_fixed_y, GlobalC::ParaO.nloc);
        ZEROS (GlobalC::LM.DHloc_fixed_z, GlobalC::ParaO.nloc);

        for(int ir=0; ir<GlobalC::pw.nrxx; ++ir)
        {
            GlobalC::pot.vr_eff1[ir] = GlobalC::pot.vr_eff(GlobalV::CURRENT_SPIN, ir);
        }

        //  should not be set zero if VNA is used.
        //  ZEROS(GlobalC::LM.DHloc_fixed_x,GlobalC::ParaO.nloc);
        //  ZEROS(GlobalC::LM.DHloc_fixed_y,GlobalC::ParaO.nloc);
        //  ZEROS(GlobalC::LM.DHloc_fixed_z,GlobalC::ParaO.nloc);
        GlobalC::UHM.GG.cal_force(GlobalC::pot.vr_eff1);


        for(int i=0; i<GlobalV::NLOCAL; i++)
        {
            const int iat = GlobalC::ucell.iwt2iat[i];
            for(int j=0; j<GlobalV::NLOCAL; j++)
            {
                //const int iat2 = GlobalC::ucell.iwt2iat[j];
                const int mu = GlobalC::ParaO.trace_loc_row[j];
                const int nu = GlobalC::ParaO.trace_loc_col[i];
                if (mu >= 0 && nu >= 0 )
                {
                    const int index = mu * GlobalC::ParaO.ncol + nu;
                    //contribution from deriv of AO's in T+VNL term

                    double dm2d2 = 2.0 * dm2d(is, index);

					if(isforce)
					{
						fvl_dphi(iat,0) -= dm2d2 * ( GlobalC::LM.DHloc_fixed_x[index] + tmpDHx[index] );
						fvl_dphi(iat,1) -= dm2d2 * ( GlobalC::LM.DHloc_fixed_y[index] + tmpDHy[index] );
						fvl_dphi(iat,2) -= dm2d2 * ( GlobalC::LM.DHloc_fixed_z[index] + tmpDHz[index] );
					}
                    if(isstress)
                    {
                        svl_dphi(0,0) += dm2d2 * GlobalC::LM.DHloc_fixed_11[index];
                        svl_dphi(0,1) += dm2d2 * GlobalC::LM.DHloc_fixed_12[index];
                        svl_dphi(0,2) += dm2d2 * GlobalC::LM.DHloc_fixed_13[index];
                        svl_dphi(1,1) += dm2d2 * GlobalC::LM.DHloc_fixed_22[index];
                        svl_dphi(1,2) += dm2d2 * GlobalC::LM.DHloc_fixed_23[index];
                        svl_dphi(2,2) += dm2d2 * GlobalC::LM.DHloc_fixed_33[index];
                    }
                    //  cout << setw(5) << iat << setw(5) << iat2 
                    //  << setw(5) << mu << setw(5) << nu
                    //  << setw(15) << GlobalC::LM.DHloc_fixed_z[index] << endl;
                }
            }
        }

//          cout << "fvl_dphi:" << endl;
//          for(int iat=0; iat<GlobalC::ucell.nat; ++iat)
//          {
//              cout << setw(5) << iat << setw(15) << fvl_dphi[iat][0] 
//              << setw(15) << fvl_dphi[iat][1]
//              << setw(15) << fvl_dphi[iat][2] << endl;
//          }


    } // end spin
    // test mohan tmp
//  test_gamma(GlobalC::LM.DHloc_fixed_x,"GlobalC::LM.DHloc_fixed_x");

    if(isstress)
    {
        for(int i=0;i<3;i++)
        {
            for(int j=0;j<3;j++)
            {
                if(i<j) svl_dphi(j,i) = svl_dphi(i,j);
				svl_dphi(i,j) /= GlobalC::ucell.omega;
            }
        }
    }

    delete[] tmpDHx;
    delete[] tmpDHy;
    delete[] tmpDHz;


    timer::tick("Force_LCAO_gamma","cal_fvl_dphi");
    return;
}



void Force_LCAO_gamma::cal_fvl_dphi(
	const std::vector<matrix> &dm2d, 
	const bool isforce, 
	const bool isstress, 
	matrix& fvl_dphi, 
	matrix& svl_dphi)
{   
    TITLE("Force_LCAO_gamma","cal_fvl_dphi");
    timer::tick("Force_LCAO_gamma","cal_fvl_dphi");

    ZEROS (GlobalC::LM.DHloc_fixed_x, GlobalC::ParaO.nloc);
    ZEROS (GlobalC::LM.DHloc_fixed_y, GlobalC::ParaO.nloc);
    ZEROS (GlobalC::LM.DHloc_fixed_z, GlobalC::ParaO.nloc);
    if(GlobalV::STRESS)
    {
        ZEROS (GlobalC::LM.DHloc_fixed_11, GlobalC::ParaO.nloc);
        ZEROS (GlobalC::LM.DHloc_fixed_12, GlobalC::ParaO.nloc);
        ZEROS (GlobalC::LM.DHloc_fixed_13, GlobalC::ParaO.nloc);
        ZEROS (GlobalC::LM.DHloc_fixed_22, GlobalC::ParaO.nloc);
        ZEROS (GlobalC::LM.DHloc_fixed_23, GlobalC::ParaO.nloc);
        ZEROS (GlobalC::LM.DHloc_fixed_33, GlobalC::ParaO.nloc);
    }


    double* tmpDHx = new double[GlobalC::ParaO.nloc];
    double* tmpDHy = new double[GlobalC::ParaO.nloc];
    double* tmpDHz = new double[GlobalC::ParaO.nloc];
    ZEROS( tmpDHx, GlobalC::ParaO.nloc );
    ZEROS( tmpDHy, GlobalC::ParaO.nloc );
    ZEROS( tmpDHz, GlobalC::ParaO.nloc );
    for(int i=0; i<GlobalC::ParaO.nloc; ++i)
    {
        tmpDHx[i] = GlobalC::LM.DHloc_fixed_x[i];
        tmpDHy[i] = GlobalC::LM.DHloc_fixed_y[i];
        tmpDHz[i] = GlobalC::LM.DHloc_fixed_z[i];
        //cout << "  GlobalC::LM.DHloc_fixed_x=" <<  GlobalC::LM.DHloc_fixed_x[i] << endl;
        //cout << "  GlobalC::LM.DHloc_fixed_y=" <<  GlobalC::LM.DHloc_fixed_y[i] << endl;
        //cout << "  GlobalC::LM.DHloc_fixed_z=" <<  GlobalC::LM.DHloc_fixed_z[i] << endl;
    }

    //calculate dVL
    //calculate <dphi | VL | phi>

    int istep = 1;
    GlobalC::pot.init_pot(istep, GlobalC::pw.strucFac);

    for(int is=0; is<GlobalV::NSPIN; ++is)
    {
        GlobalV::CURRENT_SPIN = is;

        ZEROS (GlobalC::LM.DHloc_fixed_x, GlobalC::ParaO.nloc);
        ZEROS (GlobalC::LM.DHloc_fixed_y, GlobalC::ParaO.nloc);
        ZEROS (GlobalC::LM.DHloc_fixed_z, GlobalC::ParaO.nloc);

        for(int ir=0; ir<GlobalC::pw.nrxx; ++ir)
        {
            GlobalC::pot.vr_eff1[ir] = GlobalC::pot.vr_eff(GlobalV::CURRENT_SPIN, ir);
        }

        //should not be set zero if VNA is used.
        //ZEROS(GlobalC::LM.DHloc_fixed_x,GlobalC::ParaO.nloc);
        //ZEROS(GlobalC::LM.DHloc_fixed_y,GlobalC::ParaO.nloc);
        //ZEROS(GlobalC::LM.DHloc_fixed_z,GlobalC::ParaO.nloc);
        GlobalC::UHM.GG.cal_force(GlobalC::pot.vr_eff1);

        for(int i=0; i<GlobalV::NLOCAL; i++)
        {
            const int iat = GlobalC::ucell.iwt2iat[i];
            for(int j=0; j<GlobalV::NLOCAL; j++)
            {
                //const int iat2 = GlobalC::ucell.iwt2iat[j];
                const int mu = GlobalC::ParaO.trace_loc_row[j];
                const int nu = GlobalC::ParaO.trace_loc_col[i];
                if (mu >= 0 && nu >= 0 )
                {
                    const int index = mu * GlobalC::ParaO.ncol + nu;
                    //contribution from deriv of AO's in T+VNL term

                    double dm2d2 = 2.0 * dm2d[is](nu, mu);

                    if(isforce)
					{
						fvl_dphi(iat,0) -= dm2d2 * ( GlobalC::LM.DHloc_fixed_x[index] + tmpDHx[index] );
						fvl_dphi(iat,1) -= dm2d2 * ( GlobalC::LM.DHloc_fixed_y[index] + tmpDHy[index] );
						fvl_dphi(iat,2) -= dm2d2 * ( GlobalC::LM.DHloc_fixed_z[index] + tmpDHz[index] );
					}
                    if(isstress)
                    {
                        svl_dphi(0,0) += dm2d2 * GlobalC::LM.DHloc_fixed_11[index];
                        svl_dphi(0,1) += dm2d2 * GlobalC::LM.DHloc_fixed_12[index];
                        svl_dphi(0,2) += dm2d2 * GlobalC::LM.DHloc_fixed_13[index];
                        svl_dphi(1,1) += dm2d2 * GlobalC::LM.DHloc_fixed_22[index];
                        svl_dphi(1,2) += dm2d2 * GlobalC::LM.DHloc_fixed_23[index];
                        svl_dphi(2,2) += dm2d2 * GlobalC::LM.DHloc_fixed_33[index];
                    }
                    //cout << setw(5) << iat << setw(5) << iat2 
                    //<< setw(5) << mu << setw(5) << nu
                    //<< setw(15) << GlobalC::LM.DHloc_fixed_z[index] << endl;
                }
            }
        }

        //cout << "fvl_dphi:" << endl;
        //for(int iat=0; iat<GlobalC::ucell.nat; ++iat)
        //{
        //cout << setw(5) << iat << setw(15) << fvl_dphi[iat][0] 
        //<< setw(15) << fvl_dphi[iat][1]
        //<< setw(15) << fvl_dphi[iat][2] << endl;
        //}

    } // end spin
    // test mohan tmp
    //test_gamma(GlobalC::LM.DHloc_fixed_x,"GlobalC::LM.DHloc_fixed_x");

    if(isstress)
    {
        for(int i=0;i<3;i++)
        {
            for(int j=0;j<3;j++)
            {
                if(i<j) svl_dphi(j,i) = svl_dphi(i,j);
				svl_dphi(i,j) /= GlobalC::ucell.omega;
            }
        }
    }

    delete[] tmpDHx;
    delete[] tmpDHy;
    delete[] tmpDHz;

    timer::tick("Force_LCAO_gamma", "cal_fvl_dphi");
    return;
}
