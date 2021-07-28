#include "ELEC_cbands_k.h"
#include "LOOP_elec.h"
#include "LCAO_diago.h"
#include "../src_pw/global.h"
#include "../src_pw/symmetry_rho.h"
#include "LCAO_evolve.h"
#include "dftu.h"

ELEC_cbands_k::ELEC_cbands_k(){};
ELEC_cbands_k::~ELEC_cbands_k(){};


void ELEC_cbands_k::cal_bands(const int &istep, LCAO_Hamilt &uhm)
{
	TITLE("ELEC_cbands_k","cal_bands");
	timer::tick("ELEC_cbands_k","cal_bands");

	int start_spin = -1;
	uhm.GK.reset_spin(start_spin);
	uhm.GK.allocate_pvpR();
						
	// pool parallization in future -- mohan note 2021-02-09
	for(int ik=0; ik<GlobalC::kv.nks; ik++)
	{	
		//-----------------------------------------
		//(1) prepare data for this k point.
		// copy the local potential from array.
		//-----------------------------------------
		if(GlobalV::NSPIN==2) 
		{
			GlobalV::CURRENT_SPIN = GlobalC::kv.isk[ik];
		}
		GlobalC::wf.npw = GlobalC::kv.ngk[ik];
		for(int ir=0; ir<GlobalC::pw.nrxx; ir++)
		{
			GlobalC::pot.vr_eff1[ir] = GlobalC::pot.vr_eff( GlobalV::CURRENT_SPIN, ir);
		}
		
		//--------------------------------------------
		//(2) check if we need to calculate 
		// pvpR = < phi0 | v(spin) | phiR> for a new spin.
		//--------------------------------------------
		if(GlobalV::CURRENT_SPIN == uhm.GK.get_spin() )
		{
			//GlobalV::ofs_running << " Same spin, same vlocal integration." << endl;
		}
		else
		{
			//GlobalV::ofs_running << " (spin change)" << endl;
			uhm.GK.reset_spin( GlobalV::CURRENT_SPIN );

			// if you change the place of the following code,
			// rememeber to delete the #include	
			if(GlobalV::VL_IN_H)
			{
				// vlocal = Vh[rho] + Vxc[rho] + Vl(pseudo)
				uhm.GK.cal_vlocal_k(GlobalC::pot.vr_eff1,GridT);
				// added by zhengdy-soc, for non-collinear case
				// integral 4 times, is there any method to simplify?
				if(GlobalV::NSPIN==4)
				{
					for(int is=1;is<4;is++)
					{
						for(int ir=0; ir<GlobalC::pw.nrxx; ir++)
						{
							GlobalC::pot.vr_eff1[ir] = GlobalC::pot.vr_eff( is, ir);
						}
						uhm.GK.cal_vlocal_k(GlobalC::pot.vr_eff1, GridT, is);
					}
				}
			}
		}


		if(!uhm.init_s)
    	{
    	    WARNING_QUIT("Hamilt_Linear::solve_using_cg","Need init S matrix firstly");
    	}

		//--------------------------------------------
		// (3) folding matrix, 
		// and diagonalize the H matrix (T+Vl+Vnl).
		//--------------------------------------------

		// with k points
		timer::tick("Efficience","each_k");
		timer::tick("Efficience","H_k");
		uhm.calculate_Hk(ik);

		// Effective potential of DFT+U is added to total Hamiltonian here; Quxin adds on 20201029
		if(INPUT.dft_plus_u)
		{
      vector<complex<double>> eff_pot(GlobalC::ParaO.nloc);
			dftu.cal_eff_pot_mat_complex(ik, istep, &eff_pot[0]);
      
			for(int irc=0; irc<GlobalC::ParaO.nloc; irc++)
				GlobalC::LM.Hloc2[irc] += eff_pot[irc];					
		}

		timer::tick("Efficience","H_k");

		// Peize Lin add at 2020.04.04
		if(GlobalC::restart.info_load.load_H && !GlobalC::restart.info_load.load_H_finish)
		{
			GlobalC::restart.load_disk("H", ik);
			GlobalC::restart.info_load.load_H_finish = true;
		}			
		if(GlobalC::restart.info_save.save_H)
		{
			GlobalC::restart.save_disk("H", ik);
		}

		// write the wave functions into GlobalC::LOWF.WFC_K[ik].
		timer::tick("Efficience","diago_k");
		Diago_LCAO_Matrix DLM;
		DLM.solve_complex_matrix(ik, GlobalC::LOWF.WFC_K[ik], GlobalC::LOC.wfc_dm_2d.wfc_k[ik]);
		timer::tick("Efficience","diago_k");

		timer::tick("Efficience","each_k");
	} // end k
			
	// LiuXh modify 2019-07-15*/
	if(!GlobalC::ParaO.out_hsR)
	{
		uhm.GK.destroy_pvpR();
	}

	timer::tick("ELEC_cbands_k","cal_bands");
	return;	
}

