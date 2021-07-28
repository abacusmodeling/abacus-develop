#include "../src_lcao/LOOP_ions.h"
#include "cal_r_overlap_R.h"
#include "../src_pw/global.h"
#include "write_HS.h"


void LOOP_ions::output_HS_R(void)
{
    TITLE("LOOP_ions","output_HS_R"); 
    timer::tick("LOOP_ions","output_HS_R"); 
	
	// add by jingan for out r_R matrix 2019.8.14
	if(INPUT.out_r_matrix)
	{
		cal_r_overlap_R r_matrix;
		r_matrix.init();
		r_matrix.out_r_overlap_R(GlobalV::NSPIN);
	}

    // Parameters for HR and SR output
    double sparse_threshold = 1e-10;
    bool binary = false; // output binary file

    if(GlobalV::NSPIN==1||GlobalV::NSPIN==4)
    {
        // GlobalC::UHM.calculate_STN_R();
        // GlobalC::UHM.GK.cal_vlocal_R(0);
        // GlobalC::UHM.GK.distribute_pvpR_tr();
        // HS_Matrix::save_HSR_tr(0);

        // jingan add 2021-6-4
        GlobalC::UHM.calculate_HSR_sparse(0, sparse_threshold);
        HS_Matrix::save_HSR_sparse(0, sparse_threshold, binary);
        GlobalC::UHM.destroy_all_HSR_sparse();
    }
    ///*
    else if(GlobalV::NSPIN==2)
    {
        // GlobalC::UHM.calculate_STN_R();
        // for(int ik=0; ik<GlobalC::kv.nks; ik++)
        // {
        //     if(ik==0 || ik==GlobalC::kv.nks/2)
        //     {
        //         if(GlobalV::NSPIN==2)GlobalV::CURRENT_SPIN = GlobalC::kv.isk[ik];
        //         for(int ir=0; ir<GlobalC::pw.nrxx; ir++)
        //         {
        //             GlobalC::pot.vr_eff1[ir] = GlobalC::pot.vr_eff( GlobalV::CURRENT_SPIN, ir);
        //         }
        	    	
        //         if(!GlobalV::GAMMA_ONLY_LOCAL)
        //         {
        //             if(GlobalV::VL_IN_H)
        //             {
		// 				//GlobalC::UHM.GK.cal_vlocal_k(GlobalC::pot.vrs1,GridT);
		// 				GlobalC::UHM.GK.cal_vlocal_k(GlobalC::pot.vr_eff1, GridT, GlobalV::CURRENT_SPIN);
        //             }
        //         }
        //         GlobalC::UHM.GK.cal_vlocal_R(GlobalV::CURRENT_SPIN);
        //         GlobalC::UHM.GK.distribute_pvpR_tr();
        //         HS_Matrix::save_HSR_tr(GlobalV::CURRENT_SPIN);
        //     }
        // }

        // jingan add 2021-6-4
        for(int ik = 0; ik < GlobalC::kv.nks; ik++)
        {
            if(ik == 0 || ik == GlobalC::kv.nks/2)
            {
                if(GlobalV::NSPIN == 2)
                {
                    GlobalV::CURRENT_SPIN = GlobalC::kv.isk[ik];
                }

                for(int ir = 0; ir < GlobalC::pw.nrxx; ir++)
                {
                    GlobalC::pot.vr_eff1[ir] = GlobalC::pot.vr_eff( GlobalV::CURRENT_SPIN, ir);
                }
        	    	
                if(!GlobalV::GAMMA_ONLY_LOCAL)
                {
                    if(GlobalV::VL_IN_H)
                    {
						//GlobalC::UHM.GK.cal_vlocal_k(GlobalC::pot.vrs1,GridT);
						GlobalC::UHM.GK.cal_vlocal_k(GlobalC::pot.vr_eff1, GridT, GlobalV::CURRENT_SPIN);
                    }
                }
                GlobalC::UHM.calculate_HSR_sparse(GlobalV::CURRENT_SPIN, sparse_threshold);
                HS_Matrix::save_HSR_sparse(GlobalV::CURRENT_SPIN, sparse_threshold, binary);
                GlobalC::UHM.destroy_all_HSR_sparse();
            }
        }
    }

    if(!GlobalV::GAMMA_ONLY_LOCAL) //LiuXh 20181011
    {
        GlobalC::UHM.GK.destroy_pvpR();
    } //LiuXh 20181011

    timer::tick("LOOP_ions","output_HS_R"); 
    return;
}
