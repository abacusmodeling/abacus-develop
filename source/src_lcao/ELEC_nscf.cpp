#include "ELEC_nscf.h"
#include "../src_pw/global.h"
#include "../src_ri/exx_abfs.h"
#include "../src_ri/exx_opt_orb.h"
#include "global_fp.h"
#include "ELEC_cbands_gamma.h"
#include "ELEC_cbands_k.h"
#include "../src_io/berryphase.h"
#include "../src_io/to_wannier90.h"

ELEC_nscf::ELEC_nscf(){}
ELEC_nscf::~ELEC_nscf(){}

void ELEC_nscf::nscf(LCAO_Hamilt &uhm)
{
	TITLE("ELEC_nscf","nscf");

	cout << " NON-SELF CONSISTENT CALCULATIONS" << endl;
	
	time_t time_start= std::time(NULL);

	// Peize Lin add 2018-08-14
	switch(exx_lcao.info.hybrid_type)
	{
		case Exx_Global::Hybrid_Type::HF:
		case Exx_Global::Hybrid_Type::PBE0:
		case Exx_Global::Hybrid_Type::HSE:
			exx_lcao.cal_exx_elec_nscf();
			break;
	}

	// mohan add 2021-02-09
	// in LOOP_ions, istep starts from 1,
	// then when the istep is a variable of scf or nscf,
	// istep becomes istep-1, this should be fixed in future
	int istep=0; 
	if(GlobalV::GAMMA_ONLY_LOCAL)
	{
		ELEC_cbands_gamma::cal_bands(istep, uhm);
	}
	else
	{
		ELEC_cbands_k::cal_bands(istep, uhm);
	}

	time_t time_finish=std::time(NULL);
	OUT_TIME("cal_bands",time_start, time_finish);

    GlobalV::ofs_running << " end of band structure calculation " << endl;
    GlobalV::ofs_running << " band eigenvalue in this processor (eV) :" << endl;

    for (int ik = 0; ik < GlobalC::kv.nks; ik++)
    {
        if (GlobalV::NSPIN==2)
        {
            if (ik==0) 
			{
				GlobalV::ofs_running << " spin up :" << endl;
			}
            if (ik==( GlobalC::kv.nks / 2)) 
			{
				GlobalV::ofs_running << " spin down :" << endl;
			}
        }

		GlobalV::ofs_running << " k-points" 
			<< ik+1 << "(" << GlobalC::kv.nkstot << "): " 
			<< GlobalC::kv.kvec_c[ik].x << " " << GlobalC::kv.kvec_c[ik].y << " " << GlobalC::kv.kvec_c[ik].z << endl;

        for (int ib = 0; ib < GlobalV::NBANDS; ib++)
        {			
            GlobalV::ofs_running << " spin" << GlobalC::kv.isk[ik]+1 
			<< "final_state " << ib+1 << " " 
			<< GlobalC::wf.ekb[ik][ib] * Ry_to_eV 
			<< " " << GlobalC::wf.wg(ik, ib)*GlobalC::kv.nks << endl;
        }
		GlobalV::ofs_running << endl;
    }
	
	// add by jingan in 2018.11.7
	if(GlobalV::CALCULATION == "nscf" && INPUT.towannier90)
    {
        toWannier90 myWannier(GlobalC::kv.nkstot,ucell.G);
        myWannier.init_wannier();
    }
	
	// add by jingan
	if (berryphase::berry_phase_flag && Symmetry::symm_flag == 0)
    {
    	berryphase bp;
		bp.Macroscopic_polarization();
    }

	return;
}
