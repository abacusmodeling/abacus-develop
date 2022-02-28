#include "en_solver.h"
#include "FP/KSDFT/PW/ks_scf_pw.h"
#include "FP/KSDFT/LCAO/ks_scf_lcao.h"
#include "FP/OFDFT/ofdft.h"
#include "stdio.h"
namespace ModuleEnSover
{
void En_Solver:: printag()
{
    std::cout<<tag<<std::endl;
}


//Some API to operate E_Solver
void init_esolver(En_Solver *&p_ensolver, const string use_esol)
{
     if(use_esol == "ksdft_pw")
     {
         p_ensolver = new KS_SCF_PW();
     }
    else if(use_esol == "ksdft_lcao")
     {
         p_ensolver = new KS_SCF_LCAO();
     }
     //  else if(use_esol == "sdft_pw")
    //  {
    //      p_ensolver = new KS_SCF_PW(true);
    //  }
    //  else if(use_esol == "ofdft")
    //  {
    //      p_ensolver = new OFDFT();
    //  }    
}

void clean_esolver(En_Solver *&pesolver)
{
    if(pesolver!=NULL) 
    {
        delete pesolver;
    }
}

}