#include "esolver.h"
#include "FP/KSDFT/PW/ks_scf_pw.h"
#include "FP/OFDFT/ofdft.h"
#include "stdio.h"
namespace ModuleESolver
{
void ESolver:: printag()
{
    std::cout<<tag<<std::endl;
}


//Some API to operate E_Solver
void init_esolver(ESolver *&p_ensolver, const string use_esol)
{
     if(use_esol == "ksdft_pw")
     {
         p_ensolver = new KS_SCF_PW();
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

void clean_esolver(ESolver *&pesolver)
{
    if(pesolver!=NULL) 
    {
        delete pesolver;
    }
}

}