#include "esolver_ks.h"
#include <iostream>
#include <algorithm>
#include "time.h"

//--------------Temporary----------------
#include "../module_base/global_variable.h"
#include "../src_pw/global.h"
#include "../src_pw/charge_broyden.h"
#include "../module_base/timer.h"
//---------------------------------------

namespace ModuleESolver
{

ESolver_KS::ESolver_KS()
{
    basisname = "PLEASE ADD BASISNAME FOR CURRENT ESOLVER.";
    diag_ethr = GlobalV::PW_DIAG_THR; 
    scf_thr = GlobalV::SCF_THR; 
    drho = 0.0;
    maxniter = GlobalV::SCF_NMAX;
    niter = maxniter;


}

void ESolver_KS:: hamilt2density(int istep, int iter, double ethr)
{
    //Temporarily, before HSolver is constructed, it should be overrided by
    //LCAO, PW, SDFT and TDDFT.
    //After HSolver is constructed, LCAO, PW, SDFT should delete their own
    //hamilt2density() and use:
    //this->phsol->solve(this->phamilt, this->pes, this->wf, ETHR);
    
}


void ESolver_KS:: Run(int istep, UnitCell_pseudo& cell)
{
    ModuleBase::timer:: tick("ESolver_KS","run");
    
    this->printhead(); //print the headline on the screen.
    this->beforeiter(); //Something else to do before the iter loop
    
    bool firstscf = true;
    bool conv_elec = false;
    this->niter = this->maxniter;
    for(int iter=1; iter <= this->maxniter ; ++iter)
    {
        writehead(GlobalV::ofs_running, istep, iter); 
        clock_t iterstart,iterend;
        iterstart = std::clock();
        set_ethr(istep,iter);
        eachiterinit(iter); 
        
        this->hamilt2density(istep, iter, this->diag_ethr);
        
        // double drho = this->estate.caldr2(); 
        // EState should be used after it is constructed.
        drho = GlobalC::CHR.get_drho();

        if(firstscf)
        {
            firstscf = false;
            double hsover_error = this->diag_ethr * std::max(1.0, GlobalC::CHR.nelec);
            // The error of HSolver is larger than drho, so a more precise HSolver should be excuted.
            if(hsover_error > drho)  
            {
                GlobalV::ofs_running << " Notice: Threshold on eigenvalues was too large.\n";
                ModuleBase::WARNING("scf","Threshold on eigenvalues was too large.");
                GlobalV::ofs_running << " scf_thr=" << scf_thr << " < hsover_error=" << hsover_error << std::endl;

                GlobalV::ofs_running << " Origin diag_ethr = " << this->diag_ethr << std::endl;
                this->diag_ethr =0.1 *  drho / GlobalC::CHR.nelec; 
                GlobalV::ofs_running << " New    diag_ethr = " << this->diag_ethr << std::endl;
                
                this->hamilt2density(istep, iter, this->diag_ethr);
                drho = GlobalC::CHR.get_drho();
            }   
        }
        
        if(drho < this->scf_thr)   
            conv_elec = true;
        else
        {
            //charge mixing
            //conv_elec = this->estate.mix_rho();
            GlobalC::CHR.mix_rho(iter);
        }

        // Hamilt should be used after it is constructed.
        // this->phamilt->update(conv_elec);
        updatepot(conv_elec);
        eachiterfinish(iter);  
        iterend = std::clock();
        double duration = double(iterend-iterstart) / CLOCKS_PER_SEC;
        printiter(conv_elec, iter, drho, duration, diag_ethr);
        if (conv_elec) 
        {
            this->niter = iter;
            break;
        }  
    }
    afteriter(conv_elec); 
    
    ModuleBase::timer:: tick("ESolver_KS","run");
    return;
};


void ESolver_KS:: set_ethr(int istep, int iter)
{
//It is too complex now and should be modified.
    if(iter == 1)
    {
        if (abs(this->diag_ethr - 1.0e-2) < 1.0e-10)
        {
            if (GlobalC::pot.init_chg == "file")
            {
                //======================================================
                // if you think that the starting potential is good
                // do not spoil it with a louly first diagonalization:
                // set a strict this->diag_ethr in the input file ()diago_the_init
                //======================================================
                this->diag_ethr = 1.0e-5;
            }
            else
            {
                //=======================================================
                // starting atomic potential is probably far from scf
                // don't waste iterations in the first diagonalization
                //=======================================================
                this->diag_ethr = 1.0e-2;
            }
        }
        if(GlobalV::FINAL_SCF) this->diag_ethr = 1.0e-2;
    }
    else
    {
        if (iter == 2)
        {
            this->diag_ethr = 1.e-2;
        }

		//----------------------------
		// this->diag_ethr changes in CG Method.
		// mohan update 2012-03-26
		// mohan update 2012-02-08
		//----------------------------
		if(GlobalV::BASIS_TYPE=="lcao")
		{
			this->diag_ethr = std::min( this->diag_ethr, 0.01*this->drho/ std::max(1.0, GlobalC::CHR.nelec));
		}
		// mohan update 2009-09-04
		else
		{
			this->diag_ethr = std::min( this->diag_ethr, 0.1*this->drho/ std::max(1.0, GlobalC::CHR.nelec));
			//std::cout << " new this->diag_ethr = " << this->diag_ethr << std::endl;
		}

    }
}

void ESolver_KS:: printhead()
{
    std::cout << " " << std::setw(7) << "ITER";
    if(GlobalV::NSPIN==2)
    {
        std::cout << std::setw(10) << "TMAG";
        std::cout << std::setw(10) << "AMAG";
    }
    std::cout << std::setw(15) << "ETOT(eV)";
    std::cout << std::setw(15) << "EDIFF(eV)";
    std::cout << std::setw(11) << "SCF_THR";
	std::cout << std::setw(11) << "TIME(s)" << std::endl;
}

void ESolver_KS::printiter(bool conv_elec, int iter, double drho, double duration, double ethr)
{
    GlobalC::en.print_etot(conv_elec, iter, drho, duration, ethr);
}

void ESolver_KS:: writehead(std::ofstream &ofs_running, int istep, int iter)
{
    ofs_running
        << "\n "
        << this->basisname
        << " ALGORITHM --------------- ION=" << std::setw(4) << istep + 1
        << "  ELEC=" << std::setw(4) << iter
        << "--------------------------------\n";
}

int ESolver_KS:: getniter()
{
    return this->niter;
}


}
