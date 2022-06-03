#include "esolver_fp.h"
#include "../module_base/global_variable.h"
#include "../src_pw/global.h"
namespace ModuleESolver
{   ESolver_FP::ESolver_FP()
    {
        // pw_rho = new ModuleBase::PW_Basis();
        
        pw_rho = new ModulePW::PW_Basis_Big(); 
        GlobalC::rhopw = this->pw_rho; //Temporary
        //temporary, it will be removed
        GlobalC::bigpw = static_cast<ModulePW::PW_Basis_Big*>(pw_rho);
        GlobalC::bigpw->setbxyz(INPUT.bx,INPUT.by,INPUT.bz);
    }
    ESolver_FP::~ESolver_FP()
    {
        if(pw_rho!=NULL)    delete pw_rho;
    }
    void ESolver_FP::Init(Input& inp, UnitCell_pseudo& cell)
    {
        // Initalize the plane wave basis set
        std::cout << " UNIFORM GRID DIM     : " << GlobalC::rhopw->nx << " * " << GlobalC::rhopw->ny << " * " << GlobalC::rhopw->nz << std::endl;
        std::cout << " UNIFORM GRID DIM(BIG): " << GlobalC::bigpw->nbx << " * " << GlobalC::bigpw->nby << " * " << GlobalC::bigpw->nbz << std::endl;
        if (inp.nx * inp.ny * inp.nz == 0)
            this->pw_rho->initgrids(cell.lat0, cell.latvec, inp.ecutrho, GlobalV::NPROC_IN_POOL, GlobalV::RANK_IN_POOL);
	    else
            this->pw_rho->initgrids(cell.lat0, cell.latvec, inp.nx, inp.ny, inp.nz, GlobalV::NPROC_IN_POOL, GlobalV::RANK_IN_POOL);
        this->pw_rho->initparameters(false, inp.ecutrho);
        this->pw_rho->setuptransform();
        this->pw_rho->collect_local_pw(); 
        this->pw_rho->collect_uniqgg();
        this->print_rhofft(inp, GlobalV::ofs_running);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"number of |g|", this->pw_rho->ngg);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"max |g|", this->pw_rho->gg_uniq[ this->pw_rho->ngg-1]);
	    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"min |g|", this->pw_rho->gg_uniq[0]);
    }

    void ESolver_FP::print_rhofft(Input&inp, ofstream &ofs)
    {
        ofs << "\n\n\n\n";
	    ofs << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
	    ofs << " |                                                                    |" << std::endl;
	    ofs << " | Setup plane waves:                                                 |" << std::endl;
	    ofs << " | Use the energy cutoff and the lattice vectors to generate the      |" << std::endl;
	    ofs << " | dimensions of FFT grid. The number of FFT grid on each processor   |" << std::endl;
	    ofs << " | is 'nrxx'. The number of plane wave basis in reciprocal space is   |" << std::endl;
	    ofs << " | different for charege/potential and wave functions. We also set    |" << std::endl;
	    ofs << " | the 'sticks' for the parallel of FFT.                              |" << std::endl;
	    ofs << " |                                                                    |" << std::endl;
	    ofs << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
	    ofs << "\n\n\n\n";
	    ofs << "\n SETUP THE PLANE WAVE BASIS" << std::endl;
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"energy cutoff for wavefunc (unit:Ry)",INPUT.ecutwfc);
        if(inp.nx * inp.ny * inp.nz == 0)
            GlobalV::ofs_running << " use input fft dimensions for wave functions and charge." << std::endl;
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"fft grid for wave functions", this->pw_rho->nx,this->pw_rho->ny,this->pw_rho->nz);
	    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"fft grid for charge/potential", this->pw_rho->nx,this->pw_rho->ny,this->pw_rho->nz);
	    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"fft grid division",GlobalC::bigpw->bx,GlobalC::bigpw->by,GlobalC::bigpw->bz);
	    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"big fft grid for charge/potential",GlobalC::bigpw->nbx,GlobalC::bigpw->nby,GlobalC::bigpw->nbz);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"nbxx",GlobalC::bigpw->nbxx);
	    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"nrxx",this->pw_rho->nrxx);

        GlobalV::ofs_running << "\n SETUP PLANE WAVES FOR CHARGE/POTENTIAL" << std::endl;
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"number of plane waves",this->pw_rho->npwtot);
	    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"number of sticks", this->pw_rho->nstot);

        GlobalV::ofs_running << "\n PARALLEL PW FOR CHARGE/POTENTIAL" << std::endl;
        GlobalV::ofs_running <<" "<< std::setw(8)  << "PROC"<< std::setw(15) << "COLUMNS(POT)"<< std::setw(15) << "PW" << std::endl;
        for (int i = 0; i < GlobalV::NPROC_IN_POOL ; ++i)
        {
            GlobalV::ofs_running <<" "<<std::setw(8)<< i+1 << std::setw(15) << this->pw_rho->nst_per[i] << std::setw(15) << this->pw_rho->npw_per[i] << std::endl;
        }
        GlobalV::ofs_running << " --------------- sum -------------------" << std::endl;
        GlobalV::ofs_running << " " << std::setw(8)  << GlobalV::NPROC_IN_POOL << std::setw(15) << this->pw_rho->nstot << std::setw(15) << this->pw_rho->npwtot << std::endl;

    }
}