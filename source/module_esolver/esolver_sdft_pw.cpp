#include "./esolver_sdft_pw.h"
#include "time.h"
#include <fstream>
#include <algorithm>
#include "../module_base/timer.h"

//-------------------Temporary------------------
#include "../module_base/global_variable.h"
#include "../src_pw/global.h"
#include "../src_pw/symmetry_rho.h"
//----------------------------------------------
//-----force-------------------
#include "../src_pw/sto_forces.h"
//-----stress------------------
#include "../src_pw/sto_stress_pw.h"
//---------------------------------------------------

namespace ModuleESolver
{

ESolver_SDFT_PW::ESolver_SDFT_PW()
{
    classname = "ESolver_SDFT_PW";
    basisname = "PW";
}

ESolver_SDFT_PW::~ESolver_SDFT_PW()
{
}

void ESolver_SDFT_PW::Init(Input &inp, UnitCell_pseudo &cell)
{
    ESolver_KS_PW::Init(inp, cell);
	stowf.init(GlobalC::kv.nks);
	if(INPUT.nbands_sto != 0)	Init_Sto_Orbitals(this->stowf, INPUT.seed_sto);
	else						Init_Com_Orbitals(this->stowf, GlobalC::kv);
	for (int ik =0 ; ik < GlobalC::kv.nks; ++ik)
    {
        this->stowf.shchi[ik].create(this->stowf.nchip[ik],GlobalC::wf.npwx,false);
        if(GlobalV::NBANDS > 0)
        {
            this->stowf.chiortho[ik].create(this->stowf.nchip[ik],GlobalC::wf.npwx,false);
        }
    }
	stoiter.init(GlobalC::wf.npwx, this->stowf.nchip);
}

void ESolver_SDFT_PW::beforescf(const int istep)
{
    ESolver_KS_PW::beforescf(istep);
	if(istep > 0 && INPUT.nbands_sto != 0) Update_Sto_Orbitals(this->stowf, INPUT.seed_sto);
}

void ESolver_SDFT_PW::eachiterfinish(int iter)
{
	//print_eigenvalue(GlobalV::ofs_running);
    GlobalC::en.calculate_etot();
}
void ESolver_SDFT_PW::afterscf()
{
	for(int is=0; is<GlobalV::NSPIN; is++)
    {
        std::stringstream ssc;
        std::stringstream ss1;
        ssc << GlobalV::global_out_dir << "SPIN" << is + 1 << "_CHG";
		ss1 << GlobalV::global_out_dir << "SPIN" << is + 1 << "_CHG.cube";
        GlobalC::CHR.write_rho(GlobalC::CHR.rho_save[is], is, 0, ssc.str() );//mohan add 2007-10-17
	    GlobalC::CHR.write_rho_cube(GlobalC::CHR.rho_save[is], is, ss1.str(), 3);
    }
    if(this->conv_elec)
    {
        //GlobalV::ofs_running << " convergence is achieved" << std::endl;			
        //GlobalV::ofs_running << " !FINAL_ETOT_IS " << GlobalC::en.etot * ModuleBase::Ry_to_eV << " eV" << std::endl; 
        GlobalV::ofs_running << "\n charge density convergence is achieved" << std::endl;
        GlobalV::ofs_running << " final etot is " << GlobalC::en.etot * ModuleBase::Ry_to_eV << " eV" << std::endl;
    }
    else
    {
        GlobalV::ofs_running << " convergence has NOT been achieved!" << std::endl;
    }
}

void ESolver_SDFT_PW::hamilt2density(int istep, int iter, double ethr)
{
    double *h_diag = new double[GlobalC::wf.npwx * GlobalV::NPOL];
	GlobalV::ofs_running << " "  <<setw(8) << "K-point" << setw(15) << "CG iter num" << setw(15) << "Time(Sec)"<< std::endl;
	GlobalV::ofs_running << setprecision(6) << setiosflags(ios::fixed) << setiosflags(ios::showpoint);
	for (int ik = 0;ik < GlobalC::kv.nks;ik++)
	{
		if(GlobalV::NBANDS > 0 && GlobalV::MY_STOGROUP == 0)
		{
			this->c_bands_k(ik,h_diag,istep+1,iter);
		}
		else
		{
			GlobalC::hm.hpw.init_k(ik); 
			//In fact, hm.hpw.init_k has been done in wf.wfcinit();
		}
		stoiter.stohchi.current_ik = ik;
		
#ifdef __MPI
			if(GlobalV::NBANDS > 0)
			{
				MPI_Bcast(GlobalC::wf.evc[ik].c, GlobalC::wf.npwx*GlobalV::NBANDS*2, MPI_DOUBLE , 0, PARAPW_WORLD);
				MPI_Bcast(GlobalC::wf.ekb[ik], GlobalV::NBANDS, MPI_DOUBLE, 0, PARAPW_WORLD);
			}
#endif
			stoiter.orthog(ik,this->stowf);
			stoiter.checkemm(ik,iter,this->stowf);	//check and reset emax & emin
		}

		for (int ik = 0;ik < GlobalC::kv.nks;ik++)
		{
			//init k
			if(GlobalC::kv.nks > 1) GlobalC::hm.hpw.init_k(ik);
			stoiter.stohchi.current_ik = ik;
			stoiter.sumpolyval_k(ik, this->stowf);
		}
		delete [] h_diag;
		GlobalC::en.eband  = 0.0;
        GlobalC::en.demet  = 0.0;
        GlobalC::en.ef     = 0.0;
        GlobalC::en.ef_up  = 0.0;
        GlobalC::en.ef_dw  = 0.0;
		stoiter.itermu(iter);
		//(5) calculate new charge density 
		// calculate KS rho.
		if(GlobalV::NBANDS > 0)
		{
			if(GlobalV::MY_STOGROUP == 0)
			{
				GlobalC::CHR.sum_band();
			}
			else
			{
				for(int is=0; is < GlobalV::NSPIN; is++)
				{
					ModuleBase::GlobalFunc::ZEROS(GlobalC::CHR.rho[is], GlobalC::rhopw->nrxx);
				}
			}
			MPI_Bcast(&GlobalC::en.eband,1, MPI_DOUBLE, 0,PARAPW_WORLD);
		}
		else
		{
			for(int is=0; is < GlobalV::NSPIN; is++)
			{
				ModuleBase::GlobalFunc::ZEROS(GlobalC::CHR.rho[is], GlobalC::rhopw->nrxx);
			}
		}
	// calculate stochastic rho
		stoiter.sum_stoband(this->stowf);
		

		//(6) calculate the delta_harris energy 
		// according to new charge density.
		// mohan add 2009-01-23
		//en.calculate_harris(2);

		if(GlobalV::MY_STOGROUP==0)
		{
			Symmetry_rho srho;
			for(int is=0; is < GlobalV::NSPIN; is++)
			{
				srho.begin(is, GlobalC::CHR,GlobalC::rhopw, GlobalC::Pgrid, GlobalC::symm);
			}
		}
		else
		{
			if(ModuleSymmetry::Symmetry::symm_flag)	MPI_Barrier(MPI_COMM_WORLD);
		}

		if(GlobalV::MY_STOGROUP == 0)
		{
        	GlobalC::en.deband = GlobalC::en.delta_e();
		}

}

void ESolver_SDFT_PW:: c_bands_k(const int ik, double* h_diag, const int istep, const int iter)
{
	ModuleBase::timer::tick(this->classname,"c_bands_k");
	int precondition_type = 2;
	GlobalC::hm.hpw.init_k(ik);
	
    //===========================================
    // Conjugate-Gradient diagonalization
    // h_diag is the precondition matrix
    // h_diag(1:npw) = MAX( 1.0, g2kin(1:npw) );
    //===========================================
    if (precondition_type==1)
    {
        for (int ig = 0;ig < GlobalC::wf.npw; ++ig)
		{
			double g2kin = GlobalC::wfcpw->getgk2(ik,ig) * GlobalC::wfcpw->tpiba2;
			h_diag[ig] = std::max(1.0, g2kin);
			if(GlobalV::NPOL==2) h_diag[ig+GlobalC::wf.npwx] = h_diag[ig];
		}
    }
    else if (precondition_type==2)
    {
        for (int ig = 0;ig < GlobalC::wf.npw; ig++)
		{
			double g2kin = GlobalC::wfcpw->getgk2(ik,ig) * GlobalC::wfcpw->tpiba2;
			h_diag[ig] = 1 + g2kin + sqrt( 1 + (g2kin - 1) * (g2kin - 1));
			if(GlobalV::NPOL==2) h_diag[ig+GlobalC::wf.npwx] = h_diag[ig];
		}
    }
	//h_diag can't be zero!  //zhengdy-soc
	if(GlobalV::NPOL==2)
	{
		for(int ig = GlobalC::wf.npw;ig < GlobalC::wf.npwx; ig++)
		{
			h_diag[ig] = 1.0;
			h_diag[ig+ GlobalC::wf.npwx] = 1.0;
		}
	}
	clock_t start=clock();

	//============================================================
	// diago the hamiltonian!!
	// In plane wave method, firstly using cinitcgg to diagnolize,
	// then using cg method.
	//
	// In localized orbital presented in plane wave case,
	// only using cinitcgg.
	//
	// In linear scaling method, using sparse matrix and
	// adjacent searching code and cg method to calculate the
	// eigenstates.
	//=============================================================
	double avg_iter_k = 0.0;
	GlobalC::hm.diagH_pw(istep, iter, ik, h_diag, avg_iter_k);

	GlobalC::en.print_band(ik);
	clock_t finish=clock();
	const double duration = static_cast<double>(finish - start) / CLOCKS_PER_SEC;
	GlobalV::ofs_running << " " << setw(8) 
		<< ik+1 << setw(15) 
		<< avg_iter_k << setw(15) << duration << endl;

	ModuleBase::timer::tick(this->classname,"c_bands_k");
}


void ESolver_SDFT_PW::cal_Energy(energy &en)
{
	
}

void ESolver_SDFT_PW::cal_Force(ModuleBase::matrix &force)
{
	Sto_Forces ff;
    ff.init(force, this->stowf);
}
void ESolver_SDFT_PW::cal_Stress(ModuleBase::matrix &stress)
{
	Sto_Stress_PW ss;
    ss.cal_stress(stress, this->stowf);
}


}