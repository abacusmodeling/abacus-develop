#include "wavefunc.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_hamilt_lcao/hamilt_lcaodft/wavefunc_in_pw.h"
#include "module_io/winput.h"
#include "module_base/memory.h"
#include "module_base/timer.h"
#include "module_hsolver/diago_iter_assist.h"
#include "module_hamilt_pw/hamilt_pwdft/hamilt_pw.h"
#include "module_hsolver/kernels/math_kernel_op.h"

wavefunc::wavefunc()
{
	out_wfc_pw = 0;
}

wavefunc::~wavefunc()
{
	if(GlobalV::test_deconstructor)
	{
		std::cout << " ~wavefunc()" << std::endl;
	}
	if(this->irindex == nullptr) 
	{
		delete[] this->irindex;		
		this->irindex=nullptr;
	}
}

psi::Psi<std::complex<double>>* wavefunc::allocate(const int nks)
{	
	ModuleBase::TITLE("wavefunc","allocate");

	this->npwx = GlobalC::wfcpw->npwk_max;
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"npwx",npwx);

	assert(npwx > 0);
	assert(nks > 0);

	// allocate for kinetic energy

	// if use spin orbital, do not double nks but double allocate evc and wanf2.
	int prefactor = 1;
	if(GlobalV::NSPIN==4) prefactor = GlobalV::NPOL;//added by zhengdy-soc

	const int nks2 = nks;

	psi::Psi<std::complex<double>>* psi_out = nullptr;
	if(GlobalV::CALCULATION=="nscf" && GlobalC::wf.mem_saver==1)
	{
		//initial psi rather than evc
		psi_out = new psi::Psi<std::complex<double>>(1, GlobalV::NBANDS, npwx * GlobalV::NPOL, GlobalC::kv.ngk.data());
		if(GlobalV::BASIS_TYPE=="lcao_in_pw")
		{
			wanf2[0].create(GlobalV::NLOCAL, npwx * GlobalV::NPOL);
			const size_t memory_cost = GlobalV::NLOCAL*(GlobalV::NPOL*npwx) * sizeof(std::complex<double>);
			std::cout << " Memory for wanf2 (MB): " << double(memory_cost)/1024.0/1024.0 << std::endl;
			ModuleBase::Memory::record("WF::wanf2", memory_cost) ;
		}
		const size_t memory_cost = GlobalV::NBANDS*(GlobalV::NPOL*npwx) * sizeof(std::complex<double>);
		std::cout << " MEMORY FOR PSI (MB)  : " << double(memory_cost)/1024.0/1024.0 << std::endl;
		ModuleBase::Memory::record("Psi_PW", memory_cost);
	}
	else if(GlobalV::BASIS_TYPE!="pw")
	{
		this->evc = new ModuleBase::ComplexMatrix [nks2];
		for (int ik = 0; ik < nks2; ik++)
		{
			this->evc[ik].create(GlobalV::NBANDS, npwx * GlobalV::NPOL);
		};
		const size_t memory_cost = nks2 * GlobalV::NBANDS*(GlobalV::NPOL*npwx) * sizeof(std::complex<double>);
		std::cout << " MEMORY FOR PSI (MB)  : " << double(memory_cost)/1024.0/1024.0 << std::endl;
		ModuleBase::Memory::record("Psi_PW", memory_cost);

		if((GlobalV::BASIS_TYPE=="lcao" || GlobalV::BASIS_TYPE=="lcao_in_pw") || winput::out_spillage==2)
		{//for lcao_in_pw
			this->wanf2 = new ModuleBase::ComplexMatrix [nks2];
			for (int ik = 0; ik < nks2; ik++)
			{
				this->wanf2[ik].create(GlobalV::NLOCAL, npwx * GlobalV::NPOL);
			}
			const size_t memory_cost = nks2 * GlobalV::NLOCAL*(npwx * GlobalV::NPOL) * sizeof(std::complex<double>);
			std::cout << " Memory for wanf2 (MB): " << double(memory_cost)/1024.0/1024.0 << std::endl;
			ModuleBase::Memory::record("WF::wanf2", memory_cost) ;
		}
	}
	else
	{
		//initial psi rather than evc
		psi_out = new psi::Psi<std::complex<double>>(nks2, GlobalV::NBANDS, npwx * GlobalV::NPOL, GlobalC::kv.ngk.data());
		const size_t memory_cost = nks2 * GlobalV::NBANDS*(GlobalV::NPOL*npwx) * sizeof(std::complex<double>);
		std::cout << " MEMORY FOR PSI (MB)  : " << double(memory_cost)/1024.0/1024.0 << std::endl;
		ModuleBase::Memory::record("Psi_PW", memory_cost);
	}
	return psi_out;

	//showMemStats();
}

//===================================================================
// This routine computes an estimate of the start_ wavefunctions
// from superposition of atomic wavefunctions or random wave functions.
//===================================================================
void wavefunc::wfcinit(psi::Psi<std::complex<double>>* psi_in)
{
    ModuleBase::TITLE("wavefunc","wfcinit");
    ModuleBase::timer::tick("wavefunc","wfcinit");

    this->wfcinit_k(psi_in);

    GlobalC::en.demet = 0.0;

    ModuleBase::timer::tick("wavefunc","wfcinit");
    return;
}

int wavefunc::get_starting_nw(void)const
{
    if (init_wfc == "file")
    {
		throw std::runtime_error("wavefunc::get_starting_nw. start_ wfc from file: not implemented yet! "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__)); 	// Peize Lin change 2019-05-01
        //ModuleBase::WARNING_QUIT("wfcinit_k","\n start_ wfc from file: not implemented yet!");
        //**********************************************************************
        // ... read the wavefunction into memory (if it is not done in c_bands)
        //**********************************************************************
    }
    else if (init_wfc.substr(0,6) == "atomic")
    {
        if (GlobalC::ucell.natomwfc >= GlobalV::NBANDS)
        {
            if(GlobalV::test_wf)GlobalV::ofs_running << " Start wave functions are all pseudo atomic wave functions." << std::endl;
        }
        else
        {
            if(GlobalV::test_wf)GlobalV::ofs_running << " Start wave functions are atomic + "
            << GlobalV::NBANDS - GlobalC::ucell.natomwfc
            << " random wave functions." << std::endl;
        }
        return max(GlobalC::ucell.natomwfc,  GlobalV::NBANDS);
    }
    else if (init_wfc == "random")
    {
        if(GlobalV::test_wf)GlobalV::ofs_running << " Start wave functions are all random." << std::endl;
        return GlobalV::NBANDS;
    }
    else
    {
		throw std::runtime_error("wavefunc::get_starting_nw. Don't know what to do! Please Check source code! "+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__)); 	// Peize Lin change 2019-05-01
        //ModuleBase::WARNING_QUIT("get_starting_nw","Don't know what to do! Please Check source code!");
    }
}

namespace hamilt
{

void diago_PAO_in_pw_k2(const int &ik, psi::Psi<std::complex<float>> &wvf, hamilt::Hamilt<float>* phm_in)
{
    ModuleBase::TITLE("wavefunc","diago_PAO_in_pw_k2");
    // (6) Prepare for atmoic orbitals or random orbitals
    const int starting_nw = GlobalC::wf.get_starting_nw();
    if(starting_nw == 0) return;
    assert(starting_nw > 0);
	std::vector<float> etatom(starting_nw, 0.0);

    const int nbasis = wvf.get_nbasis();
    const int nbands = wvf.get_nbands();
    const int current_nbasis = GlobalC::kv.ngk[ik];

    //special case here! use Psi(k-1) for the initialization of Psi(k)
    //this method should be tested.
    /*if(GlobalV::CALCULATION == "nscf" && GlobalC::ucell.natomwfc == 0 && ik>0)
    {
        //this is memsaver case
        if(wvf.get_nk() == 1)
        {
            return;
        }
        else
        {
            ModuleBase::GlobalFunc::COPYARRAY(&wvf(ik-1, 0, 0), &wvf(ik, 0, 0), wvf.get_nbasis()* wvf.get_nbands());
            return;
        }
    }
    */

    if( GlobalC::wf.init_wfc=="random" || ( GlobalC::wf.init_wfc.substr(0,6)=="atomic" && GlobalC::ucell.natomwfc == 0 ))
	{
		GlobalC::wf.random(wvf.get_pointer(),0,nbands,ik, GlobalC::wfcpw);

		if(GlobalV::KS_SOLVER=="cg") //xiaohui add 2013-09-02
		{
			if(phm_in!= nullptr)
			{
				hsolver::DiagoIterAssist<float>::diagH_subspace(phm_in, wvf, wvf, etatom.data());
				return;
			}
			else
			{
				ModuleBase::WARNING_QUIT("wavefunc","Hamiltonian does not exist!");
			}
		}
	}
	else if(GlobalC::wf.init_wfc.substr(0,6)=="atomic")
	{
		ModuleBase::ComplexMatrix wfcatom(starting_nw, nbasis);//added by zhengdy-soc
		if(GlobalV::test_wf)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "starting_nw", starting_nw);

		GlobalC::wf.atomic_wfc(ik, current_nbasis, GlobalC::ucell.lmax_ppwf, wfcatom, GlobalC::ppcell.tab_at, GlobalV::NQX, GlobalV::DQ);
		if( GlobalC::wf.init_wfc == "atomic+random" && starting_nw == GlobalC::ucell.natomwfc )//added by qianrui 2021-5-16
		{
			GlobalC::wf.atomicrandom(wfcatom,0,starting_nw,ik, GlobalC::wfcpw);
		}

		//====================================================
		// If not enough atomic wfc are available, complete
		// with random wfcs
		//====================================================
		GlobalC::wf.random(wfcatom.c, GlobalC::ucell.natomwfc, nbands, ik, GlobalC::wfcpw);

		// (7) Diago with cg method.
		std::vector<std::complex<float>> s_wfcatom(starting_nw * nbasis);
		castmem_z2c_h2h_op()(cpu_ctx, cpu_ctx, s_wfcatom.data(), wfcatom.c, starting_nw * nbasis);
		//if(GlobalV::DIAGO_TYPE == "cg") xiaohui modify 2013-09-02
		if(GlobalV::KS_SOLVER=="cg") //xiaohui add 2013-09-02
		{
			if(phm_in!= nullptr)
			{
				// hsolver::DiagoIterAssist<double>::diagH_subspace_init(phm_in,
				//                          wfcatom,
				//                          wvf,
				//                          etatom.data());
				hsolver::DiagoIterAssist<float>::diagH_subspace_init(phm_in,
																	s_wfcatom.data(),
																	wfcatom.nr,
																	wfcatom.nc,
																	wvf,
																	etatom.data());
				return;
			}
			else
			{
				ModuleBase::WARNING_QUIT("wavefunc","Psi does not exist!");
				//this diagonalization method is obsoleted now
				//GlobalC::hm.diagH_subspace(ik ,starting_nw, nbands, wfcatom, wfcatom, etatom.data());
			}
		}

		assert(nbands <= wfcatom.nr);
		for (int ib=0; ib<nbands; ib++)
		{
			for (int ig=0; ig<nbasis; ig++)
			{
				wvf(ib, ig) = s_wfcatom[ib * nbasis + ig];
			}
		}
	}
}

void diago_PAO_in_pw_k2(const int &ik, psi::Psi<std::complex<double>> &wvf, hamilt::Hamilt<double>* phm_in)
{
	ModuleBase::TITLE("wavefunc","diago_PAO_in_pw_k2");
	// (6) Prepare for atmoic orbitals or random orbitals
	const int starting_nw = GlobalC::wf.get_starting_nw();
	if(starting_nw == 0) return;
	assert(starting_nw > 0);
	std::vector<double> etatom(starting_nw, 0.0);

	const int nbasis = wvf.get_nbasis();
	const int nbands = wvf.get_nbands();
	const int current_nbasis = GlobalC::kv.ngk[ik];

	//special case here! use Psi(k-1) for the initialization of Psi(k)
	//this method should be tested.
	/*if(GlobalV::CALCULATION == "nscf" && GlobalC::ucell.natomwfc == 0 && ik>0)
	{
		//this is memsaver case
		if(wvf.get_nk() == 1)
		{
			return;
		}
		else
		{
			ModuleBase::GlobalFunc::COPYARRAY(&wvf(ik-1, 0, 0), &wvf(ik, 0, 0), wvf.get_nbasis()* wvf.get_nbands());
			return;
		}
	}
	*/

	if( GlobalC::wf.init_wfc=="random" || ( GlobalC::wf.init_wfc.substr(0,6)=="atomic" && GlobalC::ucell.natomwfc == 0 ))
	{
		GlobalC::wf.random(wvf.get_pointer(),0,nbands,ik, GlobalC::wfcpw);

		if(GlobalV::KS_SOLVER=="cg") //xiaohui add 2013-09-02
		{
			if(phm_in!= nullptr)
			{
				hsolver::DiagoIterAssist<double>::diagH_subspace(phm_in, wvf, wvf, etatom.data());
				return;
			}
			else
			{
				ModuleBase::WARNING_QUIT("wavefunc","Hamiltonian does not exist!");
			}
		}
	}
	else if(GlobalC::wf.init_wfc.substr(0,6)=="atomic")
	{
		ModuleBase::ComplexMatrix wfcatom(starting_nw, nbasis);//added by zhengdy-soc
		if(GlobalV::test_wf)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "starting_nw", starting_nw);
		
		GlobalC::wf.atomic_wfc(ik, current_nbasis, GlobalC::ucell.lmax_ppwf, wfcatom, GlobalC::ppcell.tab_at, GlobalV::NQX, GlobalV::DQ);
		if( GlobalC::wf.init_wfc == "atomic+random" && starting_nw == GlobalC::ucell.natomwfc )//added by qianrui 2021-5-16
		{
			GlobalC::wf.atomicrandom(wfcatom,0,starting_nw,ik, GlobalC::wfcpw);
		}

		//====================================================
		// If not enough atomic wfc are available, complete
		// with random wfcs
		//====================================================
		GlobalC::wf.random(wfcatom.c, GlobalC::ucell.natomwfc, nbands, ik, GlobalC::wfcpw);

		// (7) Diago with cg method.
		//if(GlobalV::DIAGO_TYPE == "cg") xiaohui modify 2013-09-02
		if(GlobalV::KS_SOLVER=="cg") //xiaohui add 2013-09-02
		{
			if(phm_in!= nullptr)
			{
				hsolver::DiagoIterAssist<double>::diagH_subspace_init(phm_in,
							wfcatom.c,
							wfcatom.nr,
							wfcatom.nc,
							wvf,
							etatom.data());
				return;
			}
			else
			{
				ModuleBase::WARNING_QUIT("wavefunc","Hamiltonian does not exist!");
			}
		}

		assert(nbands <= wfcatom.nr);
		for (int ib=0; ib<nbands; ib++)
		{
			for (int ig=0; ig<nbasis; ig++)
			{
				wvf(ib, ig) = wfcatom(ib, ig);
			}
		}
	} //end of atomic case
}

template <>
void diago_PAO_in_pw_k2(const psi::DEVICE_CPU* ctx, const int &ik, psi::Psi<std::complex<float>, psi::DEVICE_CPU> &wvf, hamilt::Hamilt<float, psi::DEVICE_CPU>* phm_in)
{
    diago_PAO_in_pw_k2(ik, wvf, phm_in);
}
template <>
void diago_PAO_in_pw_k2(const psi::DEVICE_CPU* ctx, const int &ik, psi::Psi<std::complex<double>, psi::DEVICE_CPU> &wvf, hamilt::Hamilt<double, psi::DEVICE_CPU>* phm_in)
{
	diago_PAO_in_pw_k2(ik, wvf, phm_in);
}

#if ((defined __CUDA) || (defined __ROCM))
template<>
void diago_PAO_in_pw_k2(const psi::DEVICE_GPU* ctx, const int &ik, psi::Psi<std::complex<float>, psi::DEVICE_GPU> &wvf, hamilt::Hamilt<float, psi::DEVICE_GPU>* phm_in)
{
    ModuleBase::TITLE("wavefunc","diago_PAO_in_pw_k2");
    // (6) Prepare for atmoic orbitals or random orbitals
    const int starting_nw = GlobalC::wf.get_starting_nw();
    if(starting_nw == 0) return;
    assert(starting_nw > 0);

    const int nbasis = wvf.get_nbasis();
    const int nbands = wvf.get_nbands();
    const int current_nbasis = GlobalC::kv.ngk[ik];

    ModuleBase::ComplexMatrix wfcatom(starting_nw, nbasis);//added by zhengdy-soc
    if(GlobalV::test_wf)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "starting_nw", starting_nw);
    if(GlobalC::wf.init_wfc.substr(0,6)=="atomic")
    {
        GlobalC::wf.atomic_wfc(ik, current_nbasis, GlobalC::ucell.lmax_ppwf, wfcatom, GlobalC::ppcell.tab_at, GlobalV::NQX, GlobalV::DQ);
        if( GlobalC::wf.init_wfc == "atomic+random" && starting_nw == GlobalC::ucell.natomwfc )//added by qianrui 2021-5-16
        {
            GlobalC::wf.atomicrandom(wfcatom,0,starting_nw,ik, GlobalC::wfcpw);
        }

        //====================================================
        // If not enough atomic wfc are available, complete
        // with random wfcs
        //====================================================
        GlobalC::wf.random(wfcatom.c, GlobalC::ucell.natomwfc, nbands, ik, GlobalC::wfcpw);
    }
    else if(GlobalC::wf.init_wfc=="random")
    {
        GlobalC::wf.random(wfcatom.c,0,nbands,ik, GlobalC::wfcpw);
    }

    // store wfcatom on the GPU
    std::complex<float> * c_wfcatom = nullptr;
    resmem_cd_op()(gpu_ctx, c_wfcatom, wfcatom.nr * wfcatom.nc);
    castmem_z2c_h2d_op()(gpu_ctx, cpu_ctx, c_wfcatom, wfcatom.c, wfcatom.nr * wfcatom.nc);

    if(GlobalV::KS_SOLVER=="cg") //xiaohui add 2013-09-02
    {
        // (7) Diago with cg method.
        if(phm_in!= nullptr)
        {
            std::vector<float> etatom(starting_nw, 0.0);
            hsolver::DiagoIterAssist<float, psi::DEVICE_GPU>::diagH_subspace_init(
                    phm_in,
                    c_wfcatom,
                    wfcatom.nr,
                    wfcatom.nc,
                    wvf,
                    etatom.data());
        }
        else
        {
            //this diagonalization method is obsoleted now
            //GlobalC::hm.diagH_subspace(ik ,starting_nw, nbands, wfcatom, wfcatom, etatom.data());
        }
    }
    else if(GlobalV::KS_SOLVER=="dav")
    {
        assert(nbands <= wfcatom.nr);
        // replace by haozhihan 2022-11-23
        hsolver::matrixSetToAnother<float, psi::DEVICE_GPU>()(
                gpu_ctx,
                nbands,
                c_wfcatom,
                wfcatom.nc,
                &wvf(0,0),
                nbasis
        );
    }
    delmem_cd_op()(gpu_ctx, c_wfcatom);
}
template<>
void diago_PAO_in_pw_k2(const psi::DEVICE_GPU* ctx, const int &ik, psi::Psi<std::complex<double>, psi::DEVICE_GPU> &wvf, hamilt::Hamilt<double, psi::DEVICE_GPU>* phm_in)
{
	ModuleBase::TITLE("wavefunc","diago_PAO_in_pw_k2");
	// (6) Prepare for atmoic orbitals or random orbitals
	const int starting_nw = GlobalC::wf.get_starting_nw();
	if(starting_nw == 0) return;
	assert(starting_nw > 0);

	const int nbasis = wvf.get_nbasis();
	const int nbands = wvf.get_nbands();
	const int current_nbasis = GlobalC::kv.ngk[ik];

	ModuleBase::ComplexMatrix wfcatom(starting_nw, nbasis);//added by zhengdy-soc
	if(GlobalV::test_wf)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "starting_nw", starting_nw);
	if(GlobalC::wf.init_wfc.substr(0,6)=="atomic")
	{
		GlobalC::wf.atomic_wfc(ik, current_nbasis, GlobalC::ucell.lmax_ppwf, wfcatom, GlobalC::ppcell.tab_at, GlobalV::NQX, GlobalV::DQ);
		if( GlobalC::wf.init_wfc == "atomic+random" && starting_nw == GlobalC::ucell.natomwfc )//added by qianrui 2021-5-16
		{
			GlobalC::wf.atomicrandom(wfcatom,0,starting_nw,ik, GlobalC::wfcpw);
		}

		//====================================================
		// If not enough atomic wfc are available, complete
		// with random wfcs
		//====================================================
		GlobalC::wf.random(wfcatom.c, GlobalC::ucell.natomwfc, nbands, ik, GlobalC::wfcpw);
	}
	else if(GlobalC::wf.init_wfc=="random")
	{
		GlobalC::wf.random(wfcatom.c,0,nbands,ik, GlobalC::wfcpw);
	}

	// store wfcatom on the GPU
	std::complex<double> *z_wfcatom = nullptr;
	resmem_zd_op()(gpu_ctx, z_wfcatom, wfcatom.nr * wfcatom.nc);
	syncmem_z2z_h2d_op()(gpu_ctx, cpu_ctx, z_wfcatom, wfcatom.c, wfcatom.nr * wfcatom.nc);

	if(GlobalV::KS_SOLVER=="cg") //xiaohui add 2013-09-02
	{
		// (7) Diago with cg method.
		if(phm_in!= nullptr)
		{
			std::vector<double> etatom(starting_nw, 0.0);
			hsolver::DiagoIterAssist<double, psi::DEVICE_GPU>::diagH_subspace_init(
					     phm_in,
                         z_wfcatom,
						 wfcatom.nr,
						 wfcatom.nc,
                         wvf,
                         etatom.data());
		}
		else
		{
			//this diagonalization method is obsoleted now
			//GlobalC::hm.diagH_subspace(ik ,starting_nw, nbands, wfcatom, wfcatom, etatom.data());
		}
	}
	else if(GlobalV::KS_SOLVER=="dav")
	{
		assert(nbands <= wfcatom.nr);
		// replace by haozhihan 2022-11-23
		hsolver::matrixSetToAnother<double, psi::DEVICE_GPU>()(
			gpu_ctx,
			nbands,
			z_wfcatom,
			wfcatom.nc,
			&wvf(0,0),
			nbasis
		);
	}
	delmem_zd_op()(gpu_ctx, z_wfcatom);
}
#endif

}//namespace hamilt

void wavefunc::wfcinit_k(psi::Psi<std::complex<double>>* psi_in)
{
	ModuleBase::TITLE("wavefunc","wfcinit_k");

	if(GlobalV::BASIS_TYPE=="pw") 
	{
		this->irindex = new int [GlobalC::wfcpw->fftnxy];
		GlobalC::wfcpw->getfftixy2is(this->irindex);
    #if defined(__CUDA) || defined(__ROCM)
    if (GlobalV::device_flag == "gpu") {
        GlobalC::wfcpw->get_ig2ixyz_k();
    }
    #endif
	}
	if(GlobalV::CALCULATION=="nscf")
	{
		return;
	}
	return;
}

//--------------------------------------------
// get the nearest unitcell positions
// that exist overlaps between two orbitals
// iw1 and iw2
//--------------------------------------------
int wavefunc::get_R(int ix, int iy, int iz)   // pengfei 2016-11-23
{
	int count;
	ModuleBase::Vector3<double> r,r1,r2;

	for(int iw1=0; iw1<GlobalV::NLOCAL; iw1++)
	{
		for(int iw2=0; iw2<GlobalV::NLOCAL; iw2++)
		{
			int it1 = iw2it(iw1); int ia1 = iw2ia(iw1);
			int it2 = iw2it(iw2); int ia2 = iw2ia(iw2);
			//std::cout <<"iw1= "<<iw1<<" iw2= "<<iw2<<" it1= "<<it1<<" ia1= "<<ia1<<" it2= "<<it2<<" ia2= "<<ia2<<std::endl;
			count = 0;

			for(int nx=-int(ix/2);nx<=int(ix/2);nx++)
			{
				for(int ny=-int(iy/2);ny<=int(iy/2);ny++)
				{
					for(int nz=-int(iz/2);nz<=int(iz/2);nz++)
					{
						//std::cout <<"count = "<<count<<std::endl;
						//std::cout<<"nx= "<<nx<<" ny= "<<ny<<" nz= "<<nz<<std::endl;
						r1.x = GlobalC::ucell.atoms[it1].tau[ia1].x * GlobalC::ucell.lat0;
						r1.y = GlobalC::ucell.atoms[it1].tau[ia1].y * GlobalC::ucell.lat0;
						r1.z = GlobalC::ucell.atoms[it1].tau[ia1].z * GlobalC::ucell.lat0;
						r2.x = (GlobalC::ucell.atoms[it2].tau[ia2].x
						+ GlobalC::ucell.latvec.e11 * nx + GlobalC::ucell.latvec.e21 * ny + GlobalC::ucell.latvec.e31 * nz) * GlobalC::ucell.lat0;
						r2.y = (GlobalC::ucell.atoms[it2].tau[ia2].y
						+ GlobalC::ucell.latvec.e12 * nx + GlobalC::ucell.latvec.e22 * ny + GlobalC::ucell.latvec.e32 * nz) * GlobalC::ucell.lat0;
						r2.z = (GlobalC::ucell.atoms[it2].tau[ia2].z
						+ GlobalC::ucell.latvec.e13 * nx + GlobalC::ucell.latvec.e23 * ny + GlobalC::ucell.latvec.e33 * nz) * GlobalC::ucell.lat0;
						r = r2 - r1;
						double distance = sqrt(r*r);

						if(distance < (GlobalC::ucell.atoms[it1].Rcut + GlobalC::ucell.atoms[it2].Rcut))
						{
							R[iw1][iw2][count].x = nx;
							R[iw1][iw2][count].y = ny;
							R[iw1][iw2][count].z = nz;
							count++;
						}
					}
				}
			}
			Rmax[iw1][iw2] = count;
		}
	}

	int NR = 0;
	for(int iw1=0;iw1<GlobalV::NLOCAL;iw1++)
	{
		for(int iw2=0;iw2<GlobalV::NLOCAL;iw2++)
		{
			if(Rmax[iw1][iw2] > NR)
			{
				NR = Rmax[iw1][iw2];
			}
		}
	}

	return NR;
}


int wavefunc::iw2it( int iw)    // pengfei 2016-11-23
{
    int ic, type;
    ic =0;
    for(int it =0; it<GlobalC::ucell.ntype; it++)
	{
        for(int ia = 0; ia<GlobalC::ucell.atoms[it].na; ia++)
        {
            for(int L=0; L<GlobalC::ucell.atoms[it].nwl+1; L++)
			{
                for(int N=0; N<GlobalC::ucell.atoms[it].l_nchi[L]; N++)
                {
                    for(int i=0; i<(2*L+1); i++)
                    {
                        if(ic == iw)
                        {
                           type = it;
                        }
                        ic++;
					}
                }
			}
        }
	}
    return type;
}

int wavefunc::iw2ia( int iw)    // pengfei 2016-11-23
{
	int ic, na;
	ic =0;
	for(int it =0; it<GlobalC::ucell.ntype; it++)
	{
		for(int ia = 0; ia<GlobalC::ucell.atoms[it].na; ia++)
		{
			for(int L=0; L<GlobalC::ucell.atoms[it].nwl+1; L++)
				for(int N=0; N<GlobalC::ucell.atoms[it].l_nchi[L]; N++)
				{
					for(int i=0; i<(2*L+1); i++)
					{
						if(ic == iw)
						{
							na = ia;
						}
						ic++;
					}
				}
		}
	}
	return na;
}

//LiuXh add a new function here,
//20180515
void wavefunc::init_after_vc(const int nks)
{
    ModuleBase::TITLE("wavefunc","init");
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"npwx",this->npwx);

    assert(this->npwx > 0);
    assert(nks > 0);
    assert(GlobalV::NBANDS > 0);

    const int nks2 = nks;
	const int nbasis = this->npwx * GlobalV::NPOL;

	if((GlobalV::BASIS_TYPE=="lcao" || GlobalV::BASIS_TYPE=="lcao_in_pw") || winput::out_spillage==2)
	{
		if(wanf2 != nullptr)delete[] wanf2;
		this->wanf2 = new ModuleBase::ComplexMatrix [nks2];
		for (int ik = 0; ik < nks2; ik++)
		{
			this->wanf2[ik].create(GlobalV::NLOCAL, nbasis);
		}
	}

    if(GlobalV::test_wf)
    {
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"psi allocation","Done");
        if(GlobalV::BASIS_TYPE=="lcao" || GlobalV::BASIS_TYPE=="lcao_in_pw")
        {
            ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"wanf2 allocation","Done");
        }
    }

    return;
}