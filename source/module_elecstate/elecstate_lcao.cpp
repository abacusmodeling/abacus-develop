#include "elecstate_lcao.h"
#include "math_tools.h"
#include "module_base/timer.h"
#include "src_lcao/grid_technique.h"

namespace elecstate
{

//for gamma_only(double case) and multi-k(complex<double> case)
template<typename T> void ElecStateLCAO::cal_dm(const ModuleBase::matrix& wg,
    const psi::Psi<T>& wfc,
    psi::Psi<T>& dm)
{
	ModuleBase::TITLE("ElecStateLCAO", "cal_dm");
	
    dm.resize( wfc.get_nk(), this->loc->ParaV->ncol, this->loc->ParaV->nrow );
    const int nbands_local = wfc.get_nbands();
    const int nbasis_local = wfc.get_nbasis();

	// dm = wfc.T * wg * wfc.conj()
	// dm[is](iw1,iw2) = \sum_{ib} wfc[is](ib,iw1).T * wg(is,ib) * wfc[is](ib,iw2).conj()
    for(int ik=0; ik<wfc.get_nk(); ++ik)
    {
        wfc.fix_k(ik);
        dm.fix_k(ik);
        // wg_wfc(ib,iw) = wg[ib] * wfc(ib,iw);
        psi::Psi<T> wg_wfc( wfc, 1 );

        int ib_global = 0;
        for(int ib_local=0; ib_local<nbands_local; ++ib_local)
        {
            while(ib_local != this->loc->ParaV->trace_loc_col[ib_global])
            {
                ++ib_global;
                if(ib_global>=wg.nc)
                {
                    ModuleBase::WARNING_QUIT("ElecStateLCAO::cal_dm", "please check trace_loc_col!");
                }
            }
            const double wg_local = wg(ik, ib_global);
            T* wg_wfc_pointer = &(wg_wfc(0, ib_local, 0));
            BlasConnector::scal( nbasis_local, wg_local, wg_wfc_pointer, 1 );
        }

        // C++: dm(iw1,iw2) = wfc(ib,iw1).T * wg_wfc(ib,iw2)
#ifdef __MPI
        psiMulPsiMpi(wg_wfc, wfc, dm, this->loc->ParaV->desc_wfc, this->loc->ParaV->desc);
#else 
        psiMulPsi(wg_wfc, wfc, dm);
#endif
    }

	return;
}

void ElecStateLCAO::psiToRho(const psi::Psi<std::complex<double>>& psi)
{
    ModuleBase::TITLE("ElecStateLCAO", "psiToRho");
    ModuleBase::timer::tick("ElecStateLCAO", "psiToRho");

    psi::Psi<std::complex<double>> dm_k_2d(psi.get_nk(), psi.get_nbasis(), psi.get_nbasis());

    ModuleBase::GlobalFunc::NOTE("Calculate the density matrix.");
    //this->loc->cal_dk_k(GlobalC::GridT);
    if (GlobalV::KS_SOLVER == "genelpa" || GlobalV::KS_SOLVER == "scalapack_gvx"
        || GlobalV::KS_SOLVER == "lapack") // Peize Lin test 2019-05-15
    {
        this->cal_dm(this->wg, psi, dm_k_2d);
    }

    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        ModuleBase::GlobalFunc::ZEROS(this->charge->rho[is], this->charge->nrxx); // mohan 2009-11-10
    }

    //------------------------------------------------------------
    // calculate the charge density on real space grid.
    //------------------------------------------------------------

    ModuleBase::GlobalFunc::NOTE("Calculate the charge on real space grid!");
    //uhm.GK.cal_rho_k(this->loc->DM_R);

    this->charge->renormalize_rho();

    ModuleBase::timer::tick("ElecStateLCAO", "psiToRho");
    return;
}

// Gamma_only case
void ElecStateLCAO::psiToRho(const psi::Psi<double>& psi)
{
    ModuleBase::TITLE("ElecStateLCAO", "psiToRho");
    ModuleBase::timer::tick("ElecStateLCAO", "psiToRho");

    if (GlobalV::KS_SOLVER == "genelpa" || GlobalV::KS_SOLVER == "scalapack_gvx"
                || GlobalV::KS_SOLVER == "lapack")
    {
        // LiuXh modify 2021-09-06, clear memory, cal_dk_gamma() not used for genelpa solver.
        // density matrix has already been calculated.
        ModuleBase::timer::tick("ElecStateLCAO", "cal_dm_2d");

        psi::Psi<double> dm_gamma_2d(psi.get_nk(), psi.get_nbasis(), psi.get_nbasis());
        // caution:wfc and dm
        this->cal_dm(this->wg, psi, dm_gamma_2d);

        ModuleBase::timer::tick("ElecStateLCAO", "cal_dm_2d");

        //this->loc->cal_dk_gamma_from_2D(); // transform dm_gamma[is].c to this->loc->DM[is]
    }

    for (int is = 0; is < GlobalV::NSPIN; is++)
    {
        ModuleBase::GlobalFunc::ZEROS(this->charge->rho[is], this->charge->nrxx); // mohan 2009-11-10
    }

    //------------------------------------------------------------
    // calculate the charge density on real space grid.
    //------------------------------------------------------------
    ModuleBase::GlobalFunc::NOTE("Calculate the charge on real space grid!");
    //uhm.GG.cal_rho(this->loc->DM);

    this->charge->renormalize_rho();

    ModuleBase::timer::tick("ElecStateLCAO", "psiToRho");
    return;
}

} // namespace elecstate