#include "hsolver_pw.h"

#include "diago_cg.h"
#include "diago_david.h"
#include "module_base/tool_quit.h"
#include "module_base/timer.h"
#include "module_elecstate/elecstate_pw.h"
#include "src_pw/global.h"

namespace hsolver
{

/*void HSolverPW::init(const PW_Basis* pbas_in)
{
    this->pbas = pbas_in;
    return;
}

void HSolverPW::update()
{
    return;
}*/

void HSolverPW::solve(hamilt::Hamilt* pHamilt, psi::Psi<std::complex<double>>& psi, elecstate::ElecState* pes, const std::string method_in, const bool skip_charge)
{
    ModuleBase::TITLE("HSolverPW", "solve");
    ModuleBase::timer::tick("HSolverPW", "solve");
    // prepare for the precondition of diagonalization
    this->precondition.resize(psi.get_nbasis());

    // select the method of diagonalization
    this->method = method_in;
    if (this->method == "cg")
    {
        if(pdiagh!=nullptr)
        {
            if(pdiagh->method != this->method)
            {
                delete[] pdiagh;
                pdiagh = new DiagoCG(precondition.data());
                pdiagh->method = this->method;
            }
        }
        else
        {
            pdiagh = new DiagoCG(precondition.data());
            pdiagh->method = this->method;
        }
    }
    else if (this->method == "dav")
    {
        DiagoDavid::PW_DIAG_NDIM = GlobalV::PW_DIAG_NDIM;
        if (pdiagh != nullptr)
        {
            if (pdiagh->method != this->method)
            {
                delete[] pdiagh;
                pdiagh = new DiagoDavid( precondition.data());
                pdiagh->method = this->method;
            }
        }
        else
        {
            pdiagh = new DiagoDavid( precondition.data());
            pdiagh->method = this->method;
        }
    }
    else
    {
        ModuleBase::WARNING_QUIT("HSolverPW::solve", "This method of DiagH is not supported!");
    }

    /// Loop over k points for solve Hamiltonian to charge density
    for (int ik = 0; ik < this->wfc_basis->nks; ++ik)
    {
        /// update H(k) for each k point
        pHamilt->updateHk(ik);

        this->updatePsiK(pHamilt, psi, ik);

        // template add precondition calculating here
        update_precondition(precondition, ik, this->wfc_basis->npwk[ik]);

        /// solve eigenvector and eigenvalue for H(k)
        double* p_eigenvalues = &(pes->ekb(ik, 0));
        this->hamiltSolvePsiK(pHamilt, psi, p_eigenvalues);
        /// calculate the contribution of Psi for charge density rho
    }

    // DiagoCG would keep 9*nbasis memory in cache during loop-k
    // it should be deleted before calculating charge
    if(this->method == "cg")
    {
        delete (DiagoCG*)pdiagh;
        pdiagh = nullptr;
    }
    //psi only should be initialed once for PW
    if(!this->initialed_psi)
    {
        this->initialed_psi = true;
    }

    if(skip_charge)
    {
        ModuleBase::timer::tick("HSolverPW", "solve");
        return;
    }
    pes->psiToRho(psi);

    ModuleBase::timer::tick("HSolverPW", "solve");
    return;
}

void HSolverPW::updatePsiK(hamilt::Hamilt* pHamilt, psi::Psi<std::complex<double>>& psi, const int ik)
{
    psi.fix_k(ik);
    if(!this->initialed_psi)
    {
        if(GlobalV::BASIS_TYPE=="pw")
        {
            // generate PAOs first, then diagonalize to get
            // inital wavefunctions.
            GlobalC::wf.diago_PAO_in_pw_k2(ik, psi, pHamilt);
        }
        else
        {
            ModuleBase::WARNING_QUIT("HSolverPW::updatePsiK", "lcao_in_pw is not supported now.");
        }
        return;
    }
}

void HSolverPW::hamiltSolvePsiK(hamilt::Hamilt* hm, psi::Psi<std::complex<double>>& psi, double* eigenvalue)
{
    pdiagh->diag(hm, psi, eigenvalue);
}

void HSolverPW::update_precondition(std::vector<double> &h_diag, const int ik, const int npw)
{
    h_diag.resize(h_diag.size(), 1.0);
    int precondition_type = 2;
    const double tpiba2 = this->wfc_basis->tpiba2;
    
    //===========================================
    // Conjugate-Gradient diagonalization
    // h_diag is the precondition matrix
    // h_diag(1:npw) = MAX( 1.0, g2kin(1:npw) );
    //===========================================
    if (precondition_type == 1)
    {
        for (int ig = 0; ig < npw; ig++)
        {
            double g2kin = this->wfc_basis->getgk2(ik,ig) * tpiba2;    
            h_diag[ig] = std::max(1.0, g2kin);
        }
    }
    else if (precondition_type == 2)
    {
        for (int ig = 0; ig < npw; ig++)
        {
            double g2kin = this->wfc_basis->getgk2(ik,ig) * tpiba2;
            h_diag[ig] = 1 + g2kin + sqrt(1 + (g2kin - 1) * (g2kin - 1));
        }
    }
    if(GlobalV::NSPIN==4)
    {
        const int size = h_diag.size();
        for (int ig = 0; ig < npw; ig++)
        {
            h_diag[ig+size/2] = h_diag[ig];
        }
    }
}

} // namespace hsolver