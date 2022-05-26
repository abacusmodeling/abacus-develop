#include "hsolver_lcao.h"

#include "diago_blas.h"
#include "diago_elpa.h"

namespace hsolver
{

template <typename T>
void HSolverLCAO::solveTemplate(hamilt::Hamilt* pHamilt, psi::Psi<T>& psi, elecstate::ElecState* pes, const std::string method_in, const bool skip_charge)
{
    // select the method of diagonalization
    this->method = method_in;
    if (this->method == "genelpa")
    {
        if (pdiagh != nullptr)
        {
            if (pdiagh->method != this->method)
            {
                delete[] pdiagh;
                pdiagh = new DiagoElpa();
                pdiagh->method = this->method;
            }
        }
        else
        {
            pdiagh = new DiagoElpa();
            pdiagh->method = this->method;
        }
    }
    else if (this->method == "scalapack_gvx")
    {
        if (pdiagh != nullptr)
        {
            if (pdiagh->method != this->method)
            {
                delete[] pdiagh;
                pdiagh = new DiagoBlas();
                pdiagh->method = this->method;
            }
        }
        else
        {
            pdiagh = new DiagoBlas();
            pdiagh->method = this->method;
        }
    }
    else
    {
        ModuleBase::WARNING_QUIT("HSolverLCAO::solve", "This method of DiagH is not supported!");
    }

    /// Loop over k points for solve Hamiltonian to charge density
    for (int ik = 0; ik < psi.get_nk(); ++ik)
    {
        /// update H(k) for each k point
        pHamilt->updateHk(ik);

        psi.fix_k(ik);

        /// solve eigenvector and eigenvalue for H(k)
        double* p_eigenvalues = &(pes->ekb(ik, 0));
        this->hamiltSolvePsiK(pHamilt, psi, p_eigenvalues);
    }

    if (this->method != "genelpa" && this->method != "scalapack_gvx")
    {
        delete pdiagh;
        pdiagh = nullptr;
    }

    //used in nscf calculation
    if(skip_charge) return;

    //calculate charge by psi
    //called in scf calculation
    pes->psiToRho(psi);
}

int HSolverLCAO::out_mat_hs = 0;
int HSolverLCAO::out_mat_hsR = 0;

void HSolverLCAO::solve(hamilt::Hamilt* pHamilt, psi::Psi<std::complex<double>>& psi, elecstate::ElecState* pes, const std::string method_in, const bool skip_charge)
{
    this->solveTemplate(pHamilt, psi, pes, method, skip_charge);
}
void HSolverLCAO::solve(hamilt::Hamilt* pHamilt, psi::Psi<double>& psi, elecstate::ElecState* pes, const std::string method_in, const bool skip_charge)
{
    this->solveTemplate(pHamilt, psi, pes, method, skip_charge);
}

void HSolverLCAO::hamiltSolvePsiK(hamilt::Hamilt* hm, psi::Psi<std::complex<double>>& psi, double* eigenvalue)
{
    pdiagh->diag(hm, psi, eigenvalue);
}

void HSolverLCAO::hamiltSolvePsiK(hamilt::Hamilt* hm, psi::Psi<double>& psi, double* eigenvalue)
{
    pdiagh->diag(hm, psi, eigenvalue);
}

} // namespace hsolver