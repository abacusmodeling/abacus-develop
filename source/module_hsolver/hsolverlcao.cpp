#include "hsolverlcao.h"

#include "src_pw/global.h"

namespace ModuleHSolver
{

void HSolverLCAO::solve(
    ModuleHamilt::Hamilt* pHamilt, 
    ModulePsi::Psi<std::complex<double>>& psi, 
    ModuleElecS::ElecState* pes
)
{
    // select the method of diagonalization
    /*if(this->method == "genelpa") pdiagh = new DiagoELPA(&(GlobalC::hm.hpw), pbas, precondition.data());
    else if(this->method == "scalapack_gvx") pdiagh = new DiagoSca(&GlobalC::hm.hpw, pbas, precondition.data());
    else ModuleBase::WARNING_QUIT("HSolverPW::solve", "This method of DiagH is not supported!");*/

    ///Loop over k points for solve Hamiltonian to charge density 
    for(int ik=0;ik<psi.get_nk();++ik)
    {
        ///update H(k) for each k point
        pHamilt->updateHk(ik);

        psi.fix_k(ik);

        ///solve eigenvector and eigenvalue for H(k)
        double *p_eigenvalues = GlobalC::wf.ekb[ik];
        this->hamiltSolvePsiK(pHamilt, psi, p_eigenvalues);
        ///calculate the contribution of Psi for charge density rho
        pes->updateRhoK(psi);
    }
}

void HSolverLCAO::hamiltSolvePsiK(ModuleHamilt::Hamilt* hm, ModulePsi::Psi<std::complex<double>>& psi, double* eigenvalue)
{
    pdiagh->diag(hm, psi, eigenvalue);
}

}