#ifndef HSOLVERLCAO_H
#define HSOLVERLCAO_H

#include "hsolver.h"

namespace ModuleHSolver
{

class HSolverLCAO : public HSolver
{
    public:

    //HSolverLCAO(const PW_Basis* pbas_in){
    //    this->pbas = pbas_in;}
        /*this->init(pbas_in);*/

    /*void init(
        const Basis* pbas
        //const Input &in,
    ) override;
    void update(//Input &in
    ) override;*/
    
    void solve(
        ModuleHamilt::Hamilt* pHamilt, 
        ModulePsi::Psi<std::complex<double>>& psi, 
        ModuleElecS::ElecState* pes) override;
    

    private:
    void hamiltSolvePsiK(ModuleHamilt::Hamilt* hm, ModulePsi::Psi<std::complex<double>>& psi, double* eigenvalue);

    //const PW_Basis* pbas;
};

}

#endif