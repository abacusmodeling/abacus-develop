#ifndef HSOLVER_H
#define HSOLVER_H

#include "diagh.h"
#include "module_hamilt/hamilt.h"
#include "module_elecstate/elecstate.h"
#include "module_psi/psi.h"


namespace ModuleHSolver
{
    
class HSolver
{
    public:
    //initialization, used in construct function or restruct a new HSolver
    virtual void init( 
        /*const Basis* pbas, //We need Basis class here, use global class for this initialization first */
        /*const Input &in, //We need new Input class here, use global variable for this initialization first */
        ElecState *pes
        )=0;
    //initialization, only be called for change some parameters only
    virtual void init(
        /*Input &in */)=0;
        
    //solve Hamiltonian to electronic density in ElecState
    virtual void solve(ModuleHamilt::Hamilt* phm, 
                ModulePsi::Psi* ppsi,
                ModuleElecS::ElecState* pes) =0;
    
    protected:
    
    DiagH* pdiagh;          // for single Hamiltonian matrix diagonal solver

    

};

}
#endif