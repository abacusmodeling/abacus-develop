#ifndef EN_SOLVER_H
#define EN_SOLVER_H

#include "../input.h"
#include "../module_cell/unitcell_pseudo.h"
#include "../src_pw/energy.h"
#include "../module_base/matrix.h"
namespace ModuleEnSover
{

class En_Solver
{
public:
    En_Solver(){
        tag = "En_Solver";
    }
    // virtual ~En_Solver() = 0;

    //virtual void Init(Input_EnSolver &inp, matrix &lattice_v)=0
    virtual void Init(Input &inp, UnitCell_pseudo &cell)=0;
   
    //virtual void Run(int istep, Atom &atom) = 0;
    virtual void Run(int istep, UnitCell_pseudo &cell) = 0;
    
    virtual void cal_Energy(energy &en) = 0; 
    virtual void cal_Force(ModuleBase::matrix &force) = 0;
    virtual void cal_Stress(ModuleBase::matrix &stress) = 0;
    virtual void printag();
    string tag;
};

void init_esolver(En_Solver* &p_ensolver, const string use_esol);
void clean_esolver(En_Solver* &pesolver);

}

#endif