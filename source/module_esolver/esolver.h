#ifndef ESOLVER_H
#define ESOLVER_H

#include "../input.h"
#include "../module_cell/unitcell_pseudo.h"
#include "../src_pw/energy.h"
#include "../module_base/matrix.h"
namespace ModuleESolver
{

class ESolver
{
// protected:
//     Atom atom;
//     ModuleBase::matrix lattice_v;
public:
    ESolver(){
        tag = "ESolver";
    }
    // virtual ~ESolver() = 0;

    //virtual void Init(Input_EnSolver &inp, matrix &lattice_v)=0
    virtual void Init(Input &inp, UnitCell_pseudo &cell)=0;

    // They shoud be add after atom class is refactored
    // virtual void UpdateLatAtom(ModuleBase::matrix &lat_in, Atom &atom_in);
    // virtual void UpdateLat(ModuleBase::matrix &lat_in);
    // virtual void UpdateAtom(Atom &atom_in);
   
    //virtual void Run(int istep, Atom &atom) = 0;
    virtual void Run(int istep, UnitCell_pseudo &cell) = 0;
    
    virtual void cal_Energy(energy &en) = 0; 
    virtual void cal_Force(ModuleBase::matrix &force) = 0;
    virtual void cal_Stress(ModuleBase::matrix &stress) = 0;
    
    //Print current classname.
    virtual void printag();
    virtual int getiter(){};
    string tag;
};

void init_esolver(ESolver* &p_ensolver, const string use_esol);
void clean_esolver(ESolver* &pesolver);

}

#endif