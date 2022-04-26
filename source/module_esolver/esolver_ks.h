#ifndef ESOLVER_KS_H
#define ESOLVER_KS_H
#include "./esolver_fp.h"
#include "string.h"
#include "fstream"
// #include "estates.h"
// #include "h2e.h"
namespace ModuleESolver
{

class ESolver_KS: public ESolver_FP
{
    public:
        ESolver_KS();
        // HSolver* phsol;
        double diag_ethr; //diag threshold
        double scf_thr; //scf threshold
        double drho;
        double niter;
        virtual void Run(int istep, UnitCell_pseudo& cell) override;
        virtual void hamilt2density(int istep, double ethr);

    protected:
        //Something to do before iter loop
        virtual void beforeiter(){}; 
        //Something to do before hamilt2density function in each iter loop.
        virtual void eachiterinit(){}; 
        //Something to do after hamilt2density function in each iter loop.
        virtual void eachiterfinish(bool){}; 
        //Something to do after the iter loop when scf is converged or comes to the max iter step.
        virtual void afteriter(bool){};
        //temporary----------------------------------------------------------
        virtual void updatepot(bool conv){};
        

    //TOOLS:
    protected:
        //Set ethr for hsolver
        void set_ethr(int istep, int iter);
        //print the headline on the screen:
        //ITER   ETOT(eV)       EDIFF(eV)      SCF_THR    "ITERTAG"    TIME(s) 
        void printhead();
        //write the headline in the running_log file
        //"basisname" ALGORITHM --------------- ION=   1  ELEC=   1--------------------------------
        void writehead(std::ofstream &ofs_running,int istep, int iter);



    protected:
        std::string basisname;
        std::string diagname;

};
}
#endif