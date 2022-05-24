//=========================================================
//AUTHOR : mohan
//DATE : 2009-09-16
//REFACTOR : Peize Lin, 2021.06.28
//=========================================================
#ifndef GINT_GAMMA_H
#define GINT_GAMMA_H
#include "gint.h"
#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "grid_technique.h"
#include "LCAO_matrix.h"
#include <omp.h>

//=========================================================
// ModuleBase::Integral On 3D Grids, different from Grid_Integral
// Feature : Matrix Elements Of Local Potential For 
// Numerical Orbitals
//=========================================================

class Gint_Gamma : public Gint
{
	public:

	Gint_Gamma();
	~Gint_Gamma();

	void cal_vlocal(Gint_inout *inout);

	// (4) calcualte the envelope function
	void cal_env(const double* wfc, double* rho);

	// (5) calculate the Mulliken charge
	void cal_mulliken(double** mulliken);

private:

    double***  DM;   //pointer to LOC.DM

	// for calculation of Mulliken charge.
	void gamma_mulliken(double** mulliken);
	// for calculation of envelope functions.
	void gamma_envelope(const double* wfc, double* rho);// mohan add 2011-07-01

    //------------------------------------------------------
    // in gint_gamma_vl.cpp 
    //------------------------------------------------------

    void vl_grid_to_2D(const int lgd, LCAO_Matrix& lm); //redistribute the Hamiltonian to 2D block format

    ///===============================
    /// Use MPI_Alltoallv to convert a grid distributed matrix
    /// to 2D - block cyclic distributed matrix.
    ///===============================
    int sender_index_size;
    int *sender_local_index;
    int sender_size;
    int *sender_size_process;
    int *sender_displacement_process;
    double* sender_buffer;

    int receiver_index_size;
    int *receiver_global_index;
    int receiver_size;
    int *receiver_size_process;
    int *receiver_displacement_process;
    double* receiver_buffer;

};

#endif
