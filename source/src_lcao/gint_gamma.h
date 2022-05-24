//=========================================================
//AUTHOR : mohan
//DATE : 2009-09-16
//REFACTOR : Peize Lin, 2021.06.28
//=========================================================
#ifndef GINT_GAMMA_H
#define GINT_GAMMA_H
#include "gint_interface.h"
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

class Gint_Gamma : public Gint_Interface
{
	public:

	Gint_Gamma();
	~Gint_Gamma();

	// the unified interface to grid integration
	void cal_gint_gamma(Gint_inout *inout);
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
    // calculate the matrix elements of Hamiltonian matrix,	
	void gint_kernel_vlocal(
		const int na_grid,
		const int grid_index,
		const double delta_r,
		double* vldr3,
		const int LD_pool,
		double* pvpR_grid_in);

	void cal_meshball_vlocal(
		const int na_grid,  						// how many atoms on this (i,j,k) grid
		const int LD_pool,
		const int*const block_iw,					// block_iw[na_grid],	index of wave functions for each block
		const int*const block_size, 				// block_size[na_grid],	number of columns of a band
		const int*const block_index,				// block_index[na_grid+1], count total number of atomis orbitals
		const bool*const*const cal_flag,			// cal_flag[GlobalC::pw.bxyz][na_grid],	whether the atom-grid distance is larger than cutoff
		const double*const vldr3,					// vldr3[GlobalC::pw.bxyz]
		const double*const*const psir_ylm,			// psir_ylm[GlobalC::pw.bxyz][LD_pool]
		const double*const*const psir_vlbr3,		// psir_vlbr3[GlobalC::pw.bxyz][LD_pool]
		double* GridVlocal);		// GridVlocal[lgd_now][lgd_now]

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
