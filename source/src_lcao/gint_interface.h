#ifndef GINT_INTERFACE
#define GINT_INTERFACE

#include "gint_tools.h"

class Gint_Interface
{
    public:

    // the unified interface to grid integration
	void cal_gint(Gint_inout *inout);

    // preparing FFT grid
    void prep_grid(
        const int &nbx_in,
        const int &nby_in,
        const int &nbz_in,
        const int &nbz_start_in,
        const int& ncxyz_in);

    protected:
    // variables related to FFT grid
 	int nbx;
	int nby;
	int nbz;
	int ncxyz;
	int nbz_start;

    //------------------------------------------------------
    // in gint_k_vl.cpp 
    //------------------------------------------------------
    // calculate the matrix elements of Hamiltonian matrix,
    // < phi_0 | Vl + Vh + Vxc | phi_R> or if the Vna is used,
    // < phi_0 | delta_Vh + Vxc | phi_R>.
    void gint_kernel_vlocal(
        const int na_grid,
        const int grid_index,
        const double delta_r,
        double* vldr3,
        const int LD_pool,
        double* pvpR_reduced);

    void cal_meshball_vlocal(
        int na_grid,
        int LD_pool,
        int grid_index, 
        int* block_size,
        int* block_index,
        int* block_iw,
        bool** cal_flag, 
        double** psir_ylm,
        double** psir_vlbr3,
        double* pvpR);

    //------------------------------------------------------
    // in gint_k_fvl.cpp 
    //------------------------------------------------------
    // calculate vl contributuion to force & stress via grid integrals
    void gint_kernel_force(
        const int na_grid,
        const int grid_index,
        const double delta_r,
        double* vldr3,
        const int LD_pool,
        double**DM_R,
        const bool isforce,
        const bool isstress,
        ModuleBase::matrix* fvl_dphi,
        ModuleBase::matrix* svl_dphi);

    void cal_meshball_force(
        const int grid_index,
        const int na_grid,  					    // how many atoms on this (i,j,k) grid
        const int*const block_size, 			    // block_size[na_grid],	number of columns of a band
        const int*const block_index,		    	// block_index[na_grid+1], count total number of atomis orbitals
        const double*const*const psir_vlbr3_DMR,	    // psir_vlbr3[GlobalC::pw.bxyz][LD_pool]
        const double*const*const dpsir_x,	    // psir_vlbr3[GlobalC::pw.bxyz][LD_pool]
        const double*const*const dpsir_y,	    // psir_vlbr3[GlobalC::pw.bxyz][LD_pool]
        const double*const*const dpsir_z,	    // psir_vlbr3[GlobalC::pw.bxyz][LD_pool]
        ModuleBase::matrix *force);

    void cal_meshball_stress(
        const int na_grid,  					    // how many atoms on this (i,j,k) grid
        const int*const block_index,		    	// block_index[na_grid+1], count total number of atomis orbitals
        const double*const*const psir_vlbr3_DMR,
        const double*const*const dpsir_xx,
        const double*const*const dpsir_xy,
        const double*const*const dpsir_xz,
        const double*const*const dpsir_yy,
        const double*const*const dpsir_yz,
        const double*const*const dpsir_zz,
        ModuleBase::matrix *stress);

    //------------------------------------------------------
    // in gint_k_rho.cpp 
    //------------------------------------------------------
    // calculate the charge density via grid integrals
    void gint_kernel_rho(
        const int na_grid,
        const int grid_index,
        const double delta_r,
        int* vindex,
        const int LD_pool,
        Gint_inout *inout);

    void cal_meshball_rho(
        const int na_grid,
        int* block_index,
        int* vindex,
        double** psir_ylm,
        double** psir_DMR,
        double* rho);

    // dimension: [GlobalC::LNNR.nnrg] 
    // save the < phi_0i | V | phi_Rj > in sparse H matrix.
    bool pvpR_alloc_flag = false;
    double** pvpR_reduced;
    
	double* pvpR_grid; //stores Hamiltonian in grid format
};

#endif