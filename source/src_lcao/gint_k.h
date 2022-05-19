#ifndef GINT_K_H
#define GINT_K_H

#include "../module_orbital/ORB_atomic_lm.h"
#include "grid_technique.h"
#include "LCAO_matrix.h"
#include "../src_pw/charge.h"
#include "gint_tools.h"

// add by jingan for map<> in 2021-12-2, will be deleted in the future
#include "../src_ri/abfs-vector3_order.h"

class Gint_k
{
    public:

    Gint_k();
    ~Gint_k();

    // preparing FFT grid
    void prep_grid(
        const int &nbx_in,
		const int &nby_in,
		const int &nbz_in,
		const int &nbz_start_in,
        const int& ncxyz_in);

    // allocate the <phi_0 | V | phi_R> matrix element.
    void allocate_pvpR(void);
    void allocate_pvpR_tr(void); //LiuXh add 2019-07-15

    // destroy the temporary <phi_0 | V | phi_R> matrix element.
    void destroy_pvpR(void);
    //LiuXh add 2019-07-15
    void destroy_pvpR_tr(void);
    void distribute_pvpR_tr(LCAO_Matrix *LM);

    // jingan add 2021-6-4, modify 2021-12-2
    void distribute_pvpR_sparseMatrix(
        const int current_spin, 
        const double &sparse_threshold, 
        const std::map<Abfs::Vector3_Order<int>,
        std::map<size_t, std::map<size_t, double>>> &pvpR_sparseMatrix,
        LCAO_Matrix *LM
    );
    void distribute_pvpR_soc_sparseMatrix(
        const double &sparse_threshold, 
        const std::map<Abfs::Vector3_Order<int>,
        std::map<size_t,
        std::map<size_t, std::complex<double>>>> &pvpR_soc_sparseMatrix,
        LCAO_Matrix *LM
    );
    void cal_vlocal_R_sparseMatrix(
        const int &current_spin,
        const double &sparse_threshold,
        LCAO_Matrix *LM);

    // reset the spin.
    void reset_spin(const int &spin_now);

    // get the spin.
    int get_spin(void)const{return spin_now;}

    //------------------------------------------------------
    // in gint_k_vl.cpp 
    //------------------------------------------------------
    // calculate the matrix elements of Hamiltonian matrix,
    // < phi_0 | Vl + Vh + Vxc | phi_R> or if the Vna is used,
    // < phi_0 | delta_Vh + Vxc | phi_R>.
    void cal_vlocal_k(const double* vrs1, const Grid_Technique &gt, const int spin=0);

    // folding the < phi_0 | V | phi_R> matrix to 
    // <phi_0i | V | phi_0j>
    // V is (Vl + Vh + Vxc) if no Vna is used,
    // and is (Vna + delta_Vh + Vxc) if Vna is used.
    void folding_vl_k(const int &ik, LCAO_Matrix* LM);

    void folding_vl_k_nc(const int &ik, LCAO_Matrix* LM);//zhengdy-soc

    //------------------------------------------------------
    // in gint_k_rho.cpp 
    //------------------------------------------------------
    // calculate the charge density via grid integrals
    void cal_rho_k(double** DM_R, Charge* chr);

    //------------------------------------------------------
    // in gint_k_env.cpp 
    //------------------------------------------------------
    // calculate the envelop function via grid integrals
    void cal_env_k(
        int ik, 
        const std::complex<double>* wfc_k,
        double* rho);

    //------------------------------------------------------
    // in gint_k_fvl.cpp 
    //------------------------------------------------------
    // calculate force & stress (many k-points).

    void cal_force_k(
        const bool isforce,
        const bool isstress,
        ModuleBase::matrix& fvl_dphi, 
        ModuleBase::matrix& svl_dphi, 
        const double* vl,
        double **DM_R);
        //mohan add 2011-06-19 initial implementation
        //zhengdy add 2016-10-18 add stress calculation
        //wenfei modify 2022-5-17 reconstruction

    private:
    
    void cal_meshball_vlocal(
        int na_grid,
        int LD_pool,
        int grid_index, 
        int* block_size,
        int* block_index,
        int* block_iw,
        bool** cal_flag, 
        int* at, 
        double** psir_ylm,
        double** psir_vlbr3,
        double* pvpR);

    void cal_meshball_force(
        const int grid_index,
        const int na_grid,  					    // how many atoms on this (i,j,k) grid
        const int*const block_size, 			    // block_size[na_grid],	number of columns of a band
        const int*const block_index,		    	// block_index[na_grid+1], count total number of atomis orbitals
        const double*const*const psir_vlbr3_DMR,	    // psir_vlbr3[GlobalC::pw.bxyz][LD_pool]
        const double*const*const dpsir_x,	    // psir_vlbr3[GlobalC::pw.bxyz][LD_pool]
        const double*const*const dpsir_y,	    // psir_vlbr3[GlobalC::pw.bxyz][LD_pool]
        const double*const*const dpsir_z,	    // psir_vlbr3[GlobalC::pw.bxyz][LD_pool]
        ModuleBase::matrix &force);

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
        ModuleBase::matrix &stress);

    private:

    //----------------------------
    // key variable 
    //----------------------------

    // variables related to FFT grid
 	int nbx;
	int nby;
	int nbz;
	int ncxyz;
	int nbz_start;   

    // used only in vlocal.
    // dimension: [GlobalC::LNNR.nnrg] 
    // save the < phi_0i | V | phi_Rj > in sparse H matrix.
    double** pvpR_reduced;
    double***** pvpR_tr; //LiuXh add 2019-07-15
    std::complex<double>***** pvpR_tr_soc; //LiuXh add 2019-07-15
    
    bool pvpR_alloc_flag;
    
    int spin_now;

};

#endif
