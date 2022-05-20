#ifndef GINT_K_H
#define GINT_K_H

#include "../module_orbital/ORB_atomic_lm.h"
#include "grid_technique.h"
#include "LCAO_matrix.h"
#include "../src_pw/charge.h"
#include "gint_tools.h"

// add by jingan for map<> in 2021-12-2, will be deleted in the future
#include "../src_ri/abfs-vector3_order.h"

namespace Gint_Tools
{
    enum class job_type{vlocal, rho, force};
}

class Gint_inout
{
    public:
    //input
        double** DM_R;
        double* vl;
        bool isforce;
        bool isstress;
    //output
        Charge* chr;
        ModuleBase::matrix* fvl_dphi;
        ModuleBase::matrix* svl_dphi;

        Gint_Tools::job_type job;

        void prep_gint_inout_rho(double **DM_R_in, Charge* chr_in, Gint_Tools::job_type job_in)
        {
            DM_R = DM_R_in;
            chr = chr_in;
            job = job_in;
        }

        void prep_gint_inout_force(double** DM_R_in, double* vl_in,
            bool isforce_in, bool isstress_in,
            ModuleBase::matrix* fvl_dphi_in,
            ModuleBase::matrix* svl_dphi_in,
            Gint_Tools::job_type job_in)
        {
            DM_R = DM_R_in;
            vl = vl_in;
            isforce = isforce_in;
            isstress = isstress_in;
            fvl_dphi = fvl_dphi_in;
            svl_dphi = svl_dphi_in;
            job = job_in;
        }
};

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

    //the unified interface
    void cal_gint_k(Gint_inout *inout);

 

    //------------------------------------------------------
    // in gint_k_vl.cpp 
    //------------------------------------------------------
    // calculate the matrix elements of Hamiltonian matrix,
    // < phi_0 | Vl + Vh + Vxc | phi_R> or if the Vna is used,
    // < phi_0 | delta_Vh + Vxc | phi_R>.
    void cal_vlocal_k(const double* vrs1, const Grid_Technique &gt, const int spin=0);


    //------------------------------------------------------
    // in gint_k_pvpr.cpp 
    //------------------------------------------------------
    // pvpR and reset_spin/get_spin : auxilliary methods
    // for calculating hamiltonian

    // reset the spin.
    void reset_spin(const int &spin_now_in){this->spin_now = spin_now_in;};
    // get the spin.
    int get_spin(void)const{return spin_now;}
 
    // allocate the <phi_0 | V | phi_R> matrix element.
    void allocate_pvpR(void);
    // destroy the temporary <phi_0 | V | phi_R> matrix element.
    void destroy_pvpR(void);

    // folding the < phi_0 | V | phi_R> matrix to 
    // <phi_0i | V | phi_0j>
    // V is (Vl + Vh + Vxc) if no Vna is used,
    // and is (Vna + delta_Vh + Vxc) if Vna is used.
    void folding_vl_k(const int &ik, LCAO_Matrix* LM);

    //------------------------------------------------------
    // in gint_k_env.cpp 
    //------------------------------------------------------
    // calculate the envelop function via grid integrals
    void cal_env_k(
        int ik, 
        const std::complex<double>* wfc_k,
        double* rho);

    //related to sparse matrix
    // jingan add 2021-6-4, modify 2021-12-2
    void distribute_pvpR_sparseMatrix(
        const int current_spin, 
        const double &sparse_threshold, 
        const std::map<Abfs::Vector3_Order<int>,
        std::map<size_t, std::map<size_t, double>>> &pvpR_sparseMatrix,
        LCAO_Matrix *LM);

    void distribute_pvpR_soc_sparseMatrix(
        const double &sparse_threshold, 
        const std::map<Abfs::Vector3_Order<int>,
        std::map<size_t,
        std::map<size_t, std::complex<double>>>> &pvpR_soc_sparseMatrix,
        LCAO_Matrix *LM);

    void cal_vlocal_R_sparseMatrix(
        const int &current_spin,
        const double &sparse_threshold,
        LCAO_Matrix *LM);

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
        Gint_inout *inout);

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
    bool pvpR_alloc_flag;
    int spin_now;

    double** pvpR_reduced;
    
};

#endif
