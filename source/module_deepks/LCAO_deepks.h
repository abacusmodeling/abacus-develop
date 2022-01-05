#ifndef LCAO_DESCRIPTOR_H
#define LCAO_DESCRIPTOR_H

#ifdef __DEEPKS

#include "../module_base/intarray.h"
#include "../module_base/complexmatrix.h"
#include "../src_pw/global.h"
#include <unordered_map>

#include "torch/script.h"
#include "torch/csrc/autograd/autograd.h"
#include "torch/csrc/api/include/torch/linalg.h"

#include "../src_lcao/LCAO_matrix.h"
#include "../module_base/intarray.h"
#include "../module_base/complexmatrix.h"
#include "../src_lcao/global_fp.h"
#include "../src_pw/global.h"
#include "../src_io/winput.h"

///
/// The LCAO_Deepks contains subroutines for implementation of the DeePKS method in atomic basis.
/// In essential, it is a machine-learned correction term to the XC potential
/// in the form of delta_V=|alpha> V(D) <alpha|, where D is a list of descriptors
/// The subroutines may be roughly grouped into 3 groups
/// 1. generation of projected density matrices pdm=sum_i,occ <phi_i|alpha><alpha|phi_i>
///    and then descriptors D=eig(pdm)
///    as well as their gradients with regard to atomic position, dD/dX
/// 2. loading the model, which requires interfaces with pytorch
/// 3. apply the correction potential, delta_V, in Kohn-Sham Hamiltonian and calculation of force, stress
/// 
/// For details of DeePKS method, you can refer to [DeePKS paper](https://pubs.acs.org/doi/10.1021/acs.jctc.0c00872).
///
///
// caoyu add 2021-03-29
// wenfei modified 2022-1-5
//
class LCAO_Deepks
{

//===============================
//DeePKS Part 1
//deals with generation of descriptors as well as labels
//===============================
public:

//===============================
//DeePKS Part 1
//deals with generation of descriptors as well as labels
//realized in LCAO_descriptor.cpp
//===============================

    explicit LCAO_Deepks();
    ~LCAO_Deepks();

    ///Allocate memory and calculate the index of descriptor in all atoms. 
    ///(only for descriptor part, not including scf)
    void init(const LCAO_Orbitals &orb,
        const int nat,
        const int ntype,
        std::vector<int> na);

    //S_alpha_mu * DM  * S_nu_beta
    ///calculate projected density matrix:
    ///pdm = sum_i,occ <phi_i|alpha1><alpha2|phi_k>
    void cal_projected_DM(const ModuleBase::matrix& dm/**< [in] density matrix*/,
        const UnitCell_pseudo &ucell,
        const LCAO_Orbitals &orb,
        Grid_Driver &GridD,
        const Parallel_Orbitals &ParaO);
    void cal_projected_DM_k(const std::vector<ModuleBase::ComplexMatrix>& dm,
        const UnitCell_pseudo &ucell,
        const LCAO_Orbitals &orb,
        Grid_Driver &GridD,
        const Parallel_Orbitals &ParaO,
        const K_Vectors &kv);

    ///Calculates descriptors
    ///which are eigenvalues of pdm in blocks of I_n_l
	void cal_descriptor(void);

//===============================
//DeePKS Part 2
//deals with application of correction dV to Hamiltonian and force
//realized in LCAO_descriptor_dV.cpp
//===============================

    //load the trained neural network model
    void load_model(const std::string& model_file);

    /// 1. Initialize the deltaV Hamiltonian matrix 
    /// 2. If FORCE, initialize the matrces for force
    void allocate_V_delta(const int nat,const int nloc, const int nks=1);
    void allocate_V_deltaR(const int nnr);

    ///add dV to the Hamiltonian matrix
    void add_v_delta(const UnitCell_pseudo &ucell,
        const LCAO_Orbitals &orb,
        Grid_Driver &GridD,
        const Parallel_Orbitals &ParaO);
    void add_v_delta_k(const UnitCell_pseudo &ucell,
        const LCAO_Orbitals &orb,
        Grid_Driver &GridD,
        const Parallel_Orbitals &ParaO,
        const int nnr_in);
    
    //calculates sum_(L0,M0) alpha<psi_i|alpha><alpha|psi_j>
    void build_v_delta_alpha(const bool& cal_deri/**< [in] 0 for 3-center intergration, 1 for its derivation*/,
        const UnitCell_pseudo &ucell,
        const LCAO_Orbitals &orb,
        Grid_Driver &GridD,
        const Parallel_Orbitals &ParaO,
        const ORB_gen_tables &UOT);

    //for gamma only, pulay and HF terms of force are calculated together
    void cal_f_delta_gamma(const ModuleBase::matrix& dm/**< [in] density matrix*/,
        const UnitCell_pseudo &ucell,
        const LCAO_Orbitals &orb,
        Grid_Driver &GridD,
        const Parallel_Orbitals &ParaO,
        const bool isstress, ModuleBase::matrix& svnl_dalpha);

    //for multi-k, pulay and HF terms of force are calculated together
    void cal_f_delta_k(const std::vector<ModuleBase::ComplexMatrix>& dm/**<[in] density matrix*/,
        const UnitCell_pseudo &ucell,
        const LCAO_Orbitals &orb,
        Grid_Driver &GridD,
        const Parallel_Orbitals &ParaO,
        const K_Vectors &kv,
        const bool isstress, ModuleBase::matrix& svnl_dalpha);

    ///calculate tr(\rho V_delta)
    void cal_e_delta_band(const std::vector<ModuleBase::matrix>& dm/**<[in] density matrix*/,
        const Parallel_Orbitals &ParaO);
    void cal_e_delta_band_k(const std::vector<ModuleBase::ComplexMatrix>& dm/**<[in] density matrix*/,
        const Parallel_Orbitals &ParaO,
        const int nks);

//============================
//DeePKS Part 3
//subroutines that deals with io as well as interface with libtorch
//realized in LCAO_descriptor_io.cpp
//============================

    ///calculate partial of energy correction to descriptors
    void cal_gedm(const int nat);

    ///calculates gradient of descriptors w.r.t atomic positions
    ///----------------------------------------------------
    ///m, n: 2*l+1
    ///v: eigenvalues of dm , 2*l+1
    ///a,b: natom 
    /// - (a: the center of descriptor orbitals
    /// - b: the atoms whose force being calculated)
    ///gvdm*gdmx->gvx
    ///----------------------------------------------------
    void cal_gvx(const int nat);

//calculate the gradient of pdm with regard to atomic positions
//d/dX D_{Inl,mm'}

    void cal_gdmx(const ModuleBase::matrix& dm,
        const UnitCell_pseudo &ucell,
        const LCAO_Orbitals &orb,
        Grid_Driver &GridD,
        const Parallel_Orbitals &ParaO);	//dD/dX, precondition of cal_gvx
    void cal_gdmx_k(const std::vector<ModuleBase::ComplexMatrix>& dm,
        const UnitCell_pseudo &ucell,
        const LCAO_Orbitals &orb,
        Grid_Driver &GridD,
        const Parallel_Orbitals &ParaO,
        const K_Vectors &kv);	//dD/dX, precondition of cal_gvx

    ///print descriptors based on LCAO basis
    void print_descriptor(const int nat);
    
    ///print the force related to\f$V_\delta\f$ for each atom
    void print_F_delta(const std::string &fname/**< [in] the name of output file*/, const UnitCell_pseudo &ucell);

	///----------------------------------------------------------------------
	///The following 4 functions save the `[dm_eig], [e_base], [f_base], [grad_vx]`
	///of current configuration as `.npy` file, when `deepks_scf = 1`.
	///After a full group of consfigurations are calculated,
    ///we need a python script to `load` and `torch.cat` these `.npy` files,
    ///and get `l_e_delta,npy` and `l_f_delta.npy` corresponding to the exact E, F data.
    /// 
    /// Unit of energy: Ry
    ///
    /// Unit of force: Ry/Bohr
    ///----------------------------------------------------------------------
	void save_npy_d(const int nat);
	void save_npy_e(const double &e/**<[in] \f$E_{base}\f$ or \f$E_{tot}\f$, in Ry*/,const std::string &e_file);
	void save_npy_f(const ModuleBase::matrix &fbase/**<[in] \f$F_{base}\f$ or \f$F_{tot}\f$, in Ry/Bohr*/,
        const std::string &f_file, const int nat);
    void save_npy_gvx(const int nat);

//========
// Extra : parallelization, also in LCAO_descriptor.cpp
//========

#ifdef __MPI
    //reduces a dim 2 array
    void allsum_deepks(
        int inlmax, //first dimension
        int ndim, //second dimension
        double** mat); //the array being reduced 
#endif

//-------------------
// public variables
//-------------------
public:
    
    ///(Unit: Ry) Correction energy provided by NN
    double E_delta = 0.0;
    ///(Unit: Ry)  \f$tr(\rho H_\delta), \rho = \sum_i{c_{i, \mu}c_{i,\nu}} \f$ (for gamma_only)
    double e_delta_band = 0.0;
    ///Correction term to the Hamiltonian matrix: \f$\langle\psi|V_\delta|\psi\rangle\f$
    double* H_V_delta;

    double* H_V_deltaR;
    std::complex<double>** H_V_delta_k;
    ///(Unit: Ry/Bohr) Total Force due to the DeePKS correction term \f$E_{\delta}\f$
    ModuleBase::matrix	F_delta;

//-------------------
// private variables
//-------------------
private:

	int lmaxd = 0;
	int nmaxd = 0;
	int inlmax = 0;

	// deep neural network module that provides corrected Hamiltonian term and
	// related derivatives.
	torch::jit::script::Module module;

    // saves <psi(0)|alpha(R)>
    std::vector<std::vector<std::unordered_map<int,std::vector<std::vector<double>>>>> nlm_save;
    
    typedef std::tuple<int,int,int,int> key_tuple;
    std::vector<std::map<key_tuple,std::unordered_map<int,std::vector<std::vector<double>>>>> nlm_save_k;

    // projected density matrix
	double** pdm;	//[tot_Inl][2l+1][2l+1]	caoyu modified 2021-05-07
	std::vector<torch::Tensor> pdm_tensor;

	// descriptors
	std::vector<torch::Tensor> d_tensor;

	//gedm:dE/dD, [tot_Inl][2l+1][2l+1]	(E: Hartree)
	std::vector<torch::Tensor> gedm_tensor;

	//gdmx: dD/dX		\sum_{mu,nu} 2*c_mu*c_nu * <dpsi_mu/dx|alpha_m><alpha_m'|psi_nu>
	double*** gdmx;	//[natom][tot_Inl][2l+1][2l+1]
	double*** gdmy;
	double*** gdmz;

	///dE/dD, autograd from loaded model(E: Ry)
	double** gedm;	//[tot_Inl][2l+1][2l+1]

    //gvx:d(d)/dX, [natom][3][natom][des_per_atom]
    torch::Tensor gvx_tensor;

    //d(d)/dD, autograd from torch::linalg::eigh
    std::vector<torch::Tensor> gevdm_vector;

    //dD/dX, tensor form of gdmx
    std::vector<torch::Tensor> gdmr_vector;

    ///size of descriptor(projector) basis set
    int n_descriptor;

	// \sum_L{Nchi(L)*(2L+1)}
	int des_per_atom;

	ModuleBase::IntArray* alpha_index;
	ModuleBase::IntArray* inl_index;	//caoyu add 2021-05-07
	int* inl_l;	//inl_l[inl_index] = l of descriptor with inl_index

//-------------------
// private functions
//-------------------
private:

//===============================
//DeePKS Part 1
//deals with generation of descriptors as well as labels
//realized in LCAO_descriptor.cpp
//===============================

// arrange index of descriptor in all atoms
	void init_index(const int ntype,
        const int nat,
        std::vector<int> na,
        const int tot_inl,
        const LCAO_Orbitals &orb);
// data structure that saves <psi|alpha>
    void allocate_nlm(const int nat);
// array for storing gdmx, used for calculating gvx
	void init_gdmx(const int nat);

	void del_gdmx(const int nat);

//============================
//DeePKS Part 3
//subroutines that deals with io as well as interface with libtorch
//realized in LCAO_descriptor_io.cpp
//============================    
    //partial of descriptors w.r.t projected dm
    //called when force=1, precondition of cal_gvx
    void cal_gvdm(const int nat);    

};

namespace GlobalC
{
    extern LCAO_Deepks ld;
}

#endif
#endif
