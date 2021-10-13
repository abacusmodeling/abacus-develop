#ifndef LCAO_DESCRIPTOR_H
#define LCAO_DESCRIPTOR_H
#endif
#ifdef __DEEPKS

#include "../module_base/intarray.h"
#include "../module_base/complexmatrix.h"
#include <torch/script.h>
#include "../src_pw/global.h"

///
/// This class computes the descriptors for each atom from LCAO basis set,
/// interfaces with pytorch to obtain the correction potential in LCAO basis,
/// and computes the forces according to the correction potential.
/// 
/// For details of DeePKS method, you can refer to [DeePKS paper](https://pubs.acs.org/doi/10.1021/acs.jctc.0c00872).
///
//
// caoyu add 2021-03-29
//
class LCAO_Descriptor
{

//-------------------
// public functions
//-------------------
public:

    explicit LCAO_Descriptor();
    ~LCAO_Descriptor();

    ///only for descriptor part, not including scf
    void init(const int lm/**< [in] max angular momentum quantum number: 'L'*/,
        const int nm/**< [in] max orbital number with the same 'L', for each 'L'*/,
        const int tot_inl/**< [in] total number of radial orbitals (sum of atoms 'I', angular number 'L' and orbital number 'N') */);

	/// calculate\f$S_{\alpha, \mu} = \langle\alpha|\phi_\mu\rangle \f$ overlap between lcao basis Phi and descriptor basis Alpha
    void build_S_descriptor(const bool &calc_deri/**< [in] 0 for \f$\langle\phi|\alpha\rangle\f$, 1 for \f$\langle\frac{d\phi}{dR}|\alpha\rangle\f$*/);


	/// 1. Load DeePKS model
    /// 2. Initialize the deltaV Hamiltonian matrix 
    /// 3. If FORCE, initialize the matrces for force
    void deepks_pre_scf(const std::string& model_file/**< [in] path of a traced model file, provided by deepks-kit*/);


    //S_alpha_mu * DM  * S_nu_beta
    ///calculate projected density matrix:
    ///\f[D^ I_{ nlmm'} = \sum_{i}\sum_{\mu, \nu}\langle\alpha^I_{nlm}|\phi_\mu\rangle c_{i,\mu}c_{i,\nu} \langle\phi_\nu|\alpha^I_{nlm' }\rangle\f]
    void cal_projected_DM(const ModuleBase::matrix& dm/**< [in] density matrix*/);
    
    ///EIGENVALUE of pdm in block of I_n_l
    void cal_descriptor(void);
    
    ///compute the descriptor for each atom
    void cal_dm_as_descriptor(const ModuleBase::matrix& dm/**< [in] density matrix*/); // mohan add 2021-08-04

    ///calculate \f$\frac{dE_\delta}{dD^I_{nlmm'}}\f$
    void cal_gedm(const ModuleBase::matrix& dm/**< [in] density matrix*/);	//need to load model in this step

    ///calculate \f$\frac{d\mathbf{d^I_{nl}}}{dX}=\frac{d\mathbf{d^I_{nl}}}{dD^I_{nlmm'}}*\frac{dD^I_{nlmm'}}{dX}\f$
    ///----------------------------------------------------
    ///m, n: 2*l+1
    ///v: eigenvalues of dm , 2*l+1
    ///a,b: natom 
    /// - (a: the center of descriptor orbitals
    /// - b: the atoms whose force being calculated)
    ///gvdm*gdmx->gvx
    ///----------------------------------------------------
    void cal_gvx(const ModuleBase::matrix &dm);
    
    ///calculate \f$\sum_{I}\sum_{nlmm'}\langle\phi_\mu|\alpha^I_{nlm}\rangle{\frac{dE}{dD^I_{nlmm'}}}\langle\alpha^I_{nlm'}|\phi_\nu\rangle\f$ (for gamma_only)
    void build_v_delta_alpha(const bool& cal_deri/**< [in] 0 for 3-center intergration, 1 for its derivation*/);
    
    ///calculate \f$\sum_{I}\sum_{nlmm'}\langle\phi_\mu|\alpha^I_{nlm}\rangle{\frac{dE}{dD^I_{nlmm'}}}\langle\alpha^I_{nlm'}|\phi_\nu\rangle\f$ (for multi-k)
    void build_v_delta_mu(const bool &cal_deri/**< [in] 0 for 3-center intergration, 1 for its derivation*/);
    
    ///compute \f$H_{\delta, \mu\nu} = \langle\phi_\mu|V_\delta|\phi_\nu\rangle\f$ 
    void cal_v_delta(const ModuleBase::matrix& dm/**< [in] density matrix*/);
    
    ///add \f$H_{\delta, \mu\nu}\f$ to the Hamiltonian matrix
    void add_v_delta(void);

    ///compute Hellmann-Feynman term of the force contribution of \f$E_\delta\f$
    void cal_f_delta_hf(const ModuleBase::matrix& dm/**< [in] density matrix*/);
    
    ///compute Pulay  term of the force contribution of \f$E_\delta\f$
    void cal_f_delta_pulay(const ModuleBase::matrix& dm/**< [in] density matrix*/);
    
    ///compute the force contribution of \f$E_\delta\f$
    void cal_f_delta(const ModuleBase::matrix& dm/**< [in] density matrix*/);


    ///print descriptors based on LCAO basis
    void print_descriptor(void);
    
    ///print the \f$H_\delta\f$ matrix in LCAO basis
    void print_H_V_delta(void);
    
    ///print the force related to\f$V_\delta\f$ for each atom
    void print_F_delta(const std::string &fname/**< [in] the name of output file*/);

	///----------------------------------------------------------------------
	///The following 3 functions save the `[dm_eig], [e_base], [f_base]`
	///of current configuration as `.npy` file, when `deepks_scf = 1`.
	///After a full group of consfigurations are calculated,
    ///we need a python script to `load` and `torch.cat` these `.npy` files,
    ///and get `l_e_delta,npy` and `l_f_delta.npy` corresponding to the exact E, F data.
    /// 
    /// Unit of energy: Ry
    ///
    /// Unit of force: Ry/Bohr
    ///----------------------------------------------------------------------
	void save_npy_d(void);
	void save_npy_e(const double &e/**<[in] \f$E_{base}\f$ or \f$E_{tot}\f$, in Ry*/, const std::string &e_file);
	void save_npy_f(const ModuleBase::matrix &fbase/**<[in] \f$F_{base}\f$ or \f$F_{tot}\f$, in Ry/Bohr*/, const std::string &f_file);
    void save_npy_gvx(void);

    ///calculate \f$tr(\rho H_\delta), \rho = \sum_i{c_{i, \mu}c_{i,\nu}} \f$ (for gamma_only)
    void cal_e_delta_band(const std::vector<ModuleBase::matrix>& dm/**<[in] density matrix*/);
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

	//density matrix: dm_gamma
	double* dm_double;
	// overlap between lcao and descriptor basis
	double** S_mu_alpha;	//[tot_Inl][GlobalV::NLOCAL][2l+1]	caoyu modified 2021-05-07

	//d(S) for f_delta:	<\psi_mu|d\alpha^I_nlm> , [tot_Inl][GlobalV::NLOCAL][2l+1]
	double** DS_mu_alpha_x;
	double** DS_mu_alpha_y;
	double** DS_mu_alpha_z;

    double* DH_V_delta_x;
    double* DH_V_delta_y;
    double* DH_V_delta_z;

    // projected density matrix
	double** pdm;	//[tot_Inl][2l+1][2l+1]	caoyu modified 2021-05-07
	std::vector<torch::Tensor> pdm_tensor;

	// descriptors
    double *d;
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

    //d(d)/dD, autograd from torch::symeig
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

	void init_index(void);	// index of descriptor in all atoms

	void set_S_mu_alpha(
		const int &iw1_all,
		const int& inl,
		const int& im,
		const double& v);

    void print_projected_DM(
		std::ofstream &ofs,
		ModuleBase::ComplexMatrix &des,
		const int &it,
		const int &ia,
		const int &l,
		const int& n);

	void set_DS_mu_alpha(
		const int& iw1_all,
		const int& inl,
		const int& im,
		const double& vx,
		const double& vy,
		const double& vz);

	void init_gdmx(void);
    void load_model(const std::string& model_file);
    
    void cal_gvdm();    //called when force=1, precondition of cal_gvx
    void cal_gdmx(const ModuleBase::matrix& dm);	//dD/dX, precondition of cal_gvx
	void del_gdmx(void);

	void getdm_double(const ModuleBase::matrix& dm);

	void cal_descriptor_tensor(void);

};
namespace GlobalC
{
extern LCAO_Descriptor ld;
}

#endif
