#ifndef LCAO_DESCRIPTOR_H
#define LCAO_DESCRIPTOR_H
#endif
#ifdef __DEEPKS

#include "../module_base/intarray.h"
#include "../module_base/complexmatrix.h"
#include <torch/script.h>
#include "../src_pw/global.h"

//------------------------------------------------------------------------------
// This class computes the descriptors for each atom from LCAO basis set,
// interfaces with pytorch to obtain the correction potential in LCAO basis,
// and computes the forces according to the correction potential
//------------------------------------------------------------------------------
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

    // index of descriptor in all atoms
    //only for descriptor part, not including scf
    void init(const int lm, const int nm, const int tot_inl);

	// cal S_alpha_mu: overlap between lcao basis Phi and descriptor basis Al
    void build_S_descriptor(const bool &calc_deri);


	// 1. Load DeePKS model
    // 2. Initialize the deltaV Hamiltonian matrix
    // 3. If FORCE, initialize the matrces for force
    void deepks_pre_scf(const std::string& model_file);


	//------------------------------------------------------------------------------
	// cal_projected_DM: pdm = S_alpha_mu * inv(Sloc) * DM * inv(Sloc) * S_nu_beta
	// cal_descriptor: EIGENVALUE of pdm in block of I_n_l
	// cal_dm_as_descriptor: compute the descriptor for each atom
	// cal_v_delta: compute <psi|deltaV|psi>
	// add_v_delta: add <psi|deltaV|psi> to the Hamiltonian matrix
	// cal_f_delta: compute the force related to deltaV, input dm is density matrix
	//------------------------------------------------------------------------------
    void cal_projected_DM(const ModuleBase::matrix &dm);
    void cal_descriptor(void);
    void cal_dm_as_descriptor(const ModuleBase::matrix& dm); // mohan add 2021-08-04

    void cal_gedm(const ModuleBase::matrix& dm);	//need to load model in this step
    void build_v_delta_alpha(const bool& cal_deri);
    void build_v_delta_mu(const bool& cal_deri);
    void cal_v_delta(const ModuleBase::matrix& dm);
    void add_v_delta(void);

    void cal_f_delta_hf(const ModuleBase::matrix& dm);
    void cal_f_delta_pulay(const ModuleBase::matrix& dm);
    void cal_f_delta(const ModuleBase::matrix& dm);


	//----------------------------------------------------------------------
	// print_descriptors: print descriptors based on LCAO basis
	// print_H_V_delta: print the deltaV matrix in LCAO basis
	// print_F_delta: print the force related to deltaV for each atom
	//----------------------------------------------------------------------
	void print_descriptor(void);
	void print_H_V_delta(void);
	void print_F_delta(const std::string& fname);

	//----------------------------------------------------------------------
	/*These 3 functions save the [dm_eig], [e_base], [f_base]
	of current configuration as .npy file, when deepks_scf = 1.
	After a full group of consfigurations are calculated,
    we need a python script to 'load' and 'torch.cat' these .npy files,
    and get l_e_delta and l_f_delta corresponding to the exact e,f data.*/
	//----------------------------------------------------------------------
	void save_npy_d(void);
	void save_npy_e(const double &ebase);	//Ry
	void save_npy_f(const ModuleBase::matrix &fbase);//Ry

    void cal_e_delta_band(const std::vector<ModuleBase::matrix>& dm);	//tr[rho*H_V_delta]
//-------------------
// public variables
//-------------------
public:

	//------------------------------------------------------
    //E_delta: in Ry, correction energy provided by NN
    //e_delta_band: tr(dm*H_V_delta)
    //H_V_delta: correction term to the Hamiltonian matrix
	//F_delta: in Ry/Bohr, force due to the correction term
	//------------------------------------------------------
    double E_delta = 0.0;
    double e_delta_band=0.0;
    double* H_V_delta;
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

	//gdmx: dD/dX		\sum_{mu,nu} 4*c_mu*c_nu * <dpsi_mu/dx|alpha_m><alpha_m'|psi_nu>
	double*** gdmx;	//[natom][tot_Inl][2l+1][2l+1]
	double*** gdmy;
	double*** gdmz;

	//dE/dD, autograd from loaded model(E: Ry)
	double** gedm;	//[tot_Inl][2l+1][2l+1]

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

    void cal_gdmx(const ModuleBase::matrix& dm);	//dD/dX
	void del_gdmx(void);

	void getdm_double(const ModuleBase::matrix& dm);

	void cal_descriptor_tensor(void);

};
namespace GlobalC
{
extern LCAO_Descriptor ld;
}

#endif
