#ifndef LOCAL_ORBITAL_CHARGE
#define LOCAL_ORBITAL_CHARGE

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/matrix.h"
#include "../module_base/complexmatrix.h"
#include "../src_parallel/parallel_common.h"
#include "grid_technique.h"
#include "src_lcao/record_adj.h"
#include "src_lcao/local_orbital_wfc.h"
#include "src_lcao/LCAO_hamilt.h"

class Local_Orbital_Charge
{

	public:

	Local_Orbital_Charge();
	~Local_Orbital_Charge();

	// mohan added 2021-02-08
    void allocate_dm_wfc(const Grid_Technique& gt,
        Local_Orbital_wfc &lowf);
    // sum bands to compute the electron charge density
	void sum_bands(LCAO_Hamilt &UHM);

	//-----------------
	// in DM_gamma.cpp
	//-----------------
	void allocate_gamma(const Grid_Technique &gt);

    void gamma_file(const Grid_Technique& gt,
        Local_Orbital_wfc &lowf);
    void cal_dk_gamma_from_2D_pub(void);


	//-----------------
	// in DM_k.cpp
	//-----------------
	void allocate_DM_k(void);

	// liaochen modify on 2010-3-23 
	// change its state from private to public
	double*** DM;	
	double** DM_R;

	// whether to printout density matrix
    static int out_dm; // output density matrix or not.

	void write_dm(const int &is, const int &iter, const std::string &fn, const int &precision);

	void read_dm(const int &is, const std::string &fn);

    
    //-----------------
	// in dm_2d.cpp
    //-----------------
    // dm stands for density matrix
    std::vector<ModuleBase::matrix> dm_gamma;			// dm_gamma[is](iw1,iw2);
    std::vector<ModuleBase::ComplexMatrix> dm_k;		// dm_k[ik](iw1,iw2);

    void init_dm_2d(void);
    
    // dm = wfc.T * wg * wfc.conj(); used in gamma_only
    void cal_dm(const ModuleBase::matrix& wg,   // wg(ik,ib), cal all dm 
        std::vector<ModuleBase::matrix>& wfc_gamma,
        std::vector<ModuleBase::matrix>& dm_gamma);

    // in multi-k,  it is dm(k)
    void cal_dm(const ModuleBase::matrix& wg,    // wg(ik,ib), cal all dm 
        std::vector<ModuleBase::ComplexMatrix>& wfc_k,
        std::vector<ModuleBase::ComplexMatrix>& dm_k);

    // dm(R) = wfc.T * wg * wfc.conj()*kphase, only used in multi-k 
    void cal_dm_R(
        std::vector<ModuleBase::ComplexMatrix>& dm_k,
        Record_adj& ra,
        double** dm2d);     //output, dm2d[NSPIN][LNNR]

    //-----------------
	// wavefunctions' pointer
    //-----------------
    Local_Orbital_wfc* LOWF;
    //-----------------
	// Parallel Variables' pointer
    //-----------------
    const Parallel_Orbitals* ParaV;

private:

	// whether the DM array has been allocated
	bool init_DM;
	// whether the DM(R) array has been allocated
	bool init_DM_R;

	void cal_dk_gamma(void);

	void cal_dk_k(const Grid_Technique &gt);

	// mohan add 2010-09-06
	int lgd_last;// sub-FFT-mesh orbitals number in previous step.
	int lgd_now;// sub-FFT-mesh orbitals number in this step.

	int nnrg_last;// sub-FFT-mesh orbtials number in previous step, with k.
	int nnrg_now; // sub-FFT-mesh orbitals number in this step, with k.

	// add by yshen on 9/22/2014
	// these variables are memory pool for DM series matrixes, 
	// so that these matrixes will be storaged continuously in the memory.
	double **DM_pool;
	
	// Buffer parameters for tranforming 2D block-cyclic distributed DM matrix 
	// to grid distributed DM matrix
    int *sender_2D_index;
    int sender_size;
    int *sender_size_process;
    int *sender_displacement_process;
    double* sender_buffer;

    int *receiver_local_index;
    int receiver_size;
    int *receiver_size_process;
    int *receiver_displacement_process;
    double* receiver_buffer;

#ifdef __MPI
    int setAlltoallvParameter(MPI_Comm comm_2D, int blacs_ctxt, int nblk);
#endif
    void cal_dk_gamma_from_2D(void);
};

#endif
