#ifndef LOCAL_ORBITAL_WFC
#define LOCAL_ORBITAL_WFC

#include "grid_technique.h"
#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_orbital/ORB_control.h" // mohan add 2021-05-24

class Local_Orbital_wfc
{
	public:

	Local_Orbital_wfc();
	~Local_Orbital_wfc();

	// used to generate density matrix: LOC.DM_R,
	// which is used to calculate the charge density. 
	// which is got after the diagonalization of 
	// std::complex Hamiltonian matrix.
	std::complex<double>*** WFC_K; // [NK, GlobalV::NBANDS, GlobalV::NLOCAL]	
	std::complex<double>* WFC_K_POOL; // [NK*GlobalV::NBANDS*GlobalV::NLOCAL]


    void allocate_k(const Grid_Technique& gt,
        std::vector<ModuleBase::ComplexMatrix> &wfc_k);

    //=========================================
    // Init Cij, make it satisfy 2 conditions:
    // (1) Unit
    // (2) Orthogonal <i|S|j>= \delta{ij}
    //=========================================
	// void init_Cij(const bool change_c = 1);

	// mohan move orb_con here, 2021-05-24 
	ORB_control orb_con;
	
	private:

	bool wfck_flag; 
	bool complex_flag;
	bool allocate_flag;

};

#endif
