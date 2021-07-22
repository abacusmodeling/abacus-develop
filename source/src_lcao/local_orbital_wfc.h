#ifndef LOCAL_ORBITAL_WFC
#define LOCAL_ORBITAL_WFC

#include "grid_technique.h"
#include "../src_pw/tools.h"
#include "../module_orbital/ORB_control.h" // mohan add 2021-05-24

class Local_Orbital_wfc
{
	public:

	Local_Orbital_wfc();
	~Local_Orbital_wfc();

	// used to generate density matrix: LOC.DM_R,
	// which is used to calculate the charge density. 
	// which is got after the diagonalization of 
	// complex Hamiltonian matrix.
	complex<double>*** WFC_K; // [NK, NBANDS, NLOCAL]	
	complex<double>* WFC_K_POOL; // [NK*NBANDS*NLOCAL]

	// augmented wave functions to 'c',
	// used to generate density matrix 
	// according to 2D data block.
	// mohan add 2010-09-26
	// daug means : dimension of augmented wave functions
	double*** WFC_GAMMA_aug; // [NSPIN, NBANDS, daug];
	complex<double>*** WFC_K_aug; // [NK, NBANDS, daug];
	int* trace_aug;
	
	// how many elements are missing. 
	int daug;

	void allocate_k(const Grid_Technique &gt);
	void set_trace_aug(const Grid_Technique &gt);
	bool get_allocate_aug_flag(void)const{return allocate_aug_flag;}

    //=========================================
    // Init Cij, make it satisfy 2 conditions:
    // (1) Unit
    // (2) Orthogonal <i|S|j>= \delta{ij}
    //=========================================
	// void init_Cij(const bool change_c = 1);
	bool get_allocate_flag(void)const{return allocate_flag;}	

	// mohan move orb_con here, 2021-05-24 
	ORB_control orb_con;
	
	private:

	bool wfck_flag; 
	bool complex_flag;
	bool allocate_flag;
	bool allocate_aug_flag;

};

#endif
