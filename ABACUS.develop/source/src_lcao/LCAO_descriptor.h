#ifndef LCAO_DESCRIPTOR_H
#define LCAO_DESCRIPTOR_H

#include "../src_global/intarray.h"
#include "../src_global/complexmatrix.h"

//caoyu add 2021-03-29
class LCAO_Descriptor
{
public:

    LCAO_Descriptor(int lm, int inlm);
    ~LCAO_Descriptor();

	// cal S_alpha_mu: overlap between lcao basis Phi and descriptor basis Al
    void build_S_descriptor(const bool &calc_deri); 

	// cal PDM: S_alpha_mu * inv(Sloc) * DM * inv(Sloc) * S_nu_beta
    void cal_projected_DM(void);

	// cal d: EIGENVALUE of PDM in block of I_n_l
    void cal_descriptor(void);
    void print_descriptor(void);

private:

	// overlap between lcao and descriptor basis
    double **S_mu_alpha;	//[Alpha.totnchi][NLOCAL*(2*lmax+1)]	caoyu modified 2021-05-07

	// projected density matrix
    double **PDM;	//[2l+1][2l+1]	caoyu modified 2021-05-07

	// descriptors
    double *d;

    int n_descriptor;

	// \sum_L{Nchi(L)*(2L+1)}
	int des_per_atom;

	const int lmaxd;

	const int inlmax;

	IntArray* mu_index;
	IntArray* inl_index;	//caoyu add 2021-05-07
	int* inl_ld;	//inl_ld[inl_index] = l_descriptor

    void init_mu_index(void);
    
	void set_S_mu_alpha(
		const int &iw1_all, 
		const int& inl,
		const int& im,
		const double& v);

    void print_projected_DM(
		ofstream &ofs, 
		ComplexMatrix &des, 
		const int &it, 
		const int &ia, 
		const int &l, 
		const int &n);
};

#endif
