#ifndef HAMILT_H
#define HAMILT_H

#include "tools.h"
#ifdef __CUDA
#include "hamilt_pw.cuh"
#else
#include "hamilt_pw.h"
#endif
class Hamilt
{
	public:

    Hamilt();
    ~Hamilt();

    void init_before_ions(void);

    void init_before_electrons(void);

    void clear_after_ions(void);

	// diagH_subspace then full space in pw
    void diagH_pw(
		const int &istep,
		const int &iter,
		const int &ik,
		const double *precondition,
		double &avg_iter);

	// generate H and S then call diagH_subspace
    void diagH_subspace(
		const int ik,
		const int nstart,
		const int nbnd,
		const ModuleBase::ComplexMatrix &psi,
		ModuleBase::ComplexMatrix &evc,
		double *en);

	// be called by diagH_subspace
    void diagH_LAPACK(
		const int n,
		const int m,
		const ModuleBase::ComplexMatrix &hc, // Hamiltonian matrix
		const ModuleBase::ComplexMatrix &sc, // overlap matrix
		const int ldh,
		double *e, // output: eigenvalues
		ModuleBase::ComplexMatrix &hvec); // output: eigenvectors

    Hamilt_PW hpw;

private:

    bool test_exit_cond( const int &ntry, const int &notconv);

};
#endif
