#ifndef HAMILT_H
#define HAMILT_H

#include "tools.h"
#include "hamilt_pw.h"

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
		const ComplexMatrix &psi,
		ComplexMatrix &evc,
		double *en);

	// be called by diagH_subspace
    void diagH_LAPACK(
		const int n,
		const int m,
		const ComplexMatrix &hc, // Hamiltonian matrix
		const ComplexMatrix &sc, // overlap matrix
		const int ldh,
		double *e, // output: eigenvalues
		ComplexMatrix &hvec); // output: eigenvectors

    Hamilt_PW hpw;

private:

    bool test_exit_cond( const int &ntry, const int &notconv);

};
#endif
