#ifndef HAMILT_H
#define HAMILT_H

#include "tools.h"
#include "../module_ORB/ORB_control.h"
#include "hamilt_pw.h"

class Hamilt
{
	public:

    Hamilt();
    ~Hamilt();

    void init_before_ions(void);

    void init_before_electrons(void);

    void clear_after_ions(void);

    void diagH_subspace(
		const int ik, 
		const int nstart,
		const int nbnd,
		const ComplexMatrix &psi,
		ComplexMatrix &evc,
		double *en);

    void diago(const int &istep,const int &iter,const int &ik,
               const double *precondition,double &avg_iter);

    void diagH_LAPACK(
		const int n,
		const int m,
		const ComplexMatrix &hc, // Hamiltonian matrix
		const ComplexMatrix &sc, // overlap matrix
		const int ldh,
		double *e, // output: eigenvalues
		ComplexMatrix &hvec); // output: eigenvectors

    Hamilt_PW hpw;
	
	// mohan update 2021-02-10
	ORB_control orb_con;

private:

    bool test_exit_cond( const int &ntry, const int &notconv);

};
#endif
