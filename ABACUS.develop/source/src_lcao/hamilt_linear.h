//==========================================================
// AUTHOR : mohan, ywcui
// Last Update: 2009-4-22
//==========================================================
#ifndef HAMILT_LINEAR_H
#define HAMILT_LINEAR_H

#include "../src_pw/tools.h"
#include "../src_pw/algorithms.h"

#include "use_overlap_table.h" // (2)
#include "local_orbital_pairs.h" // (5) pairs basis
#include "../src_external/src_pdiag/pdiag_double.h"

class Hamilt_Linear
{
public:
    Hamilt_Linear();
    ~Hamilt_Linear();

    // Generate the S(overlap),T,NL matrix.
    void set_orb_tables();
    void clear_after_ions();

    // For test
    void solve_using_lapack
    (
        const int ik,
        const ComplexMatrix &lpsi,// localized psi
        const int nstart,
        const int nbands,
        ComplexMatrix &evc,// wave functions
        double *en// band energy
    );

    // Main routine
    void solve_using_cg
    (
        const int ik,
        const ComplexMatrix &lpsi,// localized psi
        const int nstart,
        const int nbands,
        ComplexMatrix &evc,// wave functions
        double *en// band energy
    );

	Local_Orbital_Pairs pairs;
};
#endif
