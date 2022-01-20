#ifndef PARALLEL_ORBITALS_H
#define PARALLEL_ORBITALS_H

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../src_pdiag/pdiag_double.h"

class Parallel_Orbitals : public Pdiag_Double
{
    public:

    Parallel_Orbitals();
    ~Parallel_Orbitals();

    // type : Sloc(1) Hloc(2) Hloc_fixed(3)
    bool in_this_processor(const int &iw1_all, const int &iw2_all);
    
    int* trace_loc_row;
    int* trace_loc_col;
    int out_hs; // mohan add 2010-09-02
    int out_hsR; // LiuXh add 2019-07-16

    void set_trace(void);

    // use MPI_Alltoallv to convert a orbital distributed matrix to 2D-block cyclic distributed matrix
    // here are the parameters which will be used
    int sender_index_size;
    int *sender_local_index;
    int sender_size;
    int *sender_size_process;
    int *sender_displacement_process;
    double* sender_buffer;

    int receiver_index_size;
    int *receiver_global_index;
    int receiver_size;
    int *receiver_size_process;
    int *receiver_displacement_process;
    double* receiver_buffer;
};

#endif
