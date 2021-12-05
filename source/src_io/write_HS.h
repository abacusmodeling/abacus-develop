#ifndef HS_MATRIX_H
#define HS_MATRIX_H

#include "../src_pw/tools.h"

// mohan add this file 2010-09-10
namespace HS_Matrix
{
    void saving_HS(const double *Hloc, const double* Sloc, bool bit, const int &out_hs);

    void save_HS(const double *H, const double *S, bool bit);

    void save_HS_complex(const std::complex<double> *H, const std::complex<double> *S, bool bit);

    void save_HSR_tr(const int current_spin); //LiuXh add 2019-07-15

    // jingan add 2021-6-4, modify 2021-12-2
    void save_HSR_sparse(
        const double &sparse_threshold, 
        const bool &binary,  
        const std::string &SR_filename, 
        const std::string &HR_filename_up, 
        const std::string &HR_filename_down
    );
    void save_SR_sparse(
        const double &sparse_threshold, 
        const bool &binary,  
        const std::string &SR_filename
    );
    void output_single_R(std::ofstream &ofs, const std::map<size_t, std::map<size_t, double>> &XR, const double &sparse_threshold, const bool &binary);
    void output_soc_single_R(std::ofstream &ofs, const std::map<size_t, std::map<size_t, std::complex<double>>> &XR, const double &sparse_threshold, const bool &binary);

// mohan comment out 2021-02-10
// void save_HS_ccf(const int &iter, const int &Hnnz, const int *colptr_H, const int *rowind_H, 
// const double *nzval_H, const double *nzval_S, bool bit);

    void saving_HS_complex(std::complex<double> *Hloc, std::complex<double>* Sloc, bool bit, const int &out_hs); //LiuXh, 2017-03-21

	void save_HS_complex(std::complex<double> *H, std::complex<double> *S, bool bit); //LiuXh, 2017-03-21

  void save_SR_sparse(const double &sparse_threshold, const bool &binary, const std::string file);
  
}

#endif
