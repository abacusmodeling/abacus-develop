#include "module_hsolver/kernels/math_kernel_op.h"

#include <iomanip>
#include <iostream>

namespace hsolver
{

// CPU specialization of actual computation.
template <typename FPTYPE> 
struct zdot_real_op<FPTYPE, psi::DEVICE_CPU> {
  FPTYPE operator() (
      const psi::DEVICE_CPU* d,
      const int& dim,
      const std::complex<FPTYPE>* psi_L,
      const std::complex<FPTYPE>* psi_R,
      const bool reduce) 
  {
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    // qianrui modify 2021-3-14
    // Note that  ddot_(2*dim,a,1,b,1) = REAL( zdotc_(dim,a,1,b,1) )
    const FPTYPE* pL = reinterpret_cast<const FPTYPE*>(psi_L);
    const FPTYPE* pR = reinterpret_cast<const FPTYPE*>(psi_R);
    FPTYPE result = BlasConnector::dot(2 * dim, pL, 1, pR, 1);
    if (reduce) {
      Parallel_Reduce::reduce_double_pool(result);
    }
    return result;
  }
};

template <typename FPTYPE> struct vector_div_constant_op<FPTYPE, psi::DEVICE_CPU>
{
    void operator()(const psi::DEVICE_CPU* d,
                    const int dim,
                    std::complex<FPTYPE>* result,
                    const std::complex<FPTYPE>* vector,
                    const FPTYPE constant)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096/sizeof(FPTYPE))
#endif
        for (int i = 0; i < dim; i++)
        {
            result[i] = vector[i] / constant;
        }
    }
};

template <typename FPTYPE> struct vector_mul_vector_op<FPTYPE, psi::DEVICE_CPU>
{
    void operator()(const psi::DEVICE_CPU* d,
                    const int& dim,
                    std::complex<FPTYPE>* result,
                    const std::complex<FPTYPE>* vector1,
                    const FPTYPE* vector2)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096/sizeof(FPTYPE))
#endif
        for (int i = 0; i < dim; i++)
        {
            result[i] = vector1[i] * vector2[i];
        }
    }
};

template <typename FPTYPE> struct vector_div_vector_op<FPTYPE, psi::DEVICE_CPU>
{
    void operator()(const psi::DEVICE_CPU* d,
                    const int& dim,
                    std::complex<FPTYPE>* result,
                    const std::complex<FPTYPE>* vector1,
                    const FPTYPE* vector2)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096/sizeof(FPTYPE))
#endif
        for (int i = 0; i < dim; i++)
        {
            result[i] = vector1[i] / vector2[i];
        }
    }
};

template <typename FPTYPE> struct constantvector_addORsub_constantVector_op<FPTYPE, psi::DEVICE_CPU>
{
    void operator()(const psi::DEVICE_CPU* d,
                    const int& dim,
                    std::complex<FPTYPE>* result,
                    const std::complex<FPTYPE>* vector1,
                    const FPTYPE constant1,
                    const std::complex<FPTYPE>* vector2,
                    const FPTYPE constant2)
    {
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096/sizeof(FPTYPE))
#endif
        for (int i = 0; i < dim; i++)
        {
            result[i] = vector1[i] * constant1 + vector2[i] * constant2;
        }
    }
};

template <typename FPTYPE>
struct scal_op<FPTYPE, psi::DEVICE_CPU> {
    void operator()(
            const psi::DEVICE_CPU * /*ctx*/,
            const int &N,
            const std::complex<FPTYPE> *alpha,
            std::complex<FPTYPE> *X,
            const int &incx)
    {
        BlasConnector::scal(N, *alpha, X, incx);
    }
};

template <typename FPTYPE>
struct gemv_op<FPTYPE, psi::DEVICE_CPU> {
    void operator()(
            const psi::DEVICE_CPU *d,
            const char &trans,
            const int &m,
            const int &n,
            const std::complex<FPTYPE> *alpha,
            const std::complex<FPTYPE> *A,
            const int &lda,
            const std::complex<FPTYPE> *X,
            const int &incx,
            const std::complex<FPTYPE> *beta,
            std::complex<FPTYPE> *Y,
            const int &incy)
    {
        BlasConnector::gemv(trans, m, n, *alpha, A, lda, X, incx, *beta, Y, incy);
    }
};

template <typename FPTYPE>
struct axpy_op<FPTYPE, psi::DEVICE_CPU> {
    void operator()(
            const psi::DEVICE_CPU * /*ctx*/,
            const int &dim,
            const std::complex<FPTYPE> *alpha,
            const std::complex<FPTYPE> *X,
            const int &incX,
            std::complex<FPTYPE> *Y,
            const int &incY)
    {
        BlasConnector::axpy(dim, *alpha, X, incX, Y, incY);
    }
};

template <typename FPTYPE>
struct gemm_op<FPTYPE, psi::DEVICE_CPU> {
    void operator()(
            const psi::DEVICE_CPU * /*ctx*/,
            const char &transa,
            const char &transb,
            const int &m,
            const int &n,
            const int &k,
            const std::complex<FPTYPE> *alpha,
            const std::complex<FPTYPE> *a,
            const int &lda,
            const std::complex<FPTYPE> *b,
            const int &ldb,
            const std::complex<FPTYPE> *beta,
            std::complex<FPTYPE> *c,
            const int &ldc)
    {
        BlasConnector::gemm(transb, transa, n, m, k, *alpha, b, ldb, a, lda, *beta, c, ldc);
    }
};

template <typename FPTYPE> struct matrixTranspose_op<FPTYPE, psi::DEVICE_CPU>
{
    void operator()(const psi::DEVICE_CPU* d,
                    const int& row,
                    const int& col,
                    const std::complex<FPTYPE>* input_matrix,
                    std::complex<FPTYPE>* output_matrix)
    {
        std::complex<FPTYPE>* temp = nullptr;
        psi::memory::resize_memory_op<std::complex<FPTYPE>, psi::DEVICE_CPU>()(d, temp, row * col, "MTransOp");
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 4096/sizeof(FPTYPE))
#endif
        for (int j = 0; j < col; j++)
        {
            for (int i = 0; i < row; i++)
            {
                temp[j * row + i] = input_matrix[i * col + j];
            }
        }
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 4096/sizeof(FPTYPE))
#endif
        for (int i = 0; i < row * col; i++)
        {
            output_matrix[i] = temp[i];
        }
        psi::memory::delete_memory_op<std::complex<FPTYPE>, psi::DEVICE_CPU>()(d, temp);
    }
};

template <typename FPTYPE> struct matrixSetToAnother<FPTYPE, psi::DEVICE_CPU>
{
    void operator()(const psi::DEVICE_CPU* d,
                    const int& n,
                    const std::complex<FPTYPE>* A,
                    const int& LDA,
                    std::complex<FPTYPE>* B,
                    const int& LDB)
    {
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static, 4096/sizeof(FPTYPE))
#endif
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < LDA; j++)
            {
                B[i * LDB + j] = A[i * LDA + j];
            }
        }
    }
};



// Explicitly instantiate functors for the types of functor registered.
template struct scal_op<float, psi::DEVICE_CPU>;
template struct axpy_op<float, psi::DEVICE_CPU>;
template struct gemv_op<float, psi::DEVICE_CPU>;
template struct gemm_op<float, psi::DEVICE_CPU>;
template struct zdot_real_op<float, psi::DEVICE_CPU>;
template struct vector_div_constant_op<float, psi::DEVICE_CPU>;
template struct vector_mul_vector_op<float, psi::DEVICE_CPU>;
template struct vector_div_vector_op<float, psi::DEVICE_CPU>;
template struct constantvector_addORsub_constantVector_op<float, psi::DEVICE_CPU>;
template struct matrixTranspose_op<float, psi::DEVICE_CPU>;
template struct matrixSetToAnother<float, psi::DEVICE_CPU>;

template struct scal_op<double, psi::DEVICE_CPU>;
template struct axpy_op<double, psi::DEVICE_CPU>;
template struct gemv_op<double, psi::DEVICE_CPU>;
template struct gemm_op<double, psi::DEVICE_CPU>;
template struct zdot_real_op<double, psi::DEVICE_CPU>;
template struct vector_div_constant_op<double, psi::DEVICE_CPU>;
template struct vector_mul_vector_op<double, psi::DEVICE_CPU>;
template struct vector_div_vector_op<double, psi::DEVICE_CPU>;
template struct constantvector_addORsub_constantVector_op<double, psi::DEVICE_CPU>;
template struct matrixTranspose_op<double, psi::DEVICE_CPU>;
template struct matrixSetToAnother<double, psi::DEVICE_CPU>;

} // namespace hsolver