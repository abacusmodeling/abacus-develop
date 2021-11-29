#ifndef GLOBAL_H
#define GLOBAL_H

#include "../run_pw.h"
#include "../module_base/global_variable.h"
#include "../module_base/global_function.h"
#include "pw_basis.h"
#include "energy.h"
#include "VNL_in_pw.h"
#include "charge_broyden.h"
#include "potential.h"
#include "xc_type.h"
#include "hamilt.h"
#include "../src_ions/ions.h"
#include "wavefunc.h"
#include "use_fft.h"
#include "klist.h"
#include "../src_io/output.h"
#include "magnetism.h"
#include "vdwd2.h"
#include "vdwd2_parameters.h"
#include "vdwd3.h"
#include "vdwd3_parameters.h"	
#include "../src_io/restart.h" 
#include "exx_global.h"
#include "../src_lcao/exx_lip.h"
#include "../src_parallel/ft.h"

#ifdef __CUDA
static const char *_cublasGetErrorEnum(cublasStatus_t error)
{
    switch (error)
    {
        case CUBLAS_STATUS_SUCCESS:
            return "CUBLAS_STATUS_SUCCESS";
        case CUBLAS_STATUS_NOT_INITIALIZED:
            return "CUBLAS_STATUS_NOT_INITIALIZED";
        case CUBLAS_STATUS_ALLOC_FAILED:
            return "CUBLAS_STATUS_ALLOC_FAILED";
        case CUBLAS_STATUS_INVALID_VALUE:
            return "CUBLAS_STATUS_INVALID_VALUE";
        case CUBLAS_STATUS_ARCH_MISMATCH:
            return "CUBLAS_STATUS_ARCH_MISMATCH";
        case CUBLAS_STATUS_MAPPING_ERROR:
            return "CUBLAS_STATUS_MAPPING_ERROR";
        case CUBLAS_STATUS_EXECUTION_FAILED:
            return "CUBLAS_STATUS_EXECUTION_FAILED";
        case CUBLAS_STATUS_INTERNAL_ERROR:
            return "CUBLAS_STATUS_INTERNAL_ERROR";
    }
    return "<unknown>";
}

static const char *_cusolverGetErrorEnum(cusolverStatus_t error)
{
    switch (error)
    {
        case CUSOLVER_STATUS_SUCCESS:
            return "CUSOLVER_STATUS_SUCCESS";
        case CUSOLVER_STATUS_NOT_INITIALIZED:
            return "CUSOLVER_STATUS_NOT_INITIALIZED";
        case CUSOLVER_STATUS_ALLOC_FAILED:
            return "CUSOLVER_STATUS_ALLOC_FAILED";
        case CUSOLVER_STATUS_INVALID_VALUE:
            return "CUSOLVER_STATUS_INVALID_VALUE";
        case CUSOLVER_STATUS_ARCH_MISMATCH:
            return "CUSOLVER_STATUS_ARCH_MISMATCH";
        case CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED:
            return "CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED";
        case CUSOLVER_STATUS_EXECUTION_FAILED:
            return "CUSOLVER_STATUS_EXECUTION_FAILED";
        case CUSOLVER_STATUS_INTERNAL_ERROR:
            return "CUSOLVER_STATUS_INTERNAL_ERROR";
    }
    return "<unknown>";
}

static const char *_cufftGetErrorEnum(cufftResult_t error)
{
    switch (error)
    {
        case CUFFT_SUCCESS:
            return "CUFFT_SUCCESS";
        case CUFFT_INVALID_PLAN:
            return "CUFFT_INVALID_PLAN";
        case CUFFT_ALLOC_FAILED:
            return "CUFFT_ALLOC_FAILED";
        case CUFFT_INVALID_TYPE:
            return "CUFFT_INVALID_TYPE";
        case CUFFT_INVALID_VALUE:
            return "CUFFT_INVALID_VALUE";
        case CUFFT_INTERNAL_ERROR:
            return "CUFFT_INTERNAL_ERROR";
        case CUFFT_EXEC_FAILED:
            return "CUFFT_EXEC_FAILED";
        case CUFFT_SETUP_FAILED:
            return "CUFFT_SETUP_FAILED";
        case CUFFT_INVALID_SIZE:
            return "CUFFT_INVALID_SIZE";
        case CUFFT_UNALIGNED_DATA:
            return "CUFFT_UNALIGNED_DATA";
        case CUFFT_INCOMPLETE_PARAMETER_LIST:
            return "CUFFT_INCOMPLETE_PARAMETER_LIST";
        case CUFFT_INVALID_DEVICE:
            return "CUFFT_INVALID_DEVICE";
        case CUFFT_PARSE_ERROR:
            return "CUFFT_PARSE_ERROR";
        case CUFFT_NO_WORKSPACE:
            return "CUFFT_NO_WORKSPACE";
        case CUFFT_NOT_IMPLEMENTED:
            return "CUFFT_NOT_IMPLEMENTED";
        case CUFFT_LICENSE_ERROR:
            return "CUFFT_LICENSE_ERROR";
        case CUFFT_NOT_SUPPORTED:
            return "CUFFT_NOT_SUPPORTED";
    }
    return "<unknown>";
}

#define CHECK_CUDA(func)\
{\
    cudaError_t status = (func);\
    if(status != cudaSuccess)\
    {\
        printf("In File %s : CUDA API failed at line %d with error: %s (%d)\n",\
            __FILE__, __LINE__, cudaGetErrorString(status), status);\
    }\
}

#define CHECK_CUBLAS(func)\
{\
    cublasStatus_t status = (func);\
    if(status != CUBLAS_STATUS_SUCCESS)\
    {\
        printf("In File %s : CUBLAS API failed at line %d with error: %s (%d)\n",\
            __FILE__, __LINE__, _cublasGetErrorEnum(status), status);\
    }\
}

#define CHECK_CUSOLVER(func)\
{\
    cusolverStatus_t status = (func);\
    if(status != CUSOLVER_STATUS_SUCCESS)\
    {\
        printf("In File %s : CUSOLVER API failed at line %d with error: %s (%d)\n",\
            __FILE__, __LINE__, _cusolverGetErrorEnum(status), status);\
    }\
}

#define CHECK_CUFFT(func)\
{\
    cufftResult_t status = (func);\
    if(status != CUFFT_SUCCESS)\
    {\
        printf("In File %s : CUFFT API failed at line %d with error: %s (%d)\n",\
            __FILE__, __LINE__, _cufftGetErrorEnum(status), status);\
    }\
}
#endif


//==========================================================
// EXPLAIN : define "GLOBAL CLASS"
//==========================================================
namespace GlobalC
{
extern K_Vectors kv;
extern Use_FFT UFFT;
extern output out;
extern PW_Basis pw;
extern Stochastic_WF sto_wf; //qianrui add 2021-2-5
extern energy en;
extern wavefunc wf;
extern Hamilt hm;
extern Exx_Global exx_global;
extern Exx_Lip exx_lip;
extern pseudopot_cell_vnl ppcell;
}


#include "../module_symmetry/symmetry.h"
#include "../module_cell/unitcell_pseudo.h"
#include "../src_parallel/parallel_grid.h"
#include "../src_parallel/parallel_kpoints.h"
namespace GlobalC
{
extern UnitCell_pseudo ucell;
extern xcfunc xcf;
extern Charge_Broyden CHR;
extern Potential pot;
extern ModuleSymmetry::Symmetry symm;
extern Parallel_Grid Pgrid; 
extern Parallel_Kpoints Pkpoints;
extern Vdwd2_Parameters vdwd2_para;		// Peize Lin add 2021.03.09
extern Vdwd3_Parameters vdwd3_para;		// jiyy add 2021-05-02	
extern Restart restart;	// Peize Lin add 2020.04.04
}

//extern Magnetism mag;

#ifdef __LCAO
#include "../src_lcao/global_fp.h"
#endif

#endif
