#include "nonlocal_pw.h"

#include "module_base/blas_connector.h"
#include "module_base/timer.h"
#include "src_parallel/parallel_reduce.h"
#include "module_base/tool_quit.h"
#include "module_psi/include/device.h"

using hamilt::Nonlocal;
using hamilt::OperatorPW;

template<typename FPTYPE, typename Device>
Nonlocal<OperatorPW<FPTYPE, Device>>::Nonlocal(
    const int* isk_in,
    const pseudopot_cell_vnl* ppcell_in,
    const UnitCell_pseudo* ucell_in
)
{
    this->cal_type = pw_nonlocal;
    this->isk = isk_in;
    this->ppcell = ppcell_in;
    this->ucell = ucell_in;
    this->deeq = psi::device::get_device_type<Device>(this->ctx) == psi::GpuDevice ?
            this->ppcell->d_deeq :  // for GpuDevice
            this->ppcell->deeq.ptr; // for CpuDevice
    if( this->isk == nullptr || this->ppcell == nullptr || this->ucell == nullptr)
    {
        ModuleBase::WARNING_QUIT("NonlocalPW", "Constuctor of Operator::NonlocalPW is failed, please check your code!");
    }
}

template<typename FPTYPE, typename Device>
Nonlocal<OperatorPW<FPTYPE, Device>>::~Nonlocal() {
    delete_memory_op()(this->ctx, this->ps);
    delete_memory_op()(this->ctx, this->becp);
}

template<typename FPTYPE, typename Device>
void Nonlocal<OperatorPW<FPTYPE, Device>>::init(const int ik_in)
{
    this->ik = ik_in;
    // Calculate nonlocal pseudopotential vkb
	if(this->ppcell->nkb > 0) //xiaohui add 2013-09-02. Attention...
	{
		this->ppcell->getvnl(this->ik, this->ppcell->vkb);
	}

    if(this->next_op != nullptr)
    {
        this->next_op->init(ik_in);
    }
}

//--------------------------------------------------------------------------
// this function sum up each non-local pseudopotential located on each atom,
//--------------------------------------------------------------------------
template<typename FPTYPE, typename Device>
void Nonlocal<OperatorPW<FPTYPE, Device>>::add_nonlocal_pp(std::complex<FPTYPE> *hpsi_in, const std::complex<FPTYPE> *becp, const int m) const
{
    ModuleBase::timer::tick("Nonlocal", "add_nonlocal_pp");

    // number of projectors
    int nkb = this->ppcell->nkb;

    // std::complex<FPTYPE> *ps = new std::complex<FPTYPE>[nkb * m];
    // ModuleBase::GlobalFunc::ZEROS(ps, m * nkb);
    resize_memory_op()(this->ctx, this->ps, nkb * m);
    set_memory_op()(this->ctx, this->ps, 0, nkb * m);

    int sum = 0;
    int iat = 0;
    if (this->npol == 1)
    {
        const int current_spin = this->isk[this->ik];
        for (int it = 0; it < this->ucell->ntype; it++)
        {
            const int nproj = this->ucell->atoms[it].nh;
            // denghui replace 2022-10-20
            // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            nonlocal_op()(
                this->ctx,   // device context
                this->ucell->atoms[it].na, m, nproj, // four loop size
                sum, iat, current_spin, nkb,   // additional index params
                this->ppcell->deeq.getBound2(), this->ppcell->deeq.getBound3(), this->ppcell->deeq.getBound4(), // realArray operator()
                this->deeq, // array of data
                this->ps, this->becp); //  array of data
            // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            // for (int ia = 0; ia < this->ucell->atoms[it].na; ia++)
            // {
            //     // each atom has nproj, means this is with structure factor;
            //     // each projector (each atom) must multiply coefficient
            //     // with all the other projectors.
            //     for (int ib = 0; ib < m; ++ib)
            //     {
            //         for (int ip2 = 0; ip2 < nproj; ip2++)
            //         {
            //             for (int ip = 0; ip < nproj; ip++)
            //             {
            //                 this->ps[(sum + ip2) * m + ib]
            //                     += this->ppcell->deeq(current_spin, iat, ip, ip2) * this->becp[ib * nkb + sum + ip];
            //             } // end ib
            //         } // end ih
            //     } // end jh
            //     sum += nproj;
            //     ++iat;
            // } // end na
        } // end nt
    }
    else
    {
#if defined(__CUDA) || defined(__ROCM)
        ModuleBase::WARNING_QUIT("NonlocalPW", " gpu implementation of this->npol != 1 is not supported currently !!! ");
#endif
        for (int it = 0; it < this->ucell->ntype; it++)
        {
            int psind = 0;
            int becpind = 0;
            std::complex<FPTYPE> becp1 = std::complex<FPTYPE>(0.0, 0.0);
            std::complex<FPTYPE> becp2 = std::complex<FPTYPE>(0.0, 0.0);

            const int nproj = this->ucell->atoms[it].nh;
            for (int ia = 0; ia < this->ucell->atoms[it].na; ia++)
            {
                // each atom has nproj, means this is with structure factor;
                // each projector (each atom) must multiply coefficient
                // with all the other projectors.
                for (int ip = 0; ip < nproj; ip++)
                {
                    for (int ip2 = 0; ip2 < nproj; ip2++)
                    {
                        for (int ib = 0; ib < m; ib+=2)
                        {
                            psind = (sum + ip2) * m + ib;
                            becpind = ib * nkb + sum + ip;
                            becp1 = becp[becpind];
                            becp2 = becp[becpind + nkb];
                            ps[psind] += this->ppcell->deeq_nc(0, iat, ip2, ip) * becp1
                                         + this->ppcell->deeq_nc(1, iat, ip2, ip) * becp2;
                            ps[psind + 1] += this->ppcell->deeq_nc(2, iat, ip2, ip) * becp1
                                             + this->ppcell->deeq_nc(3, iat, ip2, ip) * becp2;
                        } // end ib
                    } // end ih
                } // end jh
                sum += nproj;
                ++iat;
            } // end na
        } // end nt
    }

    // use simple method.
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    // qianrui optimize 2021-3-31
    char transa = 'N';
    char transb = 'T';
    if (m == 1)
    {
        int inc = 1;
        // denghui replace 2022-10-20
        // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        gemv_op()(this->ctx,
                  transa,
                  this->npw,
                  this->ppcell->nkb,
                  &ModuleBase::ONE,
                  this->ppcell->vkb.c,
                  this->ppcell->vkb.nc,
                  this->ps,
                  inc,
                  &ModuleBase::ONE,
                  hpsi_in,
                  inc);
        // zgemv_(&transa,
        //        &this->npw,
        //        &(this->ppcell->nkb),
        //        &ModuleBase::ONE,
        //        this->ppcell->vkb.c,
        //        &this->ppcell->vkb.nc,
        //        ps,
        //        &inc,
        //        &ModuleBase::ONE,
        //        hpsi_in,
        //        &inc);
    }
    else
    {
        int npm = m;
        //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        // denghui replace 2022-10-20
        gemm_op()(
            this->ctx,
            transa,
            transb,
            this->npw,
            npm,
            this->ppcell->nkb,
            &ModuleBase::ONE,
            this->ppcell->vkb.c,
            this->ppcell->vkb.nc,
            this->ps,
            npm,
            &ModuleBase::ONE,
            hpsi_in,
            this->max_npw
        );
        // zgemm_(&transa,
        //        &transb,
        //        &this->npw,
        //        &npm,
        //        &(this->ppcell->nkb),
        //        &ModuleBase::ONE,
        //        this->ppcell->vkb.c,
        //        &this->ppcell->vkb.nc,
        //        ps,
        //        &npm,
        //        &ModuleBase::ONE,
        //        hpsi_in,
        //        &this->max_npw);
    }
    ModuleBase::timer::tick("Nonlocal", "add_nonlocal_pp");
}

template<typename FPTYPE, typename Device>
void Nonlocal<OperatorPW<FPTYPE, Device>>::act
(
    const psi::Psi<std::complex<FPTYPE>, Device>* psi_in, 
    const int n_npwx, 
    const std::complex<FPTYPE>* tmpsi_in, 
    std::complex<FPTYPE>* tmhpsi)const
{
    ModuleBase::timer::tick("Operator", "NonlocalPW");
    this->npw = psi_in->get_ngk(this->ik);
    this->max_npw = psi_in->get_nbasis() / psi_in->npol;
    this->npol = psi_in->npol;

    if (this->ppcell->nkb > 0)
    {
        //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        // qianrui optimize 2021-3-31
        int nkb = this->ppcell->nkb;
        resize_memory_op()(this->ctx, this->becp, n_npwx * nkb);
        // ModuleBase::ComplexMatrix becp(n_npwx, nkb, false);
        char transa = 'C';
        char transb = 'N';
        if (n_npwx == 1)
        {
            int inc = 1;
            // denghui replace 2022-10-20
            // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            gemv_op()(this->ctx,
                      transa,
                      this->npw,
                      nkb,
                      &ModuleBase::ONE,
                      this->ppcell->vkb.c,
                      this->ppcell->vkb.nc,
                      tmpsi_in,
                      inc,
                      &ModuleBase::ZERO,
                      this->becp,
                      inc);
            // zgemv_(&transa,
            //        &this->npw,
            //        &nkb,
            //        &ModuleBase::ONE,
            //        this->ppcell->vkb.c,
            //        &this->ppcell->vkb.nc,
            //        tmpsi_in,
            //        &inc,
            //        &ModuleBase::ZERO,
            //        this->becp,
            //        &inc);

        }
        else
        {
            int npm = n_npwx;
            //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            // denghui replace 2022-10-20
            gemm_op()(
                this->ctx,
                transa,
                transb,
                nkb,
                npm,
                this->npw,
                &ModuleBase::ONE,
                this->ppcell->vkb.c,
                this->ppcell->vkb.nc,
                tmpsi_in,
                this->max_npw,
                &ModuleBase::ZERO,
                this->becp,
                nkb
            );
            //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            // zgemm_(&transa,
            //        &transb,
            //        &nkb,
            //        &npm,
            //        &this->npw,
            //        &ModuleBase::ONE,
            //        this->ppcell->vkb.c,
            //        &this->ppcell->vkb.nc,
            //        tmpsi_in,
            //        &this->max_npw,
            //        &ModuleBase::ZERO,
            //        becp,
            //        &nkb);
            //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        }

        Parallel_Reduce::reduce_complex_double_pool(becp, nkb * n_npwx);

        this->add_nonlocal_pp(tmhpsi, becp, n_npwx);
    }
    ModuleBase::timer::tick("Operator", "NonlocalPW");
}

namespace hamilt{
template class Nonlocal<OperatorPW<double, psi::DEVICE_CPU>>;
#if ((defined __CUDA) || (defined __ROCM))
template class Nonlocal<OperatorPW<double, psi::DEVICE_GPU>>;
#endif
} // namespace hamilt