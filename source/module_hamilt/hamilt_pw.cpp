#include "hamilt_pw.h"

#include "module_base/blas_connector.h"
#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/timer.h"
#include "src_parallel/parallel_reduce.h"
#include "src_pw/global.h"
#include "src_pw/myfunc.h"

namespace hamilt
{

void HamiltPW::updateHk(const int ik)
{
    ModuleBase::TITLE("HamiltPW","updateHk");

	// mohan add 2010-09-30
	// (1) Which spin to use.
	if(GlobalV::NSPIN==2)
	{
		this->current_spin = GlobalC::kv.isk[ik];
        //GlobalV::CURRENT_SPIN = this->current_spin;
	}

	this->current_ik = ik;

	// (2) Take the local potential.
    this->current_veff = &GlobalC::pot.vr_eff(this->current_spin, 0);
	/*for (int ir=0; ir<GlobalC::rhopw->nrxx; ir++)
	{
		GlobalC::pot.vr_eff1[ir] = GlobalC::pot.vr_eff(GlobalV::CURRENT_SPIN, ir);//mohan add 2007-11-12
	}*/

	// (3) Calculate nonlocal pseudopotential vkb
	//if (GlobalC::ppcell.nkb > 0 && !LINEAR_SCALING) xiaohui modify 2013-09-02
	if(GlobalC::ppcell.nkb > 0 && (GlobalV::BASIS_TYPE=="pw" || GlobalV::BASIS_TYPE=="lcao_in_pw")) //xiaohui add 2013-09-02. Attention...
	{
		GlobalC::ppcell.getvnl(ik, GlobalC::ppcell.vkb);
	}

	// (4) The number of wave functions.
	this->current_npw = GlobalC::kv.ngk[ik];
    //GlobalC::wf.npw = this->current_npw;
    this->max_npw = GlobalC::wf.npwx;

    return;
}

void HamiltPW::hPsi(const std::complex<double> *psi_in, std::complex<double> *hpsi, const size_t size) const
{
    ModuleBase::timer::tick("HamiltPW", "h_psi");
    int m = int(size / this->max_npw / GlobalV::NPOL);
    if(int(size - m * this->max_npw * GlobalV::NPOL) != 0)m++;
    int i = 0;
    int j = 0;
    int ig = 0;
    const int ik = this->current_ik;
    const double tpiba2 = GlobalC::ucell.tpiba2;
    const int npwx = this->max_npw;
    const int npw = this->current_npw;
    const int nrxx = GlobalC::rhopw->nrxx;

    // if(GlobalV::NSPIN!=4) ModuleBase::GlobalFunc::ZEROS(hpsi, npw);
    // else ModuleBase::GlobalFunc::ZEROS(hpsi, this->max_npw * GlobalV::NPOL);//added by zhengdy-soc
    const int dmax = npwx * GlobalV::NPOL;

    //------------------------------------
    //(1) the kinetical energy.
    //------------------------------------
    std::complex<double> *tmhpsi;
    const std::complex<double> *tmpsi_in;
    if (GlobalV::T_IN_H)
    {
        tmhpsi = hpsi;
        tmpsi_in = psi_in;
        for (int ib = 0; ib < m; ++ib)
        {
            for (ig = 0; ig < npw; ++ig)
            {
                tmhpsi[ig] = GlobalC::wfcpw->getgk2(ik, ig) * tpiba2 * tmpsi_in[ig];
            }
            if (GlobalV::NSPIN == 4)
            {
                for (ig = npw; ig < npwx; ++ig)
                {
                    tmhpsi[ig] = 0;
                }
                tmhpsi += npwx;
                tmpsi_in += npwx;
                for (ig = 0; ig < npw; ++ig)
                {
                    tmhpsi[ig] = GlobalC::wfcpw->getgk2(ik, ig) * tpiba2 * tmpsi_in[ig];
                }
                for (ig = npw; ig < npwx; ++ig)
                {
                    tmhpsi[ig] = 0;
                }
            }
            tmhpsi += npwx;
            tmpsi_in += npwx;
        }
    }

    //------------------------------------
    //(2) the local potential.
    //-----------------------------------
    ModuleBase::timer::tick("HamiltPW", "vloc");
    std::complex<double> *porter = new std::complex<double>[GlobalC::wfcpw->nmaxgr];
    if (GlobalV::VL_IN_H)
    {
        tmhpsi = hpsi;
        tmpsi_in = psi_in;
        for (int ib = 0; ib < m; ++ib)
        {
            if (GlobalV::NSPIN != 4)
            {
                GlobalC::wfcpw->recip2real(tmpsi_in, porter, ik);
                for (int ir = 0; ir < nrxx; ++ir)
                {
                    porter[ir] *= this->current_veff[ir];
                }
                GlobalC::wfcpw->real2recip(porter, tmhpsi, ik, true);
            }
            else
            {
                std::complex<double> *porter1 = new std::complex<double>[nrxx];
                // fft to real space and doing things.
                GlobalC::wfcpw->recip2real(tmpsi_in, porter, ik);
                GlobalC::wfcpw->recip2real(tmpsi_in + npwx, porter1, ik);
                std::complex<double> sup, sdown;
                for (int ir = 0; ir < nrxx; ir++)
                {
                    sup = porter[ir] * (GlobalC::pot.vr_eff(0, ir) + GlobalC::pot.vr_eff(3, ir))
                          + porter1[ir]
                                * (GlobalC::pot.vr_eff(1, ir)
                                   - std::complex<double>(0.0, 1.0) * GlobalC::pot.vr_eff(2, ir));
                    sdown = porter1[ir] * (GlobalC::pot.vr_eff(0, ir) - GlobalC::pot.vr_eff(3, ir))
                            + porter[ir]
                                  * (GlobalC::pot.vr_eff(1, ir)
                                     + std::complex<double>(0.0, 1.0) * GlobalC::pot.vr_eff(2, ir));
                    porter[ir] = sup;
                    porter1[ir] = sdown;
                }
                // (3) fft back to G space.
                GlobalC::wfcpw->real2recip(porter, tmhpsi, ik, true);
                GlobalC::wfcpw->real2recip(porter1, tmhpsi + npwx, ik, true);

                delete[] porter1;
            }
            tmhpsi += dmax;
            tmpsi_in += dmax;
        }
    }
    ModuleBase::timer::tick("HamiltPW", "vloc");

    //------------------------------------
    // (3) the nonlocal pseudopotential.
    //------------------------------------
    ModuleBase::timer::tick("HamiltPW", "vnl");
    if (GlobalV::VNL_IN_H)
    {
        if (GlobalC::ppcell.nkb > 0)
        {
            //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            // qianrui optimize 2021-3-31
            int nkb = GlobalC::ppcell.nkb;
            ModuleBase::ComplexMatrix becp(GlobalV::NPOL * m, nkb, false);
            char transa = 'C';
            char transb = 'N';
            if (m == 1 && GlobalV::NPOL == 1)
            {
                int inc = 1;
                zgemv_(&transa,
                       &npw,
                       &nkb,
                       &ModuleBase::ONE,
                       GlobalC::ppcell.vkb.c,
                       &npwx,
                       psi_in,
                       &inc,
                       &ModuleBase::ZERO,
                       becp.c,
                       &inc);
            }
            else
            {
                int npm = GlobalV::NPOL * m;
                zgemm_(&transa,
                       &transb,
                       &nkb,
                       &npm,
                       &npw,
                       &ModuleBase::ONE,
                       GlobalC::ppcell.vkb.c,
                       &npwx,
                       psi_in,
                       &npwx,
                       &ModuleBase::ZERO,
                       becp.c,
                       &nkb);
            }

            Parallel_Reduce::reduce_complex_double_pool(becp.c, nkb * GlobalV::NPOL * m);

            this->add_nonlocal_pp(hpsi, becp.c, m);
        }
    }
    ModuleBase::timer::tick("HamiltPW", "vnl");
    //------------------------------------
    // (4) the metaGGA part
    //------------------------------------
    ModuleBase::timer::tick("HamiltPW", "meta");
    if (XC_Functional::get_func_type() == 3)
    {
        tmhpsi = hpsi;
        tmpsi_in = psi_in;
        for (int ib = 0; ib < m; ++ib)
        {
            for (int j = 0; j < 3; j++)
            {
                for (int ig = 0; ig < GlobalC::kv.ngk[this->current_ik]; ig++)
                {
                    double fact = GlobalC::wfcpw->getgpluskcar(ik, ig)[j] * tpiba2;
                    porter[ig] = tmpsi_in[ig] * complex<double>(0.0, fact);
                }

                GlobalC::wfcpw->recip2real(porter, porter, ik);

                for (int ir = 0; ir < nrxx; ir++)
                {
                    porter[ir] *= GlobalC::pot.vofk(GlobalV::CURRENT_SPIN, ir);
                }
                GlobalC::wfcpw->real2recip(porter, porter, ik);

                for (int ig = 0; ig < npw; ig++)
                {
                    double fact = GlobalC::wfcpw->getgpluskcar(ik, ig)[j] * tpiba2;
                    tmhpsi[ig] -= complex<double>(0.0, fact) * porter[ig];
                }
            } // x,y,z directions
        }
    }
    delete[] porter;
    ModuleBase::timer::tick("HamiltPW", "meta");
    ModuleBase::timer::tick("HamiltPW", "h_psi");
    return;
}//end hPsi

void HamiltPW::sPsi
(
    const std::complex<double> *psi,
    std::complex<double> *spsi,
    size_t size
) const
{
    for (size_t i=0; i<size; i++)
    {
        spsi[i] = psi[i];
    }
    return;
}

//--------------------------------------------------------------------------
// this function sum up each non-local pseudopotential located on each atom,
//--------------------------------------------------------------------------
void HamiltPW::add_nonlocal_pp(
	std::complex<double> *hpsi_in,
	const std::complex<double> *becp,
	const int m) const
{
    ModuleBase::timer::tick("HamiltPW","add_nonlocal_pp");

	// number of projectors
	int nkb = GlobalC::ppcell.nkb;

	std::complex<double> *ps  = new std::complex<double> [nkb * GlobalV::NPOL * m];
    ModuleBase::GlobalFunc::ZEROS(ps, GlobalV::NPOL * m * nkb);

    int sum = 0;
    int iat = 0;
    if(GlobalV::NSPIN!=4)
	{
		for (int it=0; it<GlobalC::ucell.ntype; it++)
		{
			const int nproj = GlobalC::ucell.atoms[it].nh;
			for (int ia=0; ia<GlobalC::ucell.atoms[it].na; ia++)
			{
				// each atom has nproj, means this is with structure factor;
				// each projector (each atom) must multiply coefficient
				// with all the other projectors.
				for (int ip=0; ip<nproj; ip++)
				{
					for (int ip2=0; ip2<nproj; ip2++)
					{
						for(int ib = 0; ib < m ; ++ib)
						{
							ps[(sum + ip2) * m + ib] +=
							GlobalC::ppcell.deeq(GlobalV::CURRENT_SPIN, iat, ip, ip2)
							* becp[ib * nkb + sum + ip];
						}//end ib
					}// end ih
				}//end jh
				sum += nproj;
				++iat;
			} //end na
		} //end nt
	}
	else
	{
		for (int it=0; it<GlobalC::ucell.ntype; it++)
		{
			int psind=0;
			int becpind=0;
			std::complex<double> becp1=std::complex<double>(0.0,0.0);
			std::complex<double> becp2=std::complex<double>(0.0,0.0);

			const int nproj = GlobalC::ucell.atoms[it].nh;
			for (int ia=0; ia<GlobalC::ucell.atoms[it].na; ia++)
			{
				// each atom has nproj, means this is with structure factor;
				// each projector (each atom) must multiply coefficient
				// with all the other projectors.
				for (int ip=0; ip<nproj; ip++)
				{
					for (int ip2=0; ip2<nproj; ip2++)
					{
						for(int ib = 0; ib < m ; ++ib)
						{
							psind = (sum+ip2) * 2 * m + ib * 2;
							becpind = ib*nkb*2 + sum + ip;
							becp1 =  becp[becpind];
							becp2 =  becp[becpind + nkb];
							ps[psind] += GlobalC::ppcell.deeq_nc(0, iat, ip2, ip) * becp1
								+GlobalC::ppcell.deeq_nc(1, iat, ip2, ip) * becp2;
							ps[psind +1] += GlobalC::ppcell.deeq_nc(2, iat, ip2, ip) * becp1
								+GlobalC::ppcell.deeq_nc(3, iat, ip2, ip) * becp2;
						}//end ib
					}// end ih
				}//end jh
				sum += nproj;
				++iat;
			} //end na
		} //end nt
	}

	// use simple method.
	//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	//qianrui optimize 2021-3-31
	char transa = 'N';
	char transb = 'T';
	if(GlobalV::NPOL==1 && m==1)
	{
		int inc = 1;
		zgemv_(&transa,
			&this->current_npw,
			&GlobalC::ppcell.nkb,
			&ModuleBase::ONE,
			GlobalC::ppcell.vkb.c,
			&this->max_npw,
			ps,
			&inc,
			&ModuleBase::ONE,
			hpsi_in,
			&inc);
	}
	else
	{
		int npm = GlobalV::NPOL*m;
		zgemm_(&transa,
			&transb,
			&this->current_npw,
			&npm,
			&GlobalC::ppcell.nkb,
			&ModuleBase::ONE,
			GlobalC::ppcell.vkb.c,
			&this->max_npw,
			ps,
			&npm,
			&ModuleBase::ONE,
			hpsi_in,
			&this->max_npw);
	}

	delete[] ps;
    ModuleBase::timer::tick("HamiltPW","add_nonlocal_pp");
    return;
}

} // namespace hamilt