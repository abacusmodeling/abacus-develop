#include "global.h"
#include "sto_hchi.h" 
#include "../module_base/tool_title.h"
#include "../module_base/timer.h"
#include "../src_parallel/parallel_reduce.h"


double Stochastic_hchi::Emin;
double Stochastic_hchi::Emax;



Stochastic_hchi::Stochastic_hchi()
{
	Stochastic_hchi:: Emin = INPUT.emin_sto;	
	Stochastic_hchi:: Emax = INPUT.emax_sto;
}

Stochastic_hchi::~Stochastic_hchi()
{
}

void Stochastic_hchi:: init()
{

}


void Stochastic_hchi:: hchi_reciprocal(complex<double> *chig, complex<double> *hchig, const int m)
{
	ModuleBase::timer::tick("Stochastic_hchi","hchi_reciprocal");
	
	//---------------------------------------------------

	int npwx = GlobalC::wf.npwx;
	int npw = GlobalC::wf.npw;
	int npm = GlobalV::NPOL * m;
	int inc = 1;
	//------------------------------------
	//(1) the kinetical energy.
	//------------------------------------
	complex<double> *chibg = chig;
	complex<double> *hchibg = hchig;
	if(GlobalV::T_IN_H)
	{
		for (int ib = 0; ib < m ; ++ib)
		{
			for (int ig = 0; ig < npw; ++ig)
			{
				hchibg[ig] = GlobalC::wf.g2kin[ig] * chibg[ig];
			}
			chibg += npwx;
			hchibg += npwx;
		}
	}
	
	//------------------------------------
	//(2) the local potential.
	//------------------------------------
	ModuleBase::timer::tick("Stochastic_hchi","vloc");
	if(GlobalV::VL_IN_H)
	{
		chibg = chig;
		hchibg = hchig;
		for(int ib = 0 ; ib < m ; ++ib)
		{
			ModuleBase::GlobalFunc::ZEROS( GlobalC::UFFT.porter, GlobalC::pw.nrxx);
			GlobalC::UFFT.RoundTrip( chibg, GlobalC::pot.vr_eff1, GlobalC::hm.hpw.GR_index, GlobalC::UFFT.porter );
			for (int ig = 0; ig < npw; ++ig)
			{
				hchibg[ig] += GlobalC::UFFT.porter[ GlobalC::hm.hpw.GR_index[ig] ];
			}
			chibg += npwx;
			hchibg += npwx;
		}
			
	}
	ModuleBase::timer::tick("Stochastic_hchi","vloc");


	//------------------------------------
	// (3) the nonlocal pseudopotential.
	//------------------------------------
	ModuleBase::timer::tick("Stochastic_hchi","vnl");
	if(GlobalV::VNL_IN_H)
	{
		if ( GlobalC::ppcell.nkb > 0)
		{
			int nkb = GlobalC::ppcell.nkb;
			complex<double> *becp = new complex<double>[ nkb * GlobalV::NPOL * m ];
			char transc = 'C';
			char transn = 'N';
			char transt = 'T';
			if(m==1 && GlobalV::NPOL ==1)
			{
				zgemv_(&transc, &npw, &nkb, &ModuleBase::ONE, GlobalC::ppcell.vkb.c, &GlobalC::wf.npwx, chig, &inc, &ModuleBase::ZERO, becp, &inc);
			}
			else
			{
				zgemm_(&transc,&transn,&nkb,&npm,&npw,&ModuleBase::ONE,GlobalC::ppcell.vkb.c,&GlobalC::wf.npwx,chig,&npwx,&ModuleBase::ZERO,becp,&nkb);
			}
			Parallel_Reduce::reduce_complex_double_pool( becp, nkb * GlobalV::NPOL * m);

			complex<double> *Ps  = new complex<double> [nkb * GlobalV::NPOL * m];
   			ModuleBase::GlobalFunc::ZEROS( Ps, GlobalV::NPOL * m * nkb);
			
			int sum = 0;
    		int iat = 0;
    		for (int it=0; it<GlobalC::ucell.ntype; it++)
    		{
    		    const int Nprojs = GlobalC::ucell.atoms[it].nh;
    		    for (int ia=0; ia<GlobalC::ucell.atoms[it].na; ia++)
    		    {
    		        // each atom has Nprojs, means this is with structure factor;
    		        // each projector (each atom) must multiply coefficient
    		        // with all the other projectors.
    		        for (int ip=0; ip<Nprojs; ip++)
    		        {
    		            for (int ip2=0; ip2<Nprojs; ip2++)
    		            {
							for(int ib = 0; ib < m ; ++ib)
							{
								Ps[(sum + ip2) * m + ib] += 
								GlobalC::ppcell.deeq(GlobalV::CURRENT_SPIN, iat, ip, ip2) * becp[ib * nkb + sum + ip];
							}//end ib
    		            }// end ih
    		        }//end jh 
					sum += Nprojs;
					++iat;
    		    } //end na
    		} //end nt

			if(GlobalV::NPOL==1 && m==1)
			{
				zgemv_(&transn, &npw, &nkb, &ModuleBase::ONE, 
						GlobalC::ppcell.vkb.c, &GlobalC::wf.npwx, Ps, &inc, &ModuleBase::ONE, hchig, &inc);
			}
			else
			{
				zgemm_(&transn,&transt,&npw,&npm,&nkb,&ModuleBase::ONE,
						GlobalC::ppcell.vkb.c,&GlobalC::wf.npwx,Ps,&npm,&ModuleBase::ONE,hchig,&npwx);
			}

			delete[] becp;
			delete[] Ps;
		}
	}
	ModuleBase::timer::tick("Stochastic_hchi","vnl");



	double Ebar = (Emin + Emax)/2;
	double DeltaE = (Emax - Emin)/2;
	for(int ib = 0 ; ib < m ; ++ib)
	{
		for(int ig = 0; ig < npw; ++ig)
		{
			hchig[ib*npwx+ig] = (hchig[ib*npwx+ig] - Ebar * chig[ib*npwx+ig]) / DeltaE;
		}
	}
	
	ModuleBase::timer::tick("Stochastic_hchi","hchi_reciprocal");


	return;
}
