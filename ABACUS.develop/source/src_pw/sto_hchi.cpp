#include "tools.h"
#include "global.h"
#include "sto_hchi.h" 

int Stochastic_hchi:: nrxx;
int Stochastic_hchi:: nx,Stochastic_hchi::ny,Stochastic_hchi::nz;
fftw_plan Stochastic_hchi:: pf, Stochastic_hchi::pb;
double Stochastic_hchi:: Emin, Stochastic_hchi:: Emax;
bool Stochastic_hchi:: initplan, Stochastic_hchi::ortho;
complex<double>* Stochastic_hchi:: rp_chi, * Stochastic_hchi::tmpchi2, * Stochastic_hchi::chig;


Stochastic_hchi::Stochastic_hchi()
{
	initplan = false;
	ortho = false;
	nrxx = 0;
	rp_chi = new complex<double> [1];
	tmpchi2 = new complex<double> [1];
	chig = new complex<double> [1];
}

Stochastic_hchi::~Stochastic_hchi()
{
	if(initplan)
	{
		fftw_destroy_plan(pf);
		fftw_destroy_plan(pb);
	}
	delete[] rp_chi;
	delete[] tmpchi2;
}

void Stochastic_hchi:: init()
{
	//wait for init--------------------------------------
		//nrxx
	//---------------------------------------------------
    if(nrxx != 0)
    {
        delete[] rp_chi;
        delete[] tmpchi2;
        rp_chi = new complex<double> [nrxx];
		tmpchi2 = new complex<double> [nrxx];
		pb=fftw_plan_dft_3d(nx,ny,nz,(fftw_complex *)rp_chi,(fftw_complex *)rp_chi, FFTW_BACKWARD, FFTW_MEASURE);
		pf=fftw_plan_dft_3d(nx,ny,nz,(fftw_complex *)tmpchi2,(fftw_complex *)tmpchi2, FFTW_FORWARD, FFTW_MEASURE);
		initplan = true;
    }
    else
    {
        WARNING_QUIT("Stochastic_hchi", "Number of grids should be at least one!");
    }

}

void Stochastic_hchi::orthogonal_to_psi(complex<double> *wfin, complex<double> *wfgortho)
{
	//wait for init
	int * GR_INDEX;
	double * v;

	TITLE("Stochastic_hchi","orthogonal_to_psi0");
	if(!initplan) WARNING_QUIT("Stochastic_hchi", "Please init hchi first!");

	for( int ir = 0; ir < nrxx; ++ir)
	{
		rp_chi[ir] = wfin[ir];
	}
	ZEROS(wfgortho,nrxx);
	fftw_execute(pb);
	
	delete []chig;
	chig = new complex<double> [wf.npw];
	for(int ig = 0; ig < wf.npw; ++ig)
	{
		chig[ig] = rp_chi[GR_INDEX[ig]]; 
	}

	//orthogonal part
	complex<double> sum;
	for(int iksb = 0; iksb < NBANDS; ++iksb)
	{
		complex<double> *kswf; // kswf is normalized
		sum=0;
		for(int ig = 0; ig < pw.ngmw; ++ig)
		{
			sum += conj(kswf[ig]) * chig[ig];
		}
		for(int ig = 0; ig < pw.ngmw; ++ig)
		{
			chig[ig] -= sum*kswf[ig];
		}
	}

	for(int ig = 0; ig < wf.npw; ++ig)
	{
		rp_chi[GR_INDEX[ig]] = chig[ig];
	}
	for(int ir = 0; ir < nrxx; ++ir)
	{
		wfgortho[ir] = rp_chi[ir];
	}
	ortho = true;
	return;
}




void Stochastic_hchi:: hchi(complex<double>*wfgortho, complex<double> *wfout)
{
	//wait for init--------------------------------------
		double dk1=1,dk2=1,dk3=1; double*vr;//vr= pot.vrs1 temporarily use cutoff vr.
	int * GR_INDEX;
		
	//---------------------------------------------------
	if(!initplan||!ortho) WARNING_QUIT("Stochastic_hchi", "Please init hchi first!");

	
	//------------------------------------
	//(1) the local potential.
	//------------------------------------
	if(VL_IN_H)
	{
		for(int ir = 0; ir < nrxx; ++ir)
		{
			tmpchi2[ir] = wfgortho[ir];
		}
		fftw_execute(pf);
		for(int ir = 0; ir < nrxx; ++ir)
		{
			wfout[ir] += tmpchi2[ir] * vr[ir];
		}
	}
	
	
	//------------------------------------
	//(2) the kinetical energy.
	//------------------------------------
	if(T_IN_H)
	{
		Vector3<double> gg;
		int gx,gy,gz;
		for(int ig1 = 0, i = 0; ig1 < nx; ++ig1)
		{
			for(int ig2 = 0; ig2 < ny; ++ig2)
			{
				for(int ig3 = 0; ig3 < nz; ++ig3)
				{
					gx = ig1;
					gy = ig2;
					gz = ig3;
					if(ig1 > nx/2) gx -= nx;
					if(ig2 > ny/2) gy -= ny;
					if(ig3 > nz/2) gz -= nz;
					gg.set((gx-double(nx-1)/2)*dk1, (gy-double(ny-1)/2)*dk2, (gz-double(nz-1)/2)*dk3);
					tmpchi2[i] = -gg.norm2()*wfgortho[i];
					++i;
				}
			}
		
		}
	}

	

	//------------------------------------
	// (3) the nonlocal pseudopotential.
	//------------------------------------
	if(VNL_IN_H)
	{
		if ( ppcell.nkb > 0)
		{
			complex<double> *becp = new complex<double>[ ppcell.nkb * NPOL ];
			ZEROS(becp,ppcell.nkb * NPOL);
			
			for (int i=0;i< ppcell.nkb;++i)
			{
				const complex<double>* p = &ppcell.vkb(i,0);
				const complex<double>* const p_end = p + wf.npw;
				for (int ig=0; ig< wf.npw; ++ig)
				{
					if(NSPIN!=4) becp[i] += chig[ig] * conj( p[ig] );
					else
					{
						//We didnot consider it temporarily.
					}	
				} 
			}

			//Parallel_Reduce::reduce_complex_double_pool( becp, ppcell.nkb * NPOL);
			complex<double> * Ps;
			ZEROS( Ps, ppcell.nkb * NPOL );
			int sum = 0;
    		int iat = 0;
    		// this function sum up each non-local pseudopotential located in each atom,
    		// all we need to do is put the right Dij coefficient to each becp, which
    		// is calculated before.
    		for (int it=0; it<ucell.ntype; ++it)
    		{
    		    const int Nprojs = ucell.atoms[it].nh;
    		    for (int ia=0; ia<ucell.atoms[it].na; ++ia)
    		    {
    		        // each atom has Nprojs, means this is with structure factor;
    		        // each projector (each atom) must multiply coefficient
    		        // with all the other projectors.
    		        for (int ip=0; ip<Nprojs; ++ip)
    		        {
    		            for (int ip2=0; ip2<Nprojs; ++ip2)
    		            {
							if(NSPIN!=4)
								Ps[sum+ip2] += ppcell.deeq(CURRENT_SPIN, iat, ip, ip2) * becp[sum+ip];
    		            }// end ih
    		        }//end jh
					if(NSPIN!=4) sum += Nprojs;
					++iat;
    		    } //end na
    		} //end nt

			// use simple method.
			ZEROS(chig, wf.npw);
			if(NSPIN!=4)
				for(int i=0; i<ppcell.nkb; ++i)
				{
					complex<double>* p = &ppcell.vkb(i,0);
					for(int ig=0; ig< wf.npw; ++ig)
					{
						chig[ig] = Ps[i] * p[ig];
					}
				}
			delete[] becp;
		}
	}


	//------------------------------------
	// (4) Conver (2) & (3) in Reciprocal space to Real one
	//------------------------------------
	for(int ig = 0; ig < wf.npw; ++ig)
	{
		tmpchi2[GR_INDEX[ig]] += chig[ig];
	}

	fftw_execute(pf);
	for(int i = 0; i<nrxx; ++i)
	{
		wfout[i] += tmpchi2 [i];
	}
	double Ebar = (Emin + Emax)/2;
	double DeltaE = (Emax - Emin)/2;
	for(int i = 0; i < nrxx; ++i)
	{
		wfout[i] = (wfout[i] - Ebar) / DeltaE;
	}
	
	delete [] chig;
	return;
}
