#include "tools.h"
#include "global.h"
#include "sto_hchi.h" 


Stochastic_Hchi::Stochastic_Hchi()
{
	initplan = false;
	initchi = false;
	nrxx = 0;
	tmpchi1 = new complex<double> [1];
	tmpchi2 = new complex<double> [1];
}

Stochastic_Hchi::~Stochastic_Hchi()
{
	if(initplan)
	{
		fftw_destroy_plan(pf);
	}
	delete[] tmpchi1;
	delete[] tmpchi2;
}

void Stochastic_Hchi:: init()
{
    if(nrxx != 0)
    {
        delete[] tmpchi1;
        delete[] tmpchi2;
        tmpchi1 = new complex<double> [nrxx];
		tmpchi2 = new complex<double> [nrxx];
        initchi = true;
    }
    else
    {
        WARNING_QUIT("Stochastic_Hchi", "Number of grids should be at least one!");
    }

}

void Stochastic_Hchi:: Hchi(complex<double>*wfin, complex<double> *wfout)
{
	//wait for init--------------------------------------
		double dk1=1,dk2=1,dk3=1; double*vr;
		
	//---------------------------------------------------
	if(!initchi) WARNING_QUIT("Stochastic_Hchi", "Please init Hchi first!");
	if(!initplan)
	{
		initplan=true;
		pf=fftw_plan_dft_3d(nx,ny,nz,(fftw_complex *)tmpchi1,(fftw_complex *)tmpchi2, FFTW_FORWARD, FFTW_MEASURE);
	}
	complex<double> ui(0,1);
	for(int ix = 0, i = 0; ix < nx; ++ix)
	{
		for(int iy = 0; iy < ny; ++iy)
		{
			for(int iz = 0; iz < nz; ++iz)
			{
				tmpchi1[i] = wfin[i]*exp(PI*(double(nx)/(nx-1)*ix+double(ny)/(ny-1)*iy+double(nz)/(nz-1)*iz)*ui);
				++i;
			}
		}
		
	}
	
	//------------------------------------
	//(1) the kinetical energy.
	//------------------------------------
	if(T_IN_H)
	{
		fftw_execute(pf);
		Vector3<double> gg;
		for(int ig1 = 0, i = 0; ig1 < nx; ++ig1)
		{
			for(int ig2 = 0; ig2 < ny; ++ig2)
			{
				for(int ig3 = 0; ig3 < nz; ++ig3)
				{
					gg.set((ig1-double(nx-1)/2)*dk1, (ig2-double(ny-1)/2)*dk2, (ig3-double(nz-1)/2)*dk3);
					tmpchi2[i] *= -gg.norm2();
					++i;
				}
			}
		
		}
	}

	//------------------------------------
	//(2) the local potential.
	//------------------------------------
	if(VL_IN_H)
	{
		for(int ir = 0; ir < nrxx; ++ir)
		{
			tmpchi1[ir]*=vr[ir];
		}
	}

	//------------------------------------
	// (3) the nonlocal pseudopotential.
	//------------------------------------
}

void Stochastic_Hchi::orthogonal_to_psi()
{
	TITLE("Stochastic_Hchi","orthogonal_to_psi0");


	return;
}
