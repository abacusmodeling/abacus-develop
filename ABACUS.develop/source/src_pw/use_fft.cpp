#include "use_fft.h"
#include "global.h"

Use_FFT::Use_FFT()
{
	porter = new complex<double>[pw.nrxx];	
}

Use_FFT::~Use_FFT()
{
	delete[] porter;	
}

void Use_FFT::allocate(void)
{
    TITLE("Use_FFT","allocate");

    delete[] porter;
    porter = new complex<double>[pw.nrxx];
    Memory::record("Use_FFT","porter",pw.nrxx,"complexmatrix");

    return;
}

void Use_FFT::RoundTrip(
    const complex<double> *psi,
    const double *vr,
    const int *fft_index,
    complex<double> *psic)
{
	// (1) set value
    for (int ig=0; ig< wf.npw; ig++)
    {
        psic[ fft_index[ig]  ] = psi[ig];
    }

	// (2) fft to real space and doing things.
    pw.FFT_wfc.FFT3D( psic, 1);
    for (int ir=0; ir< pw.nrxx; ir++)
    {
        psic[ir] *= vr[ir];
    }

	// (3) fft back to G space.
    pw.FFT_wfc.FFT3D( psic, -1);
    return;
}

void Use_FFT::ToRealSpace(const int &is, const ComplexMatrix &vg, double *vr)
{
	// (1) set value
    ZEROS( porter, pw.nrxx );
    for (int ig=0; ig<pw.ngmc; ig++)
    {
        porter[ pw.ig2fftc[ig] ] = vg(is, ig);
    }

	// (2) fft and get value
    pw.FFT_chg.FFT3D(porter, 1);
    for (int ir = 0; ir < pw.nrxx; ir++)
    {
        vr[ir] = porter[ir].real();
    }
    return;
}

void Use_FFT::ToRealSpace_psi(const int &ik, const complex<double> *psig, complex<double> *psir)
{
	ZEROS(psir, pw.nrxx);
	for(int ig=0; ig<wf.npw; ig++)
	{
		psir[ pw.ig2fftc[ wf.igk(ik,ig) ] ] = psig[ig];
	}
	// 1: to real space.
	pw.FFT_wfc.FFT3D(psir, 1);

	return;
}


void Use_FFT::ToRealSpace_psi(const int &ik, const int &ib, const ComplexMatrix &evc, complex<double> *psir)
{
	// (1) set value
    ZEROS( psir, pw.nrxx );
    for (int ig=0; ig<wf.npw; ig++)
    {
        psir[ pw.ig2fftc[ wf.igk(ik,ig) ] ] = evc(ib, ig);
    }

	// (2) fft and get value
    pw.FFT_wfc.FFT3D(psir, 1);
    return;
}



void Use_FFT::ToRealSpace(const int &is, const ComplexMatrix &vg, matrix &vr)
{
	// (1) set value
    ZEROS( porter, pw.nrxx);
    for (int ig=0; ig<pw.ngmc; ig++)
    {
        porter [pw.ig2fftc[ig]] = vg(is,ig);
    }

	// (2) fft and get value
    pw.FFT_chg.FFT3D(porter, 1);
    for (int ir = 0; ir < pw.nrxx; ir++)
    {
        vr(is,ir) = porter[ir].real();
    }
    return;
}


// Fourer transform of vg, 
// then put vg into vr.
void Use_FFT::ToRealSpace(const complex<double> *vg, double *vr)
{
    ZEROS( porter, pw.nrxx);
    for (int ig=0; ig<pw.ngmc; ig++)
    {
        porter[pw.ig2fftc[ig]] = vg[ig];
    }
    pw.FFT_chg.FFT3D(porter, 1);
    for (int ir = 0; ir < pw.nrxx; ir++)
    {
        vr[ir] = porter[ir].real();
    }
    return;
}

void Use_FFT::ToReciSpace(const double* vr, complex<double> *vg)
{
	ZEROS( porter, pw.nrxx);
	for (int ir=0; ir<pw.nrxx; ir++)
	{
		porter[ir] = complex<double>(vr[ir], 0.0);
	}
	pw.FFT_chg.FFT3D( porter, -1);
	for (int ig=0; ig<pw.ngmc; ig++)
	{
		vg[ig] = porter[ pw.ig2fftc[ig] ];
	}
	return;
}

