#ifndef USE_FFT_H
#define USE_FFT_H

#include "tools.h"

class Use_FFT
{
	public:

	Use_FFT();
	~Use_FFT();

	complex<double> *porter;

	void allocate(void);

	//---------------------------------------------------------------------
	
	// From G space to real space. ComplexMatrix
	void ToRealSpace(const int &is, const ComplexMatrix &vg, double *vr);
	void ToRealSpace_psi(const int &ik, const complex<double> *psig, complex<double> *psir);
	void ToRealSpace_psi(const int &ik, const int &ib, const ComplexMatrix &evc, complex<double> *psir);
	
	// From G space to real space. charge/MLWF
	void ToRealSpace(const complex<double> *vg, double *vr);

	// From G space to real space. wave functions.
	void ToRealSpace(const complex<double> *vg, complex<double> *vr);

	void ToRealSpace(const int &is, const ComplexMatrix &vg, matrix &v);

	//---------------------------------------------------------------------

	// From real space to G space.
	void ToReciSpace(const double* vr, complex<double> *vg);
	

	//---------------------------------------------------------------------
	
	void RoundTrip(
	    const complex<double> *psi,
		const double *vr,
		const int *_index,
		complex<double> *psic);

	private:


};

#endif
