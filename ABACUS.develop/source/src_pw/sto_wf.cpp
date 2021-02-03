#include "sto_wf.h"
#include "global.h"

Stochastic_WF::Stochastic_WF()
{
    chi  = new ComplexMatrix[1];
}

Stochastic_WF::~Stochastic_WF()
{ 
    delete[] chi;
}

void Stochastic_WF::allocate_chi()
{


}


// update the function to generate random wave functions as stochastic wave functions.
void Stochastic_WF::random(ComplexMatrix &chi)
{
	int ik=0; // only supports gamma-only case now

    const int ng = kv.ngk[ik];

    for (int iw = 0 ;iw < this->nchi;iw++)
    {
        for (int ig = 0;ig < ng;ig++)
        {
  //          const double rr = std::rand();
  //          const double arg= TWO_PI * std::rand();
  //          Vector3<double> v3 = kv.kvec_c[ik] + pw.gcar[this->igk(ik, ig)];
  //          psi(iw,ig) = complex<double>(rr * cos(arg), rr * sin(arg)) / (v3 * v3 + 1.0);
        }
    }
    return;
}

