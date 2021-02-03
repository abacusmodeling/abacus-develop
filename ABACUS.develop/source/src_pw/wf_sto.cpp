#include "wf_sto.h"
#include "global.h"

WF_Stochastic::WF_Stochastic()
{
    evc  = new ComplexMatrix[1];
}

WF_Stochastic::~WF_Stochastic()
{ 
    delete[] evc;
}


// update the function to generate random wave functions as stochastic wave functions.

void WF_Stochastic::random(ComplexMatrix &psi,const int iw_start,const int iw_end,const int ik)const
{
/*
    assert(iw_start >= 0);
    assert(psi.nr >= iw_end);
    const int ng = kv.ngk[ik];
    for (int iw = iw_start ;iw < iw_end;iw++)
    {
        for (int ig = 0;ig < ng;ig++)
        {
            const double rr = std::rand();
            const double arg= TWO_PI * std::rand();
            Vector3<double> v3 = kv.kvec_c[ik] + pw.gcar[this->igk(ik, ig)];
            psi(iw,ig) = complex<double>(rr * cos(arg), rr * sin(arg)) / (v3 * v3 + 1.0);
        }
        if(NPOL==2)for (int ig = wf.npwx;ig < wf.npwx + ng;ig++)
        {
            const double rr = std::rand();
            const double arg= TWO_PI * std::rand();
            Vector3<double> v3 = kv.kvec_c[ik] + pw.gcar[this->igk(ik, ig-wf.npwx)];
            psi(iw,ig) = complex<double>(rr * cos(arg), rr * sin(arg)) / (v3 * v3 + 1.0);
        }
    }
*/
    return;
}

