#include "global.h"
#include "sto_iter.h"
#include "occupy.h" 

double Stochastic_Iter:: mu;
double Stochastic_Iter:: Emin;
double Stochastic_Iter:: Emax;

Stochastic_Iter::Stochastic_Iter()
{
    mu = 0;
    spolyv = new double [1];
}

Stochastic_Iter::~Stochastic_Iter()
{
}
void Stochastic_Iter:: init()
{
    nchip = STO_WF.nchip;
    //wait for init
    targetne = 0;
    stoche.init();
    stohchi.init();
    delete [] spolyv;
    int norder = stoche.norder;
    spolyv = new double [norder];
    ZEROS(spolyv,norder);

}

void Stochastic_Iter:: itermu()
{
    //orthogonal part
    int nkk=1;// We temporarily use gamma k point.
    for(int ik = 0; ik < nkk; ++ik)
    {
        for(int ichi = 0; ichi < nchip; ++ichi)
        {
            complex<double> * p0 = &STO_WF.chi0[ik](ichi,0);
            complex<double> * pg = &STO_WF.chig[ik](ichi,0);
            stohchi.orthogonal_to_psi(p0,pg);
        }
    }

    sumpolyval();
    double dnedmu = caldnedmu();
    double ne1 = calne();
    double mu1 = mu;
    double mu2 = (targetne - ne1) / dnedmu;
    double Dne=abs(targetne - ne1);
    double ne2;
    double mu3;
    while(Dne > th_ne)
    {
        mu = mu2;
        ne2 = calne();
        mu3 = mu2;
        mu2 = (targetne - ne1) * (mu2 - mu1) / (ne2 - ne1) + mu1;
        mu1 = mu3;
        ne1 = ne2;
        Dne = abs(targetne - ne2);
    }

    //wait for init 
    double *rho;
    mu = mu1;
    calrho(rho);
}

void Stochastic_Iter:: sumpolyval()
{
    int norder = stoche.norder;
    //wait for init
    int nkk;
    int nrxx = stohchi.nrxx;

    for(int ik = 0; ik < nkk; ++ik)
    {
        for(int ichi = 0; ichi < nchip; ++ichi)
        {
            complex<double> * pg = &STO_WF.chig[ik](ichi,0);
            stoche.calpolyval(stohchi.hchi, nrxx, pg);
            for(int ior = 0; ior < norder; ++ior)
            {
                spolyv[ior] += stoche.polyvalue[ior].real();
            }
        }
    }
    return;
}

double Stochastic_Iter:: caldnedmu()
{
    stoche.calcoef(this->dfddmu);
    int norder = stoche.norder;
    double dnedmu = 0;
    for(int ior = 0; ior < norder; ++ior)
    {
        dnedmu += stoche.coef[ior] * spolyv[ior];
    }

    //wait for init
    double *en;

    //number of electrons in KS orbitals
    for(int iksb = 0; iksb < NBANDS; ++iksb)
    {
        dnedmu += fd(en[iksb]);
    }

    dnedmu *= 2;
    return dnedmu;
}

double Stochastic_Iter::calne()
{  
    stoche.calcoef(this->nfd);
    double ne = 0;
    int norder = stoche.norder;
    for(int ior = 0; ior < norder; ++ior)
    {
        ne += stoche.coef[ior] * spolyv[ior];
    }
    

    //wait for init
    double *en;

    //number of electrons in KS orbitals
    for(int iksb = 0; iksb < NBANDS; ++iksb)
    {
        ne += fd(en[iksb]);
    }

    ne *= 2;
    return ne;
}

void Stochastic_Iter::calrho( double * rho)
{  
    stoche.calcoef(this->nroot_fd);
    int nkk=1;// We temporarily use gamma k point.
    int nrxx = stohchi.nrxx;
    double ne = 0;
    complex<double> * out = new complex<double> [nrxx];
    for(int ik = 0; ik < nkk; ++ik)
    {
        for(int ichi = 0; ichi < nchip; ++ichi)
        {
            complex<double> * pg = &STO_WF.chig[ik](ichi,0);
            stoche.calresult(Stochastic_hchi::hchi, nrxx, pg, out);
            for(int ir = 0; ir < nrxx; ++ir)
            {
                rho[ir] += norm(out[ir]);
            }
        }
    }

    //wait for init
    double *rhoks;

    //number of electrons in KS orbitals
    for(int ir = 0 ; ir < nrxx; ++ir)
    {
        rho[ir] += rhoks [ir];
    }
    delete [] out;
    return;
}

double Stochastic_Iter:: dfddmu(double e)
{
    double expc = exp((e - mu) / (Occupy::gaussian_parameter * Ry_to_eV));
    return expc / pow(1 + expc, 2) / (Occupy::gaussian_parameter * Ry_to_eV);
}

double Stochastic_Iter:: ndfddmu(double e)
{
    double Ebar = (Emin + Emax)/2;
	double DeltaE = (Emax - Emin)/2;
    double expc = exp((e * DeltaE + Ebar - mu) / (Occupy::gaussian_parameter * Ry_to_eV));
    return expc / pow(1 + expc, 2) / (Occupy::gaussian_parameter * Ry_to_eV);
}

double Stochastic_Iter:: root_fd(double e)
{
    return sqrt(1 / (1 + exp((e - mu) / (Occupy::gaussian_parameter * Ry_to_eV))));
}

double Stochastic_Iter:: nroot_fd(double e)
{
    double Ebar = (Emin + Emax)/2;
	double DeltaE = (Emax - Emin)/2;
    return sqrt(1 / (1 + exp((e * DeltaE + Ebar - mu) / (Occupy::gaussian_parameter * Ry_to_eV))));
}

double Stochastic_Iter:: fd(double e)
{
    return 1 / (1 + exp((e - mu) / (Occupy::gaussian_parameter * Ry_to_eV)));
}

double Stochastic_Iter:: nfd(double e)
{
    double Ebar = (Emin + Emax)/2;
	double DeltaE = (Emax - Emin)/2;
    return 1 / (1 + exp((e * DeltaE + Ebar - mu) / (Occupy::gaussian_parameter * Ry_to_eV)));
}
