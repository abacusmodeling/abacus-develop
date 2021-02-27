#include "tools.h"
#include "global.h"
#include "sto_iter.h"
#include "occupy.h" 

double Stochastic_Iter:: mu;
double Stochastic_Iter:: mu0;
double Stochastic_Iter:: Emin;
double Stochastic_Iter:: Emax;

Stochastic_Iter::Stochastic_Iter()
{
    mu0 = 0;
    spolyv = new double [1];
}

Stochastic_Iter::~Stochastic_Iter()
{
}
void Stochastic_Iter:: init()
{
    nchip = STO_WF.nchip;
    //wait for init
    targetne = ucell.nelec;
    stoche.init();
    stohchi.init();
    stohchi.get_GRA_index();
    delete [] spolyv;
    int norder = stoche.norder;
    spolyv = new double [norder];

}

void Stochastic_Iter:: itermu()
{
    

    int nkk=1;// We temporarily use gamma k point.
    //orthogonal part
    if(NBANDS > 0)
    {
        for(int ik = 0; ik < nkk; ++ik)
        {
            for(int ichi = 0; ichi < nchip; ++ichi)
            {
                complex<double> * p0 = &STO_WF.chi0[ik](ichi,0);
                complex<double> * pchi = &STO_WF.chiortho[ik](ichi,0);
                stohchi.orthogonal_to_psi_real(p0,pchi,ik);
            }
        }
    }

    //if(NBANDS > 0)
    //{
    //    Emin = wf.ekb[0][0];}
    //Emax = pw.ecutwfc * 20-250;
    Emax = 750;
    Emin = -10;
    Stochastic_hchi:: Emin = this->Emin;
    Stochastic_hchi:: Emax = this->Emax;
     
    sumpolyval();
    mu = mu0 - 2;
    double ne1 = calne();
    double mu1 = mu;
    mu = mu0 + 2;
    double ne2 = calne();
    double mu2 = mu;
    th_ne = 1e-6;
    double Dne = th_ne + 1;
    double ne3;
    double mu3;
    
    mu = 0;
    ne1 = calne();
    //test the domin of mu
    /*for(mu = -5; mu<5;mu+=0.2)
    {
        ne1 = calne();
        cout<<"mu: "<<mu<<" ; ne: "<<ne1<<endl;
    }*/
    while(ne1 > targetne)
    {
        mu1 -= 2;
        mu = mu1;
        ne1 = calne();
        cout<<"Reset mu1 from "<<mu1+2<<" to "<<mu1<<endl;
    }
    while(ne2 < targetne)
    {
        mu2 += 2;
        mu = mu2;
        ne2 = calne();
        cout<<"Reset mu2 form "<<mu2-2<<" to "<<mu2<<endl;
    }

    int count = 0;
    while(Dne > th_ne)
    {
        mu3 = (mu2 + mu1) / 2;
        mu = mu3;
        ne3 = calne();
        if(ne3 < targetne)
        {
            ne1 = ne3;
            mu1 = mu3;
        }
        else if(ne3 > targetne)
        {
            ne2 = ne3;
            mu2 = mu3;
        }
        Dne = abs(targetne - ne3);
        //cout<<setw(20)<<"targetne"<<setw(20)<<"ne"<<setw(20)<<"mu"<<setw(20)<<"Dne"<<endl;
        //cout<<setw(20)<<targetne<<setw(20)<<ne3<<setw(20)<<mu3<<setw(20)<<Dne<<endl;
        count++;
        if(count > 100)
        {
            cout<<"Fermi energy cannot be converged."<<endl;
        }
    }
    cout<<"Converge fermi energy = "<<mu<<" Ry in "<<count<<" steps."<<endl;

    mu = mu0 = mu3;
    
    //Set wf.wg 
    if(NBANDS > 0)
    {
        for(int ikk = 0; ikk < nkk; ++ikk)
        {
            double *en=wf.ekb[ikk];
            for(int iksb = 0; iksb < NBANDS; ++iksb)
            {
                wf.wg(ikk,iksb) = fd(en[iksb])*kv.wk[ikk];
            }
        }
    }
}

void Stochastic_Iter:: sumpolyval()
{
    int norder = stoche.norder;
    //wait for init
    int nkk=1;
    
    int nrxx = stohchi.nrxx;
    ZEROS(spolyv, norder);
    
    complex<double> * pchi;
    for(int ik = 0; ik < nkk; ++ik)
    {
        for(int ichi = 0; ichi < nchip; ++ichi)
        {
            if(NBANDS > 0)
            {
                pchi = &STO_WF.chiortho[ik](ichi,0);
            }  
            else
            {
                pchi = &STO_WF.chi0[ik](ichi,0);
            }
            stoche.calpolyval(stohchi.hchi_real, nrxx, pchi);
            for(int ior = 0; ior < norder; ++ior)
            {
                spolyv[ior] += stoche.polyvalue[ior].real();
            }
        }
    }

    return;
}


double Stochastic_Iter::calne()
{  
    //wait for init
    int nkk = 1;

    stoche.calcoef(this->nfd);
    int norder = stoche.norder;
    double totne = 0;
    KS_ne = 0;
    for(int ikk = 0; ikk < nkk; ++ikk)
    {
        for(int ior = 0; ior < norder; ++ior)
        {
            totne += stoche.coef[ior] * spolyv[ior] * kv.wk[ikk];
        }
        //cout<<endl;
        if(NBANDS > 0)
        {
            double *en=wf.ekb[ikk];
            //number of electrons in KS orbitals
            for(int iksb = 0; iksb < NBANDS; ++iksb)
            {
                KS_ne += fd(en[iksb]) * kv.wk[ikk];
            }
        }
            
    }
    totne += KS_ne;
    
    return totne;
}

void Stochastic_Iter::sum_stoband()
{  
    //stoche.calcoef(this->nfd);
    stoche.calcoef(this->nroot_fd);

    int nkk=1;// We temporarily use gamma k point.
    int nrxx = stohchi.nrxx;
    
    complex<double> * out = new complex<double> [nrxx];
    complex<double> * hout = new complex<double> [nrxx];

    double dr3 = nrxx / ucell.omega;
    double Ebar = (Emin + Emax)/2;
	double DeltaE = (Emax - Emin)/2;
   
    double out2, tmpne;
    double sto_ne;
    double sto_ekband=0;
    cout<<"Eband 1 "<<en.eband<<endl;
    complex<double> * pchi;
    for(int ik = 0; ik < nkk; ++ik)
    {
        sto_ekband =0;
        double *rho = CHR.rho[ik];
        tmpne = 0;
        for(int ichi = 0; ichi < nchip; ++ichi)
        {
            if(NBANDS > 0)
            {
                pchi = &STO_WF.chiortho[ik](ichi,0);
            }  
            else
            {
                pchi = &STO_WF.chi0[ik](ichi,0);
            }   
            stoche.calresult(stohchi.hchi_real, nrxx, pchi, out);
            stohchi.hchi_real(out,hout);
            for(int ir = 0; ir < nrxx; ++ir)
            {
                out2 = norm(out[ir]);
                tmpne = out2 * kv.wk[ik];
                sto_ekband += real(conj(out[ir]) * hout[ir]) * DeltaE + Ebar * out2;
                rho[ir] += tmpne * dr3; 
                sto_ne += tmpne;
                //rho[ir] += real(conj(pchi[ir]) * out[ir]) * kv.wk[ik] * nrxx / ucell.omega;
            }
        }
        en.eband += sto_ekband * kv.wk[ik];
    }
    cout<<"Eband 2 "<<en.eband<<endl;
    cout<<"Renormalize rho from ne = "<<sto_ne+KS_ne<<" to targetne = "<<targetne<<endl;
    double factor = targetne / (sto_ne + KS_ne);
    for(int ik = 0; ik < nkk; ++ik)
    {
        for(int ir = 0; ir < nrxx ; ++ir)
        {
            CHR.rho[ik][ir] *= factor;
        }
    }
    


    /*double sum = 0;
	for(int ir = 0; ir < pw.nrxx; ++ir )
	{
		sum += CHR.rho[0][ir];
	}
	cout<<"Total nele: "<<sum*ucell.omega/pw.nrxx<<endl;*/
    
    delete [] out;
    delete [] hout;
    return;
}

double Stochastic_Iter:: root_fd(double e)
{
    double e_mu = (e - mu) / Occupy::gaussian_parameter ;
    if(e_mu > 200)
        return 0;
    else
        return 1 / sqrt(1 + exp(e_mu));
}

double Stochastic_Iter:: nroot_fd(double e)
{
    double Ebar = (Emin + Emax)/2;
	double DeltaE = (Emax - Emin)/2;
    double ne_mu = (e * DeltaE + Ebar - mu) / Occupy::gaussian_parameter ;
    if(ne_mu > 200)
        return 0;
    else
        return 1 / sqrt(1 + exp(ne_mu));
}

double Stochastic_Iter:: fd(double e)
{
    double e_mu = (e - mu) / Occupy::gaussian_parameter ;
    if(e_mu > 100)
        return 0;
    else
        return 1 / (1 + exp(e_mu));
}

double Stochastic_Iter:: nfd(double e)
{
    double Ebar = (Emin + Emax)/2;
	double DeltaE = (Emax - Emin)/2;
    double ne_mu = (e * DeltaE + Ebar - mu) / Occupy::gaussian_parameter ;
    if(ne_mu > 100)
        return 0;
    else
        return 1 / (1 + exp(ne_mu));
}
 void Stochastic_Iter:: test()
 {
     //=====================test============================
    /*
    complex<double> *in = new complex<double> [pw.nrxx];
    
    complex<double> *chig1 = new complex<double> [wf.npw];
    complex<double> *chig2 = new complex<double> [wf.npw];
    ZEROS(in,pw.nrxx);
    ZEROS(in2,pw.nrxx);*/

    //---------------------------------------------------
    /*//test hermit property of  hchi matrix
    Emin = -1;
    Emax = 1;
    Stochastic_hchi:: Emin = this -> Emin;
    Stochastic_hchi:: Emax = this -> Emax;
    complex<double> *out = new complex<double> [pw.nrxx];
    complex<double> *in2 = new complex<double> [pw.nrxx];
    cout<<"------------------------------------"<<endl;
    complex<double> cij,cji;
    double dc;
    for(int i = 0 ; i < 300 ; ++i)
    {
        if( i % 10  == 0)
            cout<<"We are testing "<<i+1 <<" rows..."<<endl;
        for(int j = i+1; j < 300 ; ++j)
        {
            in2[j] = 1;
            stohchi.hchi_real(in2,out);
            cij = out[i];
            in2[j] = 0;
            in2[i] = 1;
            stohchi.hchi_real(in2,out);
            cji = out[j];
            in2[i] = 0;
            dc = abs(conj(cij)-cji);
            if(dc > 1e-6)
            {
                cout<<"(i,j) = ("<<i+1<<" , "<<j+1<<") ; cij = "<<cij<<" ; cji = "<<cji<<endl;
            }

        }
    }
    cout<<"------------------------------------"<<endl;
    delete[] out;
    delete[] in2;*/
    //---------------------------------------------------
    
    //-------------------------------------------------------
    /*//compare hchi and h_psi
    Emin = -1;
    Emax = 1;
    Stochastic_hchi:: Emin = this -> Emin;
    Stochastic_hchi:: Emax = this -> Emax;
    int * GRA_index = stohchi.GRA_index;
    ZEROS(in,pw.nrxx);
    complex<double> *kswf = &wf.evc[0](0,0);
    for ( int ig = 0 ; ig< wf.npw; ++ig)
    {
        in[GRA_index[ig]] = kswf [ig];
    }
    fftw_plan pp=fftw_plan_dft_3d(pw.nx,pw.ny,pw.nz,(fftw_complex *)in,(fftw_complex *)in2, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(pp);
    stohchi.hchi_real(in2,out);
    fftw_plan pp2=fftw_plan_dft_3d(pw.nx,pw.ny,pw.nz,(fftw_complex *)out,(fftw_complex *)out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(pp2);
    for(int i  = 0; i <wf.npw;++i)
    {
        chig1[i] = out[GRA_index[i]] / pw.nrxx; 
    }

    hm.hpw.h_psi( kswf , chig2);
    for(int i = 0; i<wf.npw;++i)
    {
        if(i % 100 == 0)
        cout<<kswf[i]<<" "<<in[GRA_index[i]]<<" "<<chig1[i]<<" "<<chig2[i]<<endl;
    }
    cout<<endl;*/
    //---------------------------------------------------
    
    //---------------------------------------------------
    //compare eigen energy
    /*
    int * GRA_index = stohchi.GRA_index;
    complex<double> *chigout = new complex<double> [wf.npw];
    complex<double> *wave = new complex<double> [pw.nrxx];
    complex<double> *waveout = new complex<double> [pw.nrxx];
    Emax = 1000;
    Emin = 0;
    Stochastic_hchi:: Emin = this->Emin;
    Stochastic_hchi:: Emax = this->Emax;
    double Ebar = (Emax + Emin)/2;
    double DeltaE = (Emax - Emin)/2;
    fftw_plan pp=fftw_plan_dft_3d(pw.nx,pw.ny,pw.nz,(fftw_complex *)wave,(fftw_complex *)wave, FFTW_BACKWARD, FFTW_ESTIMATE);
    for(int ib = 0 ; ib < NBANDS ; ++ib)
    {
        complex<double> *kswf = &wf.evc[0](ib,0);
        hm.hpw.h_psi( kswf , chigout);
        double energy = 0;
        double norm1 =0;
        ZEROS(wave,pw.nrxx);
        for(int ig = 0 ; ig < wf.npw ; ++ig)
        {
            energy += real(conj(kswf[ig]) * chigout[ig]);
            norm1 += norm (kswf[ig]);
            wave[GRA_index[ig]] = kswf[ig];
        }
        fftw_execute(pp);
        stohchi.hchi_real(wave,waveout);
        double energy2 = 0;
        double norm2 =0;
        for(int ir = 0 ; ir < pw.nrxx ; ++ir)
        {
            energy2 += real(conj(wave[ir]) * waveout[ir]) * DeltaE + Ebar * norm(wave[ir]);
            norm2 += norm(wave[ir]);
        }


        cout<<"BAND "<<ib+1<<" "<<energy/norm1<<"  "<<energy2/norm2<<endl;
    }
    delete []chigout;
    delete []wave;
    delete []waveout;*/
    //---------------------------------------------------

    //---------------------------------------------------
    //test ne
    /*
    int * GRA_index = stohchi.GRA_index;
    complex<double> *wave = new complex<double> [pw.nrxx];
    complex<double> *waveout = new complex<double> [pw.nrxx];
    Emax = 750;
    Emin = -100;
    Stochastic_hchi:: Emin = this->Emin;
    Stochastic_hchi:: Emax = this->Emax;
    mu = en.ef;
    stoche.calcoef(this->nfd);
    fftw_plan pp=fftw_plan_dft_3d(pw.nx,pw.ny,pw.nz,(fftw_complex *)wave,(fftw_complex *)wave, FFTW_BACKWARD, FFTW_ESTIMATE);
    for(int ib = 0 ; ib < NBANDS ; ++ib)
    {
        ZEROS(wave,pw.nrxx);
        complex<double> *kswf = &wf.evc[0](ib,0);
        for(int ig = 0 ; ig < wf.npw ; ++ig)
        {
            wave[GRA_index[ig]] = kswf[ig];
        }
        fftw_execute(pp);
        double ne =0;
        double norm1 = 0;
        stoche.calresult(stohchi.hchi_real, pw.nrxx, wave, waveout);
        for(int ir = 0 ; ir < pw.nrxx ; ++ir)
        {   
            ne += real(conj(wave[ir]) * waveout[ir]);
            norm1 += norm(wave[ir]);
        }
        cout<<"Ne of Band "<<ib+1<<" is "<<ne/norm1 * kv.wk[0]<<endl;
    }
    delete []wave;
    delete []waveout;*/


    //=====================================================
 }
