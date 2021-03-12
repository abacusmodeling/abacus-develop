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
    delete[] spolyv;
}
void Stochastic_Iter:: init()
{
    nchip = STO_WF.nchip;
    stotype = STO_WF.stotype;
    //wait for init
    targetne = ucell.nelec;
    stoche.init();
    stohchi.init();
    stohchi.get_GRA_index();
    delete [] spolyv;
    int norder = stoche.norder;
    spolyv = new double [norder];

}

void Stochastic_Iter:: orthog()
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
                if(stotype == "pw")
                {
                    stohchi.orthogonal_to_psi_reciprocal(p0,pchi,ik);
                }
                else
                {
                    stohchi.orthogonal_to_psi_real(p0,pchi,ik);
                }   
            }
        }
    }
}

void Stochastic_Iter:: checkemm(int &iter)
{
    if(iter == 1)
    {
        Stochastic_hchi:: Emin = STO_WF.emin_sto;
        Stochastic_hchi:: Emax = STO_WF.emax_sto;
    }
    else if(iter > 5)
        return;
    
    int norder = stoche.norder;
    int ndim;
    if(stotype == "pw")
        ndim = wf.npw;
    else
        ndim = pw.nrxx;
        
    //wait for init
    
    complex<double> * pchi;
    int ntest = 2;
    if (nchip < ntest) ntest = nchip;

    for(int ichi = 0; ichi < ntest; ++ichi)
    {
        if(NBANDS > 0)
        {
            pchi = &STO_WF.chiortho[0](ichi,0);
        }  
        else
        {
            pchi = &STO_WF.chi0[0](ichi,0);
        }

        while(1)
        {
            bool converge;
            if(stotype == "pw")
                converge = stoche.checkconverge(stohchi.hchi_reciprocal, ndim, pchi, Stochastic_hchi:: Emax, Stochastic_hchi:: Emin, 5);
            else
                converge = stoche.checkconverge(stohchi.hchi_real, ndim, pchi, Stochastic_hchi:: Emax, Stochastic_hchi:: Emin, 5);
            if(!converge)
                cout<<"New Emax "<<Stochastic_hchi:: Emax<<" ; new Emin "<<Stochastic_hchi:: Emin<<endl;
            else
                break;
        }
    }

    Emax = Stochastic_hchi:: Emax;
    Emin = Stochastic_hchi:: Emin;
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE, &Emax, 1, MPI_DOUBLE, MPI_MAX , MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &Emin, 1, MPI_DOUBLE, MPI_MIN , MPI_COMM_WORLD);
#endif
    Stochastic_hchi:: Emin = Emin;
    Stochastic_hchi:: Emax = Emax;
    
}

void Stochastic_Iter:: itermu( int &iter) 
{
    int nkk=1;// We temporarily use gamma k point.
    double dmu;
    if(iter == 1)
    {
        dmu = 2;
        th_ne = 0.1 * DRHO2 * ucell.nelec;
        cout<<"th_ne "<<th_ne<<endl;
    }
    else
    {
        dmu = 0.1;
        th_ne = DRHO2 * 1e-2 * ucell.nelec;
    }
    
    sumpolyval();
    mu = mu0 - dmu;
    double ne1 = calne();
    double mu1 = mu;

    mu = mu0 + dmu;
    double ne2 = calne();
    double mu2 = mu;
    double Dne = th_ne + 1;
    double ne3;
    double mu3;
    
    //test the domin of mu
    /*for(mu = -5; mu<5;mu+=0.2)
    {
        ne3 = calne();
        cout<<"mu: "<<mu<<" ; ne: "<<ne3<<endl;
    }
    exit(0);*/
    while(ne1 > targetne)
    {
        ne2 = ne1;
        mu2 = mu1;
        mu1 -= dmu;
        mu = mu1;
        ne1 = calne();
        cout<<"Reset mu1 from "<<mu1+dmu<<" to "<<mu1<<endl;
        dmu *= 2;
    }
    while(ne2 < targetne)
    {
        ne1 = ne2;
        mu1 = mu2;
        mu2 += dmu;
        mu = mu2;
        ne2 = calne();
        cout<<"Reset mu2 form "<<mu2-dmu<<" to "<<mu2<<endl;
        dmu *= 2;
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
        if(count > 60)
        {
            cout<<"Fermi energy cannot be converged. Set THNE to "<<th_ne<<endl;
            th_ne *= 1e4;
        }
    }
    cout<<"Converge fermi energy = "<<mu<<" Ry in "<<count<<" steps."<<endl;

    en.ef = mu = mu0 = mu3;
    
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

    return;
}

void Stochastic_Iter:: sumpolyval()
{
    int norder = stoche.norder;
    //wait for init
    int nkk=1;
    
    int nrxx = stohchi.nrxx;
    int npw = wf.npw;
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
            
            if(stotype == "pw")
            {
                stoche.calpolyval(stohchi.hchi_reciprocal, npw, pchi);
            }
            else
            {
                stoche.calpolyval(stohchi.hchi_real, nrxx, pchi);
            }
            //LapackConnector::axpy(norder,1,stoche.polyvalue,1,spolyv,1);
            for(int ior = 0; ior < norder; ++ior)
            {
                spolyv[ior] += stoche.polyvalue[ior];
            }
        }
    }
    return;
}


double Stochastic_Iter::calne()
{  
    //wait for init
    int nkk = 1;

    stoche.calcoef(nfd);
    int norder = stoche.norder;
    double totne = 0;
    double sto_ne = 0;
    KS_ne = 0;
    for(int ikk = 0; ikk < nkk; ++ikk)
    {
        double stok_ne = LapackConnector::dot(norder,stoche.coef,1,spolyv,1);
        //double stok_ne = 0;
        //for(int ior = 0; ior < norder; ++ior)
        //{
        //    stok_ne += stoche.coef[ior] * spolyv[ior];
        //    //cout<<stoche.coef[ior]<<" "; 
        //}
        //cout<<endl;
        //cout<<"last term "<<stoche.coef[norder-1] * spolyv[norder-1]<<endl;
        if(NBANDS > 0)
        {
            double *en=wf.ekb[ikk];
            //number of electrons in KS orbitals
            for(int iksb = 0; iksb < NBANDS; ++iksb)
            {
                KS_ne += fd(en[iksb]) * kv.wk[ikk];
            }
        }
        sto_ne += stok_ne * kv.wk[ikk]; 
    }

#ifdef __MPI
        MPI_Allreduce(MPI_IN_PLACE, &sto_ne, 1, MPI_DOUBLE, MPI_SUM , MPI_COMM_WORLD);
#endif

    totne = KS_ne + sto_ne;
    return totne;
}

void Stochastic_Iter::sum_stoband()
{  
    int nkk=1;// We temporarily use gamma k point.
    int nrxx = stohchi.nrxx;
    int npw = wf.npw;
    int norder = stoche.norder;

    //cal demet
    stoche.calcoef(this->nfdlnfd);
    double stodemet = 0;

    for(int ikk = 0; ikk < nkk; ++ikk)
    {
        double stok_demet = LapackConnector::dot(norder,stoche.coef,1,spolyv,1);
        //double stok_demet=0;
        //for(int ior = 0; ior < norder; ++ior)
        //{
        //    stok_demet += stoche.coef[ior] * spolyv[ior];
        //}
        //cout<<"last term "<<stoche.coef[norder-1] * spolyv[norder-1]<<endl;
        //cout<<endl;
        if(NBANDS > 0)
        {
            double *enb=wf.ekb[ikk];
            //number of electrons in KS orbitals
            for(int iksb = 0; iksb < NBANDS; ++iksb)
            {
                en.demet += fdlnfd(enb[iksb]) * kv.wk[ikk];
            }
        }
        stodemet += stok_demet * kv.wk[ikk];
    }

    //cal eband & rho
    stoche.calcoef(this->nroot_fd);
    //stoche.calcoef(this->nfd);
    
    complex<double> * out, *hout;
    double *sto_rho = new double [nrxx];
    complex<double> *sto_rhog;
    if(stotype == "pw")
    {
        out = new complex<double> [npw];
        hout = new complex<double> [npw];
        sto_rhog = new complex<double> [npw];
    }
    else
    {
        out = new complex<double> [nrxx];
        hout = new complex<double> [nrxx];
    }

    double dr3 = ucell.omega / nrxx;
    double Ebar = (Emin + Emax)/2;
	double DeltaE = (Emax - Emin)/2;
   
    double tmprho, tmpne;
    complex<double> outtem;
    double sto_ne = 0;
    double sto_eband = 0;
    ZEROS(sto_rho, nrxx);

    complex<double> * pchi;
    complex<double>* porter = UFFT.porter;
    int* GRA_index = stohchi.GRA_index;
    double out2;
    
    for(int ik = 0; ik < nkk; ++ik)
    {
        double stok_eband =0;
        if(stotype == "pw")
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
                
                stoche.calresult(stohchi.hchi_reciprocal, npw, pchi, out);
                
                stohchi.hchi_reciprocal(out,hout);
                ZEROS( porter, pw.nrxx );
                for(int ig = 0; ig < npw; ++ig)
                {
                    outtem = out[ig];
                    stok_eband += real(conj(outtem) * hout[ig]) * DeltaE + Ebar * norm(outtem);
                    porter[ GRA_index[ig] ] = outtem;
                }
                pw.FFT_wfc.FFT3D(UFFT.porter, 1);
                for(int ir = 0 ; ir < nrxx ; ++ir)
                {
                    tmprho = norm(porter[ir]) * kv.wk[ik] / ucell.omega;
                    sto_rho[ir] += tmprho;
                    sto_ne += tmprho * dr3;
                }

            }
            sto_eband += stok_eband * kv.wk[ik];
        }
        else
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
                stoche.calresult(stohchi.hchi_real, nrxx, pchi, out);
                stohchi.hchi_real(out,hout);
                for(int ir = 0; ir < nrxx; ++ir)
                {
                    out2 = norm(out[ir]);
                    //out2 = real(conj(pchi[ir]) * out[ir]);
                    tmpne = out2 * kv.wk[ik];
                    stok_eband += real(conj(out[ir]) * hout[ir]) * DeltaE + Ebar * out2;
                    sto_rho[ir] += tmpne; 
                    sto_ne += tmpne;
                } 
            }
            sto_eband += stok_eband * kv.wk[ik];
        }
           
    }

    double *f_storho = new double [nrxx];

#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE,&stodemet,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(sto_rho,f_storho,nrxx,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&sto_eband,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&sto_ne,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
    en.eband += sto_eband;
    en.demet += stodemet;
    en.demet *= Occupy::gaussian_parameter;

    
    cout.precision(12);
    cout<<"Normalize rho from ne = "<<sto_ne+KS_ne<<" to targetne = "<<targetne<<endl;


    double factor;
    if(stotype == "pw")
    {
        if(abs(sto_ne) > 1e-20)
            factor = (targetne - KS_ne) / sto_ne;
        else
            factor = 1;
    }
    else
    {
        if(abs(sto_ne) > 1e-20)
            factor = (targetne - KS_ne) / sto_ne / dr3;
        else
            factor = 1 / dr3;
    }
    
    
    
    for(int is = 0 ; is < 1; ++is)
    {
        //LapackConnector::axpy(nrxx, factor,sto_rho,1,CHR.rho[is],1);
        for(int ir = 0; ir < nrxx ; ++ir)
        {
#ifdef __MPI
            CHR.rho[is][ir] += f_storho[ir] * factor;
#else
            CHR.rho[is][ir] += sto_rho[ir] * factor;
#endif 
        }
    }



    if(stotype == "pw")
    {
        delete [] sto_rhog;
    }
    delete [] f_storho;
    delete [] sto_rho;
    delete [] out;
    delete [] hout;
    return;
}

double Stochastic_Iter:: root_fd(double e)
{
    double e_mu = (e - mu) / Occupy::gaussian_parameter ;
    if(e_mu > 72)
        return 0;
    else
        return 1 / sqrt(1 + exp(e_mu));
}

double Stochastic_Iter:: nroot_fd(double e)
{
    double Ebar = (Emin + Emax)/2;
	double DeltaE = (Emax - Emin)/2;
    double ne_mu = (e * DeltaE + Ebar - mu) / Occupy::gaussian_parameter ;
    if(ne_mu > 72)
        return 0;
    else
        return 1 / sqrt(1 + exp(ne_mu));
}

double Stochastic_Iter:: fd(double e)
{
    double e_mu = (e - mu) / Occupy::gaussian_parameter ;
    if(e_mu > 36)
        return 0;
    else
        return 1 / (1 + exp(e_mu));
}

double Stochastic_Iter:: nfd(double e)
{
    double Ebar = (Emin + Emax)/2;
	double DeltaE = (Emax - Emin)/2;
    double ne_mu = (e * DeltaE + Ebar - mu) / Occupy::gaussian_parameter ;
    if(ne_mu > 36)
        return 0;
    else
        return 1 / (1 + exp(ne_mu));
}

double Stochastic_Iter:: fdlnfd(double e)
{
    double e_mu = (e - mu) / Occupy::gaussian_parameter ;
    if(e_mu > 36)
        return 0;
    else if(e_mu < -36)
        return 0;
    else
    {
        double f = 1 / (1 + exp(e_mu));
        return (f * log(f) + (1.0-f) * log(1.0-f)); 
    }
}

double Stochastic_Iter:: nfdlnfd(double e)
{
    double Ebar = (Emin + Emax)/2;
	double DeltaE = (Emax - Emin)/2;
    double ne_mu = (e * DeltaE + Ebar - mu) / Occupy::gaussian_parameter ;
    if(ne_mu > 36)
        return 0;
    else if(ne_mu < -36)
        return 0;
    else
    {
        double f = 1 / (1 + exp(ne_mu));
        return f * log(f) + (1-f) * log(1-f); 
    }
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
    //compare hchi_reciprocal and h_psi
    /*Emin = -1;
    Emax = 1;
    Stochastic_hchi:: Emin = this -> Emin;
    Stochastic_hchi:: Emax = this -> Emax;
    complex<double> *chig1 = new complex<double> [wf.npw];
    complex<double> *chig2 = new complex<double> [wf.npw];
    complex<double> *kswf = &wf.evc[0](0,0);
    
    
    stohchi.hchi_reciprocal(kswf,chig1);
    
    hm.hpw.h_psi( kswf , chig2);
    for(int i = 0; i<wf.npw;++i)
    {
        if(i % 100 == 0)
            cout<<kswf[i]<<" "<<chig1[i]<<" "<<chig2[i]<<endl;
    }
    cout<<endl;*/
    //---------------------------------------------------

     //-------------------------------------------------------
    /*//compare hchi_real and h_psi
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
