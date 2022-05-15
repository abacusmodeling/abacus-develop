#include "global.h"
#include "sto_iter.h"
#include "occupy.h"
#include "diago_cg.h" 
#include "../module_base/tool_title.h"
#include "../module_base/tool_quit.h"
#include "../module_base/timer.h"
#include "../src_parallel/parallel_reduce.h"

double Stochastic_Iter:: mu=0;
double Stochastic_Iter:: mu0=0;
double Stochastic_Iter:: Emin=0;
double Stochastic_Iter:: Emax=0;
double Stochastic_Iter:: fwhm=0.03;
double Stochastic_Iter:: targ_e=0;

Stochastic_Iter::Stochastic_Iter()
{
    spolyv = new double [1];
    change = false;
}

Stochastic_Iter::~Stochastic_Iter()
{
    delete[] spolyv;
}

void Stochastic_Iter::init(const int dim, int* nchip_in)
{
    nchip = nchip_in;
    targetne = GlobalC::CHR.nelec;
    stoche.init( dim, INPUT.nche_sto );
    stohchi.init();
    delete [] spolyv;
    int norder = stoche.norder;
    spolyv = new double [norder];
    Emin = INPUT.emin_sto;
    Emax = INPUT.emax_sto;
}

void Stochastic_Iter::orthog(const int& ik, Stochastic_WF& stowf)
{
    ModuleBase::TITLE("Stochastic_Iter","orthog");
    //orthogonal part
    if(GlobalV::NBANDS > 0)
    {
	    const int nchipk=stowf.nchip[ik];
	    const int npw = GlobalC::wf.npw;
	    const int npwx = GlobalC::wf.npwx;
        std::complex<double> *wfgin = stowf.chi0[ik].c, *wfgout = stowf.chiortho[ik].c;
	    for(int ig = 0 ; ig < npwx * nchipk; ++ig)
	    {
		    wfgout[ig] = wfgin[ig];
	    }

	    //orthogonal part
	    complex<double> *sum = new complex<double> [GlobalV::NBANDS * nchipk];
	    char transC='C';
	    char transN='N';
    
	    //sum(b<NBANDS, a<nchi) = < psi_b | chi_a >
	    zgemm_(&transC, &transN, &GlobalV::NBANDS, &nchipk, &npw, &ModuleBase::ONE, 
                GlobalC::wf.evc[ik].c, &npwx, wfgout, &npwx, &ModuleBase::ZERO, sum, &GlobalV::NBANDS);
	    Parallel_Reduce::reduce_complex_double_pool(sum, GlobalV::NBANDS * nchipk);
    
	    //psi -= psi * sum
	    zgemm_(&transN, &transN, &npw, &nchipk, &GlobalV::NBANDS, &ModuleBase::NEG_ONE, 
                GlobalC::wf.evc[ik].c, &npwx, sum, &GlobalV::NBANDS, &ModuleBase::ONE, wfgout, &npwx);
	    delete[] sum;
    }
}

void Stochastic_Iter::checkemm(const int& ik, int &iter, Stochastic_WF& stowf)
{
    ModuleBase::TITLE("Stochastic_Iter","checkemm");
    if(iter > 5)
	{
        return;
	}
        
    int norder = stoche.norder;
    complex<double> * pchi;
    int ntest = 1;

    if (nchip[ik] < ntest) 
	{
		ntest = nchip[ik];
	}

    for(int ichi = 0; ichi < ntest; ++ichi)
    {
        if(GlobalV::NBANDS > 0)
        {
            pchi = &stowf.chiortho[ik](ichi,0);
        }  
        else
        {
            pchi = &stowf.chi0[ik](ichi,0);
        }
        while(1)
        {
            bool converge;
            converge = stoche.checkconverge(
				stohchi.hchi_reciprocal, 
				pchi, 
				Stochastic_hchi::Emax, 
				Stochastic_hchi::Emin, 
				5);

            if(!converge)
			{
                change = true;
			}
            else
			{
                break;
			}
        }
    }
    if(ik == GlobalC::kv.nks-1)
    {
        Emax = Stochastic_hchi:: Emax;
        Emin = Stochastic_hchi:: Emin;

#ifdef __MPI
        MPI_Allreduce(MPI_IN_PLACE, &Emax, 1, MPI_DOUBLE, MPI_MAX , MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &Emin, 1, MPI_DOUBLE, MPI_MIN , MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &change, 1, MPI_CHAR, MPI_LOR , MPI_COMM_WORLD);
#endif
        Stochastic_hchi:: Emin = Emin;
        Stochastic_hchi:: Emax = Emax;
        if(change)
	    {
	    	GlobalV::ofs_running<<"New Emax "<<Stochastic_hchi:: Emax<<" ; new Emin "<<Stochastic_hchi:: Emin<<std::endl;
	    }
        change = false;
    }
}

void Stochastic_Iter::itermu(int &iter) 
{
    ModuleBase::TITLE("Stochastic_Iter","itermu");
    ModuleBase::timer::tick("Stochastic_Iter","itermu");
    double dmu;
    if(iter == 1)
    {
        dmu = 2;
        th_ne = 0.1 * GlobalV::SCF_THR * GlobalC::CHR.nelec;
        // std::cout<<"th_ne "<<th_ne<<std::endl;
    }
    else
    {
        dmu = 0.1;
        th_ne = GlobalV::SCF_THR * 1e-2 * GlobalC::CHR.nelec;
    }
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
        std::cout<<"Reset mu1 from "<<mu1+dmu<<" to "<<mu1<<std::endl;
        dmu *= 2;
    }
    while(ne2 < targetne)
    {
        ne1 = ne2;
        mu1 = mu2;
        mu2 += dmu;
        mu = mu2;
        ne2 = calne();
        // cout<<"Reset mu2 from "<<mu2-dmu<<" to "<<mu2<<endl;
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
            std::cout<<"Fermi energy cannot be converged. Set THNE to "<<th_ne<<std::endl;
            th_ne *= 1e1;
            if(th_ne > 1e1) ModuleBase::WARNING_QUIT("Stochastic_Iter",
                                         "Cannot converge feimi energy. Please retry with different random number");
        }
    }
    GlobalV::ofs_running<<"Converge fermi energy = "<<mu<<" Ry in "<<count<<" steps."<<std::endl;
    //precision check
    double tmpre;
    tmpre = stoche.coef[stoche.norder-1] * spolyv[stoche.norder-1];
    MPI_Allreduce(MPI_IN_PLACE, &tmpre, 1, MPI_DOUBLE, MPI_SUM , MPI_COMM_WORLD);
    GlobalV::ofs_running<<"Chebyshev Precision: "<<abs(tmpre/targetne)*1e9<<"E-09"<<std::endl;
    if(tmpre/targetne > GlobalV::SCF_THR )
    {
        stringstream ss;
        ss<<tmpre/targetne;
        string fractxt,tartxt,itertxt;
        ss>>fractxt;
        ss.clear();
        ss<<GlobalV::SCF_THR;
        ss>>tartxt;
        ss.clear();
        ss<<iter+1;
        ss>>itertxt;
        string warningtxt = "Iter "+itertxt+": (Chebyshev error = "+fractxt+" > threshold = "+tartxt+" ) Please add more expansion terms for Chebychev expansion.";
        ModuleBase::WARNING("Stochastic_Chebychev", warningtxt);
    }


    GlobalC::en.ef = mu = mu0 = mu3;
    
    //Set wf.wg 
    if(GlobalV::NBANDS > 0)
    {
        for(int ikk = 0; ikk < GlobalC::kv.nks; ++ikk)
        {
            double *en = GlobalC::wf.ekb[ikk];
            for(int iksb = 0; iksb < GlobalV::NBANDS; ++iksb)
            {
                GlobalC::wf.wg(ikk,iksb) = fd(en[iksb])*GlobalC::kv.wk[ikk];
            }
        }
    }
    ModuleBase::timer::tick("Stochastic_Iter","itermu");
    return;
}

void Stochastic_Iter:: sumpolyval_k(const int& ik, Stochastic_WF& stowf)
{
    ModuleBase::TITLE("Stochastic_Iter","sumpolyval_k");
    ModuleBase::timer::tick("Stochastic_Iter","sumpolyval");
    int norder = stoche.norder;
    if(ik==0)   ModuleBase::GlobalFunc::ZEROS(spolyv, norder);
    std::complex<double> * pchi;
    if(GlobalV::NBANDS > 0)  pchi = stowf.chiortho[ik].c; 
    else            pchi = stowf.chi0[ik].c;
    stoche.calpolyval(stohchi.hchi_reciprocal, pchi, nchip[ik]);
    for(int i = 0 ; i < norder ; ++i)
    {
        spolyv[i] += stoche.polyvalue[i] * GlobalC::kv.wk[ik];
    }
    ModuleBase::timer::tick("Stochastic_Iter","sumpolyval");
    return;
}


double Stochastic_Iter::calne()
{  
    ModuleBase::timer::tick("Stochastic_Iter","calne");

    stoche.calcoef(nfd);
    int norder = stoche.norder;
    double totne = 0;
    KS_ne = 0;
    double sto_ne = BlasConnector::dot(norder,stoche.coef,1,spolyv,1);
    if(GlobalV::NBANDS > 0)
    {
        for(int ikk = 0; ikk < GlobalC::kv.nks; ++ikk)
        {
            double *en=GlobalC::wf.ekb[ikk];
            for(int iksb = 0; iksb < GlobalV::NBANDS; ++iksb)
            {
                KS_ne += fd(en[iksb]) * GlobalC::kv.wk[ikk];
            }
        }
    }
    KS_ne /= GlobalV::NPROC_IN_POOL;
	MPI_Allreduce(MPI_IN_PLACE, &KS_ne, 1, MPI_DOUBLE, MPI_SUM , STO_WORLD);

#ifdef __MPI
        MPI_Allreduce(MPI_IN_PLACE, &sto_ne, 1, MPI_DOUBLE, MPI_SUM , MPI_COMM_WORLD);
#endif

    totne = KS_ne + sto_ne;
    ModuleBase::timer::tick("Stochastic_Iter","calne");
    return totne;
}

void Stochastic_Iter::sum_stoband(Stochastic_WF& stowf)
{  
    ModuleBase::TITLE("Stochastic_Iter","sum_stoband");
    ModuleBase::timer::tick("Stochastic_Iter","sum_stoband");
    int nrxx = GlobalC::pw.nrxx;
    int npwx = GlobalC::wf.npwx;
    int norder = stoche.norder;

    //cal demet
    stoche.calcoef(this->nfdlnfd);
    double stodemet = BlasConnector::dot(norder,stoche.coef,1,spolyv,1);

    if(GlobalV::NBANDS > 0)
    {
        for(int ikk = 0; ikk < GlobalC::kv.nks; ++ikk)
        {
            double *enb=GlobalC::wf.ekb[ikk];
            //number of electrons in KS orbitals
            for(int iksb = 0; iksb < GlobalV::NBANDS; ++iksb)
            {
                GlobalC::en.demet += fdlnfd(enb[iksb]) * GlobalC::kv.wk[ikk];
            }
        }
    }
    GlobalC::en.demet /= GlobalV::NPROC_IN_POOL;
	MPI_Allreduce(MPI_IN_PLACE, &GlobalC::en.demet, 1, MPI_DOUBLE, MPI_SUM , STO_WORLD);

    //cal eband
    stoche.calcoef(this->nxfd);
    double sto_eband = BlasConnector::dot(norder,stoche.coef,1,spolyv,1);
        //double sto_eband=0;



    //cal rho
    stoche.calcoef(this->nroot_fd);
    
    complex<double> * out, *hout;
    double *sto_rho = new double [nrxx];
    //int npwall = npwx * nchip;

    double dr3 = GlobalC::ucell.omega / GlobalC::pw.ncxyz;
    double Ebar = (Emin + Emax)/2;
	double DeltaE = (Emax - Emin)/2;
   
    double tmprho, tmpne;
    complex<double> outtem;
    double sto_ne = 0;
    ModuleBase::GlobalFunc::ZEROS(sto_rho, nrxx);

    complex<double> * pchi;
    complex<double>* porter = GlobalC::UFFT.porter;
    double out2;

    double *ksrho;
    if(GlobalV::NBANDS > 0 && GlobalV::MY_STOGROUP==0 )
    {
        ksrho = new double [nrxx];
        ModuleBase::GlobalFunc::DCOPY(GlobalC::CHR.rho[0],ksrho,nrxx);
        ModuleBase::GlobalFunc::ZEROS(GlobalC::CHR.rho[0],nrxx);
    }
    
    for(int ik = 0; ik < GlobalC::kv.nks; ++ik)
    {
        //init k
        if(GlobalC::kv.nks > 1) GlobalC::hm.hpw.init_k(ik);
        stoche.ndmin = GlobalC::wf.npw;

        int npw = GlobalC::kv.ngk[ik];
        double stok_eband;
        out = stowf.shchi[ik].c;
        if(GlobalV::NBANDS > 0)
            pchi = stowf.chiortho[ik].c;
        else
            pchi = stowf.chi0[ik].c;
        
        stoche.calfinalvec(stohchi.hchi_reciprocal, pchi, out, nchip[ik]);
            //hout = new complex<double> [npwall];
            //stohchi.hchi_reciprocal(out,hout,nchip);
            //stok_eband = Diago_CG::ddot_real(npwall, out, hout,false) * DeltaE 
            //            +  Diago_CG::ddot_real(npwall, out, out,false) * Ebar;
            //sto_eband += stok_eband * kv.wk[ik];
        std::complex<double> *tmpout = out;
        for(int ichi = 0; ichi < nchip[ik] ; ++ichi)
        {
            ModuleBase::GlobalFunc::ZEROS( porter, GlobalC::pw.nrxx );
            for(int ig = 0; ig < npw; ++ig)
            {
                porter[ GlobalC::pw.ig2fftw[GlobalC::wf.igk(ik, ig)] ] = tmpout[ig];
            }
            GlobalC::pw.FFT_wfc.FFT3D(GlobalC::UFFT.porter, 1);
            for(int ir = 0 ; ir < nrxx ; ++ir)
            {
                GlobalC::CHR.rho[0][ir] += norm(porter[ir]) * GlobalC::kv.wk[ik];
            }
            tmpout+=npwx;
        }
            //delete [] hout;
    }
    //for(int ip = 0 ; ip < NPOOL ; ++ip)
    //{
    //    MPI_Barrier(MPI_COMM_WORLD);
    //    if(ip == MY_POOL)
    //    {
    //        if(RANK_IN_POOL==0)
    //        {
    //            cout.clear();
    //            for(int ik = 0; ik < kv.nks; ++ik)
    //            {
    //                complex<double> *tmpout = STO_WF.shchi[ik].c;
    //                for(int ichi = 0; ichi < nchip[ik] ; ++ichi)
    //                {
    //                    for(int ig = 0; ig < kv.ngk[ik]; ++ig)
    //                    {
    //                        if(ig%100==0) cout<<tmpout[ig]<<" ";
    //                    }
    //                    cout<<endl;
    //                    tmpout+=npwx;
    //                }
    //            }
    //            if(MY_RANK!=0) cout.setstate(ios::failbit);
    //        }
    //    }
    //}
    GlobalC::CHR.rho_mpi();
    for(int ir = 0; ir < nrxx ; ++ir)
    {
        tmprho = GlobalC::CHR.rho[0][ir] / GlobalC::ucell.omega;
        sto_rho[ir] = tmprho;
        sto_ne += tmprho;
    }
    sto_ne *= dr3;



#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE,&stodemet,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&sto_eband,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&sto_ne,1,MPI_DOUBLE,MPI_SUM,POOL_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&sto_ne,1,MPI_DOUBLE,MPI_SUM,PARAPW_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,sto_rho,nrxx,MPI_DOUBLE,MPI_SUM,PARAPW_WORLD);
#endif
    GlobalC::en.eband += sto_eband;
    GlobalC::en.demet += stodemet;
    GlobalC::en.demet *= Occupy::gaussian_parameter;

    cout.precision(12);
    GlobalV::ofs_running<<"Renormalize rho from ne = "<<sto_ne+KS_ne<<" to targetne = "<<targetne<<endl;

    double factor;
    if(abs(sto_ne) > 1e-20)
        factor = (targetne - KS_ne) / sto_ne;
    else
        factor = 1;

    if(GlobalV::MY_STOGROUP==0)
    {
        if(GlobalV::NBANDS > 0)
            ModuleBase::GlobalFunc::DCOPY(ksrho,GlobalC::CHR.rho[0],nrxx);
        else
            ModuleBase::GlobalFunc::ZEROS(GlobalC::CHR.rho[0],nrxx);
    }
    
    
    if(GlobalV::MY_STOGROUP == 0)
    for(int is = 0 ; is < 1; ++is)
    {
        for(int ir = 0; ir < nrxx ; ++ir)
        {
            GlobalC::CHR.rho[is][ir] += sto_rho[ir] * factor;
        }
    }

    
    delete [] sto_rho;
    ModuleBase::timer::tick("Stochastic_Iter","sum_stoband");
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

double Stochastic_Iter:: nxfd(double rawe)
{
    double Ebar = (Emin + Emax)/2;
	double DeltaE = (Emax - Emin)/2;
    double e = rawe * DeltaE + Ebar;
    double ne_mu = (e - mu) / Occupy::gaussian_parameter ;
    if(ne_mu > 40)
        return 0;
    else
        return e / (1 + exp(ne_mu));
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



//  void Stochastic_Iter:: test(const int & ik)
//  {
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

    int m = 5;
    complex<double> *chig1 = new complex<double> [wf.npwx*m];
    complex<double> *chig2 = new complex<double> [wf.npwx*m];
    complex<double> *hchig = new complex<double> [wf.npwx*m];
    //complex<double> *kswf = &wf.evc[ik](0,0);
    complex<double> *kswf = STO_WF.chi0[ik].c;
    
    stohchi.hchi_reciprocal(kswf,chig1,m);
    hm.hpw.h_psi( kswf , chig2, m);
    cout<<"==================ik="<<ik<<"====================================="<<endl;
    if(MY_RANK==0)
    for(int ib = 0 ; ib < m ; ++ib)
    {
        cout<<wf.npw<<" "<<wf.npwx<<endl;
        for(int i = 0; i<wf.npw;++i)
        {
            if(i % 100 == 0 )
                cout<<kswf[i+ib*wf.npwx]<<" "<<chig1[i+ib*wf.npwx]<<" "<<chig2[i+ib*wf.npwx]<<endl;
        }
        cout<<"-------------------------------------------------------------"<<endl;
    }
    cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
    MPI_Barrier(MPI_COMM_WORLD);

    if(MY_RANK==1)
    {
        cout.clear();
        
    for(int ib = 0 ; ib < m ; ++ib)
    {
        cout<<wf.npw<<" "<<wf.npwx<<endl;
        for(int i = 0; i<wf.npw;++i)
        {
            if(i % 100 == 0 )
                cout<<kswf[i+ib*wf.npwx]<<" "<<chig1[i+ib*wf.npwx]<<" "<<chig2[i+ib*wf.npwx]<<endl;
        }
        cout<<"-------------------------------------------------------------"<<endl;
    }
        cout.setstate(ios::failbit);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(MY_RANK==3)
    {
        cout.clear();
        
    for(int ib = 0 ; ib < m ; ++ib)
    {
        cout<<wf.npw<<" "<<wf.npwx<<endl;
        for(int i = 0; i<wf.npw;++i)
        {
            if(i % 100 == 0 )
                cout<<kswf[i+ib*wf.npwx]<<" "<<chig1[i+ib*wf.npwx]<<" "<<chig2[i+ib*wf.npwx]<<endl;
        }
        cout<<"-------------------------------------------------------------"<<endl;
    }
        cout.setstate(ios::failbit);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    cout<<"================================================================"<<endl;
    cout<<endl;
    delete[] chig1;
    delete[] chig2;
    if(ik == kv.nks-1) exit(0);*/
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

    //-------------------------------------------------------------
	//test orthogonal
    /*int npw=wf.npw;
    char transC='C';
    char transN='N';
    int nchip = STO_WF.nchip;
    complex<double> *wave = new complex<double> [nchip*npw];
    for(int i = 0; i < nchip*npw; ++i)
    {
        wave[i]=STO_WF.chi0[0].c[i];
    }
    complex<double> *sum = new complex<double> [nchip * nchip];
	zgemm_(&transC, &transN, &nchip, &nchip, &npw, &ONE, STO_WF.chi0[0].c, &npw, wave, &npw, &ZERO, sum, &nchip);
	Parallel_Reduce::reduce_complex_double_pool(sum, nchip * nchip);
	if(MY_RANK!=0) cout.clear();
    double abs2 = 0;
    for(int i=0;i<nchip * nchip;++i)
    {
        if(i%nchip==int(i%nchip)) continue;
        abs2+=norm(sum[i]);
    }
	cout<<abs2/nchip<<endl;
    delete [] sum;
    delete [] wave;
	if(MY_RANK!=0) cout.setstate(ios::failbit);*/
	//-------------------------------------------------------------

    
    //=====================================================
//  }
