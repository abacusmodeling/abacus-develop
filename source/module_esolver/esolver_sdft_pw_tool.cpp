#include "./esolver_sdft_pw.h"
#include "module_hsolver/diago_iter_assist.h"
#include "module_hsolver/hsolver_pw_sdft.h"
#include "module_base/timer.h"
#include "module_base/constants.h"
#include "module_base/vector3.h"
#include "module_base/complexmatrix.h"
#include "module_base/global_variable.h"
#include "module_base/global_function.h"
#include "src_pw/global.h"

#define TWOSQRT2LN2 2.354820045030949 //FWHM = 2sqrt(2ln2) * \sigma
#define FACTOR 1.839939223835727e7
namespace ModuleESolver
{

void ESolver_SDFT_PW::sKG(const int nche_KG, const double fwhmin, const double wcut, 
                          const double dw_in, const int times)
{
    ModuleBase::timer::tick(this->classname,"sKG");
    cout<<"Calculating conductivity...."<<endl;
    //------------------------------------------------------------------
    // hbar = 6.62607015e-34/2pi
    // e    = 1.6021766208e-19
    // a    = 5.2917721067e-11
    // m    = 9.1093897e-31
    // 1 ha   = hbar^2/m/a^2/e  = 27.21136857564 eV
    // 1 ry   = hbar^2/2m/a^2/e = 13.60568428782 eV = 2.17987092759e-18 J
    // 1 t(ry^-1) = hbar/ry/e    = 4.837771834548454e-17 s
    // factor = hbar*e^2/a^5/m^2*t^2  = 1.839939223835727e+07  (1 a.u. = 1.84e7 Sm^-1)
    // 1 a.u. = factor*(2.17987092759e-18)^2/e^2 = 3.40599696130e+09 Wm^-1
    // k = 1.380649e-23
    // e/k = 11604.518026 , 1 eV = 11604.5 K
    //------------------------------------------------------------------
    double factor = FACTOR;
    int nw = ceil(wcut/dw_in);
    double dw =  dw_in / ModuleBase::Ry_to_eV; //converge unit in eV to Ry 
    double sigma = fwhmin / TWOSQRT2LN2 / ModuleBase::Ry_to_eV;
    double dt = ModuleBase::PI/(dw*nw)/times ; //unit in a.u., 1 a.u. = 4.837771834548454e-17 s
    int nt = ceil(sqrt(20)/sigma/dt);
    cout<<"nw: "<<nw<<" ; dw: "<<dw*ModuleBase::Ry_to_eV<<" eV"<<endl;
    cout<<"nt: "<<nt<<" ; dt: "<<dt<<" a.u.(ry^-1)"<<endl;
    assert(nw >= 1);
    assert(nt >= 1);

    ModuleBase::Chebyshev<double> che(this->nche_sto);
    ModuleBase::Chebyshev<double> chet(nche_KG);
    const int npwx = GlobalC::wf.npwx;
    const double tpiba = GlobalC::ucell.tpiba;
    Stochastic_Iter& stoiter = ((hsolver::HSolverPW_SDFT*)phsol)->stoiter;
    Stochastic_hchi& stohchi = stoiter.stohchi;
    const int nk = GlobalC::kv.nks;

    //------------------------------
    //      Convergence test
    //------------------------------
    bool change = false;
    for (int ik = 0;ik < nk; ++ik)
	{
        int ntest = 2;
        if (this->stowf.nchip[ik] < ntest) 
	    {
	    	ntest = this->stowf.nchip[ik];
	    }

        this->phami->updateHk(ik);
        const int npw = GlobalC::kv.ngk[ik];

        complex<double> *pchi = new complex<double> [npw*4*ntest];
        //\sqrt{f}|\chi>
        for(int i = 0 ; i < ntest ; ++i)
            for(int ig = 0 ; ig < npw ; ++ig)
                pchi[i*npw+ig] = this->stowf.shchi[ik](i,ig);
        //(G+k)\sqrt{f}|\chi>
        for(int a = 0 ; a < 3 ; ++a)
        {
            for(int i = 0 ; i < ntest ; ++i)
            {
                for(int ig = 0 ; ig < npw ; ++ig)
                {
                    ModuleBase::Vector3<double> v3 = GlobalC::wfcpw->getgpluskcar(ik,i);
                    pchi[((a+1)*ntest + i)*npw + ig] = this->stowf.shchi[ik](i,ig)*v3[a] * tpiba;
                }
            }
        }
        for(int i = 0; i < ntest * 4 ; ++i)
        {
            while(1)
            {
                bool converge;
                if(this->nche_sto > nche_KG)
                {
                    converge= che.checkconverge(&stohchi, &Stochastic_hchi::hchi_reciprocal, 
	        	    	&pchi[i*npw], npw, stohchi.Emax, stohchi.Emin, 5.0);
                }
                else
                {
                    converge= chet.checkconverge(&stohchi, &Stochastic_hchi::hchi_reciprocal, 
	        	    	&pchi[i*npw], npw, stohchi.Emax, stohchi.Emin, 5.0);
                }
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

        if(ik == nk-1)
        {
            stoiter.stofunc.Emax = stohchi.Emax;
            stoiter.stofunc.Emin = stohchi.Emin;
#ifdef __MPI
            MPI_Allreduce(MPI_IN_PLACE, &stoiter.stofunc.Emax, 1, MPI_DOUBLE, MPI_MAX , MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &stoiter.stofunc.Emin, 1, MPI_DOUBLE, MPI_MIN , MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &change, 1, MPI_CHAR, MPI_LOR , MPI_COMM_WORLD);
#endif
            stohchi.Emin = stoiter.stofunc.Emin;
            stohchi.Emax = stoiter.stofunc.Emax;
            if(change)
	        {
	        	GlobalV::ofs_running<<"New Emax "<<stohchi.Emax<<" ; new Emin "<<stohchi.Emin<<std::endl;
	        }
            change = false;
        }
        delete[] pchi;
    }

    //------------------------------------------------------------------
    //                    Calculate
    //------------------------------------------------------------------
    const double mu = GlobalC::en.ef;
    stoiter.stofunc.mu = mu;
    double * ct11 = new double[nt];
    double * ct12 = new double[nt];
    double * ct22 = new double[nt];
    ModuleBase::GlobalFunc::ZEROS(ct11,nt);
    ModuleBase::GlobalFunc::ZEROS(ct12,nt);
    ModuleBase::GlobalFunc::ZEROS(ct22,nt);
    double *gxyz = new double[npwx*3];

    stoiter.stofunc.t = dt;
    chet.calcoef_pair(&stoiter.stofunc, &Sto_Func<double>::ncos, &Sto_Func<double>::n_sin);
    cout<<"Chebyshev precision: "<<abs(chet.coef_complex[nche_KG-1]/chet.coef_complex[0])<<endl;
    ofstream cheofs("Chebycoef");
    for(int i  = 0 ; i < nche_KG ; ++i)
    {
        cheofs<<setw(5)<<i<<setw(20)<<abs(chet.coef_complex[i]/chet.coef_complex[0])<<endl;
    }
    cheofs.close();
    
    //KS
    int ksbandper = GlobalV::NBANDS / GlobalV::NSTOGROUP;
    int startband = ksbandper * GlobalV::MY_STOGROUP;
    if(GlobalV::MY_STOGROUP < GlobalV::NBANDS % GlobalV::NSTOGROUP)
    {  
        ++ksbandper;
        startband += GlobalV::MY_STOGROUP;
    }
    else
    {
        startband += GlobalV::NBANDS % GlobalV::NSTOGROUP;
    }
    ModuleBase::ComplexMatrix kswf(ksbandper,npwx,false);
    double *ksen = new double[ksbandper];   

    ModuleBase::timer::tick(this->classname,"kloop");
    for (int ik = 0;ik < nk;++ik)
	{
        if(nk > 1) 
        {
            this->phami->updateHk(ik);
        }
        const int npw = GlobalC::kv.ngk[ik];

        int nchip = this->stowf.nchip[ik];
        int totbands = nchip + ksbandper;
        complex<double> *leftv3 = new complex<double>[npwx*totbands*3];
        complex<double> *jleftv3 = new complex<double>[npwx*totbands*3];
        complex<double> *hleftv = new complex<double>[npwx*nchip];
        complex<double> *prightv3 = new complex<double>[npwx*totbands*3];
        complex<double> *jrightv3 = new complex<double>[npwx*totbands*3];
        complex<double> *rightv = new complex<double>[npwx*totbands];
        

        for(int ib = 0; ib < ksbandper ; ++ib)
        {
            for(int ig = 0 ; ig < npw ; ++ig)
            {
                kswf(ib,ig) = this->psi[0](ik,ib+startband,ig);
            }
            ksen[ib] = this->pelec->ekb(ik, ib+startband);
        }
        for (int ig = 0;ig < npw;ig++)
        {
            ModuleBase::Vector3<double> v3 = GlobalC::wfcpw->getgpluskcar(ik,ig);
            gxyz[ig] = v3.x * tpiba;
            gxyz[ig+npwx] = v3.y * tpiba;
            gxyz[ig+2*npwx] = v3.z * tpiba;
        }
        //-----------------------------------------------------------
        //               sto conductivity
        //-----------------------------------------------------------
        //j|left> = [(Hp+pH)/2 - mu*p]|left>
        this->phami->hPsi(this->stowf.shchi[ik].c, hleftv, nchip*npwx); //hleftv is a transitional array to get jleftv3
        for(int a = 0 ; a < 3 ; ++a)
        {
            for(int ib = 0; ib < nchip ; ++ib) // sto orbitals
		    {
		    	for (int ig = 0; ig < npw; ++ig)
		    	{
		    		leftv3[a*totbands*npwx + ib*npwx + ig] = gxyz[a*npwx + ig] * this->stowf.shchi[ik](ib,ig);
		    	}
		    }
            for(int ib = nchip ; ib < totbands ; ++ib) // KS orbital
            {
                for (int ig = 0; ig < npw; ++ig)
		    	{
		    		leftv3[a*totbands*npwx + ib*npwx + ig] = gxyz[a*npwx + ig] * kswf(ib-nchip,ig);
		    	}
            }
        }
        this->phami->hPsi(leftv3, jleftv3, totbands*3*npwx); //leftv3 is a transitional array to get jleftv3
        for(int a = 0 ; a < 3 ; ++a)
        {
            for (int ib = 0; ib < nchip ; ++ib) // sto orbitals
		    {
		    	for (int ig = 0; ig < npw; ++ig)
		    	{
                    int iaibig = a*totbands*npwx + ib*npwx + ig;
                    int ibig = ib *npwx + ig;
                    double ga = gxyz[a*npwx + ig];
                    jleftv3[iaibig] = (jleftv3[iaibig] + hleftv[ibig] * ga ) / 2.0 - mu * this->stowf.shchi[ik](ib,ig) * ga;
                }
            }
            for (int ib = nchip; ib < totbands ; ++ib) //KS orbitals
		    {
		    	for (int ig = 0; ig < npw; ++ig)
		    	{
                    int iaibig = a*totbands*npwx + ib*npwx + ig; 
                    double ga = gxyz[a*npwx + ig];
                    complex<double> wfikib = kswf(ib-nchip,ig);
                    jleftv3[iaibig] = (jleftv3[iaibig] + ga * ksen[ib-nchip] * wfikib ) / 2.0 - mu * wfikib * ga;
                }
            }
        }

        //p|left> Note: p^+=p, px=-id/dx  p|k>=k|k>
        for(int a = 0 ; a < 3 ; ++a)
        {
            for (int ib = 0; ib < nchip ; ++ib) // sto orbitals
		    {
		    	for (int ig = 0; ig < npw; ++ig)
		    	{
		    		leftv3[a*totbands*npwx + ib*npwx + ig] = gxyz[a*npwx + ig] * this->stowf.shchi[ik](ib,ig);
		    	}
		    }
            for(int ib = nchip; ib < totbands; ++ib) // KS orbitals
            {
                for (int ig = 0; ig < npw; ++ig)
		    	{
		    		leftv3[a*totbands*npwx + ib*npwx + ig] = gxyz[a*npwx + ig] * kswf(ib-nchip,ig);
		    	}
            }
        }

        ModuleBase::timer::tick(this->classname,"chebyshev");
        //(1-f)|left>
        che.calcoef_real(&stoiter.stofunc,&Sto_Func<double>::n_fd);
        che.calfinalvec_real(&stohchi, &Stochastic_hchi::hchi_reciprocal, leftv3, leftv3, npw, npwx, totbands*3);
        che.calfinalvec_real(&stohchi, &Stochastic_hchi::hchi_reciprocal, jleftv3, jleftv3, npw, npwx, totbands*3);
        ModuleBase::timer::tick(this->classname,"chebyshev");

                
        cout<<"ik="<<ik<<": ";
        ct11[0] = 0;
        for (int it = 1 ;it < nt ; ++it)
        {
            if(it%20==0) cout<<it<<" ";
            ModuleBase::timer::tick(this->classname,"chebyshev");
            //exp(-iH*dt/h)|left> Note: exp(iH*dt/h)^+=exp(-iHt/h)   
            chet.calfinalvec_complex(&stohchi, &Stochastic_hchi::hchi_reciprocal, leftv3, leftv3, npw, npwx, totbands*3);
            chet.calfinalvec_complex(&stohchi, &Stochastic_hchi::hchi_reciprocal, jleftv3, jleftv3, npw, npwx, totbands*3);


            //exp(-iH*dt/h)|right>
            if(it == 1) chet.calfinalvec_complex(&stohchi, &Stochastic_hchi::hchi_reciprocal, this->stowf.shchi[ik].c, rightv, npw, npwx, nchip);
            else        chet.calfinalvec_complex(&stohchi, &Stochastic_hchi::hchi_reciprocal, rightv, rightv, npw, npwx, nchip);
            ModuleBase::timer::tick(this->classname,"chebyshev");
            for(int ib = nchip ; ib < totbands ; ++ib)
            {
                for (int ig = 0; ig < npw; ++ig)
		    	{
		    		rightv[ib*npwx + ig] = exp(-ModuleBase::IMAG_UNIT * ksen[ib-nchip] * double(it) * dt) * kswf(ib-nchip,ig);
		    	}
            }

            //-i[(Hp+pH)/2 - mu*p]|right>
            this->phami->hPsi(rightv, hleftv, nchip*npwx); //hleftv is a transitional array to get jleftv3
            for(int a = 0 ; a < 3 ; ++a)
            {
                for (int ib = 0; ib < totbands ; ++ib) //KS + sto
		        {
		        	for (int ig = 0; ig < npw; ++ig)
		        	{
		        		prightv3[a*totbands*npwx + ib*npwx + ig] = gxyz[a*npwx + ig] * rightv[ib*npwx + ig];
		        	}
		        }
            }
            this->phami->hPsi(prightv3, jrightv3, totbands*3*npwx); //prightv3 is a transitional array to get jleftv3
            for(int a = 0 ; a < 3 ; ++a)
            {
                for (int ib = 0; ib < nchip ; ++ib) //sto orbitals
		        {
		        	for (int ig = 0; ig < npw; ++ig)
		        	{
                        int iaibig = a*totbands*npwx + ib*npwx + ig;
                        int ibig = ib *npwx + ig;
                        double ga = gxyz[a*npwx + ig];
                        jrightv3[iaibig] = -ModuleBase::IMAG_UNIT * ((jrightv3[iaibig] + hleftv[ibig] * ga ) / double(2) - mu * rightv[ib*npwx + ig] * ga);
                    }
                }
                for (int ib = nchip; ib < totbands ; ++ib) //KS orbitals
		        {
		        	for (int ig = 0; ig < npw; ++ig)
		        	{
                        int iaibig = a*totbands*npwx + ib*npwx + ig; 
                        double ga = gxyz[a*npwx + ig];
                        complex<double> wfikib = rightv[ib*npwx + ig];
                        jrightv3[iaibig] = -ModuleBase::IMAG_UNIT * ((jrightv3[iaibig] + ga * ksen[ib-nchip] * wfikib ) / double(2) - mu * wfikib * ga);
                    }
                }
            }

            //-ip|right>
            for(int a = 0 ; a < 3 ; ++a)
            {      
                for (int ib = 0; ib < totbands ; ++ib) //sto + KS
		        {
		    	    for (int ig = 0; ig < npw; ++ig)
		    	    {
		    		    prightv3[a*totbands*npwx + ib*npwx + ig] = -ModuleBase::IMAG_UNIT * gxyz[a*npwx + ig] * rightv[ib*npwx + ig];
		    	    }
		        }
            }

            //Im<left|p|right>=Re<left|-ip|right>
            for(int a = 0 ; a < 3 ; ++a)
            {
                for(int ib = 0 ; ib < nchip ; ++ib)
                {
                    ct11[it] += ModuleBase::GlobalFunc::ddot_real(npw,&leftv3[a*totbands*npwx + ib*npwx],&prightv3[a*totbands*npwx + ib*npwx],false) * GlobalC::kv.wk[ik] / 2,0;
                    ct12[it] -= ModuleBase::GlobalFunc::ddot_real(npw,&jleftv3[a*totbands*npwx + ib*npwx],&prightv3[a*totbands*npwx + ib*npwx],false) * GlobalC::kv.wk[ik] / 2.0;
                    ct22[it] += ModuleBase::GlobalFunc::ddot_real(npw,&jleftv3[a*totbands*npwx + ib*npwx],&jrightv3[a*totbands*npwx + ib*npwx],false) * GlobalC::kv.wk[ik] / 2.0;
                }
                for(int ib = nchip ; ib < totbands ; ++ib)
                {
                    double ei = ksen[ib-nchip];
                    double fi = stoiter.stofunc.fd(ei);
                    ct11[it] += ModuleBase::GlobalFunc::ddot_real(npw,&leftv3[a*totbands*npwx + ib*npwx],&prightv3[a*totbands*npwx + ib*npwx],false) * GlobalC::kv.wk[ik] / 2.0 * fi;
                    ct12[it] -= ModuleBase::GlobalFunc::ddot_real(npw,&jleftv3[a*totbands*npwx + ib*npwx],&prightv3[a*totbands*npwx + ib*npwx],false) * GlobalC::kv.wk[ik] / 2.0 * fi;
                    ct22[it] += ModuleBase::GlobalFunc::ddot_real(npw,&jleftv3[a*totbands*npwx + ib*npwx],&jrightv3[a*totbands*npwx + ib*npwx],false) * GlobalC::kv.wk[ik] / 2.0 * fi;
                }
            }
        }
        cout<<endl;

        delete [] leftv3;
        delete [] jleftv3;
        delete [] hleftv;
        delete [] prightv3;
        delete [] jrightv3;
        delete [] rightv;
    }
    ModuleBase::timer::tick(this->classname,"kloop");
    delete[] gxyz;
    delete[] ksen;
#ifdef __MPI
    MPI_Allreduce(MPI_IN_PLACE,ct11,nt,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,ct12,nt,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,ct22,nt,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
   
    //------------------------------------------------------------------
    //                    Output
    //------------------------------------------------------------------
    if(GlobalV::MY_RANK==0)
    {
        calcondw(nt,dt,fwhmin,wcut,dw_in,ct11,ct12,ct22);
    }
    delete[] ct11;
    delete[] ct12;
    delete[] ct22;
    // delete[] cterror;
    ModuleBase::timer::tick(this->classname,"sKG");
}

void ESolver_SDFT_PW::calcondw(const int nt,const double dt,const double fwhmin,const double wcut,const double dw_in,double*ct11,double*ct12,double *ct22)
{
    double factor = FACTOR;
    int nw = ceil(wcut/dw_in);
    double dw =  dw_in / ModuleBase::Ry_to_eV; //converge unit in eV to Ry 
    double sigma = fwhmin / TWOSQRT2LN2 / ModuleBase::Ry_to_eV;
    ofstream ofscond("je-je.txt");
    ofscond<<setw(8)<<"#t(a.u.)"<<setw(15)<<"c11(t)"<<setw(15)<<"c12(t)"<<setw(15)<<"c22(t)"<<setw(15)<<"decay"<<endl;
	for(int it = 0; it < nt; ++it)
	{
		ofscond <<setw(8)<<(it)*dt<<setw(15)<<-2*ct11[it]<<setw(15)<<-2*ct12[it]<<setw(15)<<-2*ct22[it]<<setw(15)<<exp(-double(1)/2*sigma*sigma*pow((it)*dt,2))<<endl;
	}
    ofscond.close();
    double * cw11 = new double [nw];
    double * cw12 = new double [nw];
    double * cw22 = new double [nw];
    double * kappa = new double [(int)ceil(wcut/dw_in)];
    ModuleBase::GlobalFunc::ZEROS(cw11,nw);
    ModuleBase::GlobalFunc::ZEROS(cw12,nw);
    ModuleBase::GlobalFunc::ZEROS(cw22,nw);
    for(int iw = 0 ; iw < nw ; ++iw )
    {
        for(int it = 0 ; it < nt ; ++it)
        {
            cw11[iw] += -2 * ct11[it] * sin( -(iw+0.5) * dw * it *dt) * exp(-double(1)/2*sigma*sigma*pow((it)*dt,2) ) / (iw+0.5) /dw * dt ;
            cw12[iw] += -2 * ct12[it] * sin( -(iw+0.5) * dw * it *dt) * exp(-double(1)/2*sigma*sigma*pow((it)*dt,2) ) / (iw+0.5) /dw * dt ;
            cw22[iw] += -2 * ct22[it] * sin( -(iw+0.5) * dw * it *dt) * exp(-double(1)/2*sigma*sigma*pow((it)*dt,2) ) / (iw+0.5) /dw * dt ;
        }
    }
    ofscond.open("Onsager.txt");
    ofscond<<setw(8)<<"## w(eV) "<<setw(20)<<"sigma (Sm^-1)"<<setw(20)<<"kappa (W(mK)^-1)"<<setw(20)<<"L12/e (Am^-1)"<<setw(20)<<"L22/e^2 (Wm^-1)"<<endl;
    for(int iw = 0; iw < nw; ++iw)
	{
        cw11[iw] *= double(2)/3/GlobalC::ucell.omega* factor; //unit in Sm^-1
        cw12[iw] *= double(2)/3/GlobalC::ucell.omega* factor * 2.17987092759e-18/1.6021766208e-19; //unit in Am^-1
        cw22[iw] *= double(2)/3/GlobalC::ucell.omega* factor * pow(2.17987092759e-18/1.6021766208e-19,2); //unit in Wm^-1
        kappa[iw] = (cw22[iw]-pow(cw12[iw],2)/cw11[iw])/stoiter.stofunc.tem/ModuleBase::Ry_to_eV/11604.518026;
	    ofscond <<setw(8)<<(iw+0.5)*dw * ModuleBase::Ry_to_eV <<setw(20)<<cw11[iw] <<setw(20)<<kappa[iw]<<setw(20)<<cw12[iw] <<setw(20)<<cw22[iw]<<endl;
	}
    cout<<setprecision(6)<<"DC electrical conductivity: "<<cw11[0] - (cw11[1] - cw11[0]) * 0.5<<" Sm^-1"<<endl;
    cout<<setprecision(6)<<"Thermal conductivity: "<<kappa[0] - (kappa[1] - kappa[0]) * 0.5<<" Wm^-1"<<endl;;
    ofscond.close();
    
    
    delete[] cw11;
    delete[] cw12;
    delete[] cw22;
    delete[] kappa;

}

}