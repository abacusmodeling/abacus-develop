#include "./sto_stress_pw.h"
#include "../module_base/timer.h"
#include "global.h"

void Sto_Stress_PW::cal_stress(ModuleBase::matrix& sigmatot, const psi::Psi<complex<double>>* psi_in, Stochastic_WF& stowf)
{
	ModuleBase::TITLE("Sto_Stress_PW","cal_stress");
	ModuleBase::timer::tick("Sto_Stress_PW","cal_stress");    

	sigmatot.create(3,3);
	ModuleBase::matrix sigmaxc(3,3);
	ModuleBase::matrix sigmahar(3,3);
	ModuleBase::matrix sigmakin(3,3);
	ModuleBase::matrix sigmaloc(3,3);
	ModuleBase::matrix sigmanl(3,3);
	ModuleBase::matrix sigmaewa(3,3);
	ModuleBase::matrix sigmaxcc(3,3);

	//kinetic contribution
	sto_stress_kin(sigmakin, psi_in, stowf);
	
	//hartree contribution
	stress_har(sigmahar, GlobalC::rhopw, 1);

    //ewald contribution
    stress_ewa(sigmaewa, GlobalC::rhopw, 1);

    //xc contribution: add gradient corrections(non diagonal)
    for(int i=0;i<3;++i)
	{
       sigmaxc(i,i) = - (GlobalC::en.etxc - GlobalC::en.vtxc) / GlobalC::ucell.omega;
    }
    stress_gga(sigmaxc);

    //local contribution
    stress_loc(sigmaloc, GlobalC::rhopw, 1);
    
    //nlcc
    stress_cc(sigmaxcc, GlobalC::rhopw, 1);
   
    //nonlocal
	sto_stress_nl(sigmanl,psi_in, stowf);


    for(int ipol=0;ipol<3;++ipol)
	{
        for(int jpol=0;jpol<3;++jpol)
		{
			sigmatot(ipol,jpol) = sigmakin(ipol,jpol) 
								+ sigmahar(ipol,jpol) 
								+ sigmanl(ipol,jpol) 
								+ sigmaxc(ipol,jpol) 
								+ sigmaxcc(ipol,jpol) 
								+ sigmaewa(ipol,jpol)
								+ sigmaloc(ipol,jpol);
        }
    }
    
	if(ModuleSymmetry::Symmetry::symm_flag)                          
	{
		GlobalC::symm.stress_symmetry(sigmatot,GlobalC::ucell);
	}

	bool ry = false;
	this->printstress_total(sigmatot, ry);

	if(GlobalV::TEST_STRESS) 
	{               
		GlobalV::ofs_running << "\n PARTS OF STRESS: " << endl;
		GlobalV::ofs_running << setiosflags(ios::showpos);
		GlobalV::ofs_running << setiosflags(ios::fixed) << setprecision(10) << endl;
		this->print_stress("KINETIC    STRESS",sigmakin,GlobalV::TEST_STRESS,ry);
		this->print_stress("LOCAL    STRESS",sigmaloc,GlobalV::TEST_STRESS,ry);
		this->print_stress("HARTREE    STRESS",sigmahar,GlobalV::TEST_STRESS,ry);
		this->print_stress("NON-LOCAL    STRESS",sigmanl,GlobalV::TEST_STRESS,ry);
		this->print_stress("XC    STRESS",sigmaxc,GlobalV::TEST_STRESS,ry);
		this->print_stress("EWALD    STRESS",sigmaewa,GlobalV::TEST_STRESS,ry);
		this->print_stress("NLCC    STRESS",sigmaxcc,GlobalV::TEST_STRESS,ry);
		this->print_stress("TOTAL    STRESS",sigmatot,GlobalV::TEST_STRESS,ry);
	}
	ModuleBase::timer::tick("Sto_Stress_PW","cal_stress");
	return;
    
}

void Sto_Stress_PW::sto_stress_kin(ModuleBase::matrix& sigma, const psi::Psi<complex<double>>* psi_in, Stochastic_WF& stowf)
{
	ModuleBase::TITLE("Sto_Stress_PW","cal_stress");
	ModuleBase::timer::tick("Sto_Stress_PW","cal_stress");
	double **gk;
	gk=new double* [3];
	ModuleBase::matrix s_kin(3,3,true);

		
	int npwx = GlobalC::wf.npwx;
	gk[0]= new double[npwx]; 
	gk[1]= new double[npwx];
	gk[2]= new double[npwx];
	double factor=ModuleBase::TWO_PI/GlobalC::ucell.lat0;
	int nksbands = psi_in->get_nbands();
	if(GlobalV::MY_STOGROUP != 0) nksbands = 0;

	for(int ik=0;ik<GlobalC::kv.nks;++ik)
	{
		const int nstobands = stowf.nchip[ik];
		const int nbandstot = nstobands + nksbands;
		const int npw = GlobalC::kv.ngk[ik];
		for(int i=0;i<npw;++i)
		{
			gk[0][i]=GlobalC::wfcpw->getgpluskcar(ik,i)[0]*factor;
			gk[1][i]=GlobalC::wfcpw->getgpluskcar(ik,i)[1]*factor;
			gk[2][i]=GlobalC::wfcpw->getgpluskcar(ik,i)[2]*factor;
		}

		//kinetic contribution

		for(int l=0;l<3;++l)
		{
			for(int m=0;m<l+1;++m)
			{
				for(int ibnd=0;ibnd<nbandstot;++ibnd)
				{
					if(ibnd < nksbands)
					{
						for(int i=0;i<npw;++i)
						{
							std::complex<double> p = psi_in->operator()(ik, ibnd, i);
							double np = p.real() * p.real() + p.imag() * p.imag();
							s_kin(l,m) +=
								GlobalC::wf.wg(ik, ibnd)*gk[l][i]*gk[m][i] * np;
						}
					}
					else
					{
						for(int i=0;i<npw;++i)
						{
							std::complex<double> p = stowf.shchi[ik](ibnd-nksbands, i);
							double np = p.real() * p.real() + p.imag() * p.imag();
							s_kin(l,m) += GlobalC::kv.wk[ik]*gk[l][i]*gk[m][i] * np;
						}
					}
				}
			}
		}
		   
	}
		

	for(int l=0;l<3;++l)
	{
		for(int m=0;m<l;++m)
		{
			s_kin(m,l)=s_kin(l,m);
		}
	}


	for(int l=0;l<3;++l)
	{
		for(int m=0;m<3;++m)
		{
			s_kin(l,m) *= ModuleBase::e2/GlobalC::ucell.omega;
		}
	}
	

	Parallel_Reduce::reduce_double_all( s_kin.c, 9 );

	for(int l=0;l<3;++l)
	{
		for(int m=0;m<3;++m)
		{
			sigma(l,m) = s_kin(l,m);
		}
	}
	//do symmetry
	if(ModuleSymmetry::Symmetry::symm_flag)                          
	{
		GlobalC::symm.stress_symmetry(sigma,GlobalC::ucell);
	}
	delete[] gk[0];
	delete[] gk[1];
	delete[] gk[2];
	delete[] gk;
	ModuleBase::timer::tick("Sto_Stress_PW","cal_stress");
		
	return;
}

void Sto_Stress_PW::sto_stress_nl(ModuleBase::matrix& sigma, const psi::Psi<complex<double>>* psi_in, Stochastic_WF& stowf){
	ModuleBase::TITLE("Sto_Stress_Func","stres_nl");
	ModuleBase::timer::tick("Sto_Stress_Func","stres_nl");
	
	const int nkb = GlobalC::ppcell.nkb;
	if(nkb == 0) 
	{
		ModuleBase::timer::tick("Stress_Func","stres_nl");
		return;
	}

	ModuleBase::matrix sigmanlc(3,3,true);
	int* nchip = stowf.nchip;
	const int npwx = GlobalC::wf.npwx;
	int nksbands = psi_in->get_nbands();
	if(GlobalV::MY_STOGROUP != 0) nksbands = 0;
	

	// vkb1: |Beta(nkb,npw)><Beta(nkb,npw)|psi(nbnd,npw)>
	ModuleBase::ComplexMatrix vkb1( nkb, npwx );
	ModuleBase::ComplexMatrix vkb0[3];
	for(int i=0;i<3;++i){
		vkb0[i].create(nkb, npwx);
	}
	ModuleBase::ComplexMatrix vkb2( nkb, npwx );
	
    for (int ik = 0;ik < GlobalC::kv.nks;ik++)
    {
		const int nstobands = stowf.nchip[ik];
		const int nbandstot = nstobands + nksbands;
		const int npw = GlobalC::kv.ngk[ik];
		// dbecp: conj( -iG * <Beta(nkb,npw)|psi(nbnd,npw)> )
		ModuleBase::ComplexMatrix dbecp( nbandstot, nkb);
		ModuleBase::ComplexMatrix becp( nbandstot, nkb);
		  
		if (GlobalV::NSPIN==2) GlobalV::CURRENT_SPIN = GlobalC::kv.isk[ik];
		// generate vkb
		if (GlobalC::ppcell.nkb > 0)
		{
			GlobalC::ppcell.getvnl(ik, GlobalC::ppcell.vkb);
		}

		// get becp according to wave functions and vkb
		// important here ! becp must set zero!!
		// vkb: Beta(nkb,npw)
		// becp(nkb,nbnd): <Beta(nkb,npw)|psi(nbnd,npw)>
		becp.zero_out();
		char transa = 'C';
        char transb = 'N';
		psi_in[0].fix_k(ik);
		//KS orbitals
		int npmks = GlobalV::NPOL * nksbands;
		zgemm_(&transa,&transb,&nkb,&npmks,&npw,&ModuleBase::ONE,
				GlobalC::ppcell.vkb.c,&npwx,
            	psi_in->get_pointer(),&npwx,
            	&ModuleBase::ZERO,becp.c,&nkb);
		//stochastic orbitals
		int npmsto = GlobalV::NPOL * nstobands;
		zgemm_(&transa,&transb,&nkb,&npmsto,&npw,&ModuleBase::ONE,
				GlobalC::ppcell.vkb.c,&npwx,
            	stowf.shchi[ik].c,&npwx,
            	&ModuleBase::ZERO,&becp(nksbands,0),&nkb);
		
		Parallel_Reduce::reduce_complex_double_pool( becp.c, becp.size);
			
		for(int i=0;i<3;++i) {
			get_dvnl1(vkb0[i],ik,i);
		}
				
		get_dvnl2(vkb2,ik);

		ModuleBase::Vector3<double> qvec;
		double* qvec0[3];
		qvec0[0] = &(qvec.x);
		qvec0[1] = &(qvec.y);
		qvec0[2] = &(qvec.z);
				
		for (int ipol = 0; ipol<3; ++ipol)
		{
			for(int jpol = 0; jpol < ipol+1; ++jpol)
			{
				dbecp.zero_out();    
				vkb1.zero_out();
				for (int i = 0;i < nkb;++i)
				{
					std::complex<double>* pvkb0i = &vkb0[ipol](i, 0);
					std::complex<double>* pvkb0j = &vkb0[jpol](i, 0);
					std::complex<double>* pvkb1 = &vkb1(i, 0);
					// third term of dbecp_noevc
					for (int ig = 0; ig < npw; ++ig) 
					{
						qvec = GlobalC::wfcpw->getgpluskcar(ik, ig);
						pvkb1[ig] += 0.5 * qvec0[ipol][0] * pvkb0j[ig] +
									0.5 * qvec0[jpol][0] * pvkb0i[ig];
					} // end ig	  
				}//end nkb
 

				ModuleBase::ComplexMatrix dbecp_noevc(nkb, GlobalC::wf.npwx, true);
				for (int i = 0; i < nkb; ++i) 
				{
					std::complex<double>* pdbecp_noevc = &dbecp_noevc(i, 0);
					std::complex<double>* pvkb = &vkb1(i, 0);
					// first term
					for (int ig = 0; ig < npw;++ig) 
					{
						pdbecp_noevc[ig] -= 2.0 * pvkb[ig];
					}
					// second termi
					if (ipol == jpol)
					{
						pvkb = &GlobalC::ppcell.vkb(i, 0);
						for (int ig = 0; ig < npw;++ig) 
						{
							pdbecp_noevc[ig] -= pvkb[ig];
						}
					}
					// third term
					pvkb = &vkb2(i,0);
					for (int ig = 0; ig < npw;++ig) 
					{
						qvec =	GlobalC::wfcpw->getgpluskcar(ik, ig);
						double qm1;
						if(qvec.norm2() > 1e-16) qm1 = 1.0 / qvec.norm(); 
						else qm1 = 0; 
						pdbecp_noevc[ig] -= 2.0 * pvkb[ig] * qvec0[ipol][0] * 
							qvec0[jpol][0] * qm1 *	GlobalC::ucell.tpiba;
					} // end ig
				}     // end i

//              //KS orbitals
				zgemm_(&transa,&transb,&nkb,&npmks,&npw,&ModuleBase::ONE,
						dbecp_noevc.c,&npwx,
        		    	psi_in->get_pointer(),&npwx,
        		    	&ModuleBase::ZERO,dbecp.c,&nkb);
				//stochastic orbitals
				zgemm_(&transa,&transb,&nkb,&npmsto,&npw,&ModuleBase::ONE,
						dbecp_noevc.c,&npwx,
        		    	stowf.shchi[ik].c,&npwx,
        		    	&ModuleBase::ZERO,&dbecp(nksbands, 0),&nkb);
				
				for (int ib=0; ib<nbandstot; ++ib)
				{
					double fac;
					if(ib < nksbands)
						fac = GlobalC::wf.wg(ik, ib);
					else
						fac = GlobalC::kv.wk[ik];
					int iat = 0;
					int sum = 0;
					for (int it=0; it<GlobalC::ucell.ntype; ++it)
					{
						const int Nprojs = GlobalC::ucell.atoms[it].nh;
						for (int ia=0; ia<GlobalC::ucell.atoms[it].na; ++ia)
						{
							for (int ip=0; ip<Nprojs; ++ip)
							{
								double ps = GlobalC::ppcell.deeq(GlobalV::CURRENT_SPIN, iat, ip, ip) ;
								const int inkb = sum + ip;
								//out<<"\n ps = "<<ps;

							 
								const double dbb = ( conj( dbecp( ib, inkb) ) * becp( ib, inkb) ).real();
								sigmanlc(ipol, jpol) -= ps * fac * dbb;
							 
							}//end ip
							++iat;        
							sum+=Nprojs;
						}//ia
					} //end it
				} //end band
			}//end jpol
		}//end ipol
	}// end ik

	for(int l=0;l<3;++l)
	{
		for(int m=0;m<3;++m)
		{
			if(m>l) 
			{
				sigmanlc(l,m) = sigmanlc(m,l);
			}
		}
	}
	// sum up forcenl from all processors
	Parallel_Reduce::reduce_double_all( sigmanlc.c, 9);

        
	for (int ipol = 0; ipol<3; ++ipol)
	{
		for(int jpol = 0; jpol < 3; ++jpol)
		{
			sigmanlc(ipol, jpol) *= 1.0 / GlobalC::ucell.omega;
		}
	}
	
	for (int ipol = 0; ipol<3; ++ipol)
	{
		for(int jpol = 0; jpol < 3; ++jpol)
		{
			sigma(ipol,jpol) = sigmanlc(ipol,jpol);
		}
	}
	//do symmetry
	if(ModuleSymmetry::Symmetry::symm_flag)                          
	{
		GlobalC::symm.stress_symmetry(sigma,GlobalC::ucell);
	}
	
	//  this->print(ofs_running, "nonlocal stress", stresnl);
	ModuleBase::timer::tick("Sto_Stress_Func","stres_nl");
	return;
}
