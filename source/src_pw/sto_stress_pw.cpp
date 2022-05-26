#include "./sto_stress_pw.h"
#include "vdwd2.h"
#include "vdwd3.h"
#include "../module_base/timer.h"
#include "global.h"

void Sto_Stress_PW::cal_stress(ModuleBase::matrix& sigmatot, Stochastic_WF& stowf)
{
	ModuleBase::TITLE("Sto_Stress_PW","cal_stress");
	ModuleBase::timer::tick("Sto_Stress_PW","cal_stress");    

	sigmatot.create(3,3);
	ModuleBase::matrix sigmaxc(3,3);
	ModuleBase::matrix sigmahar(3,3);
	ModuleBase::matrix sigmakin(3,3);
	ModuleBase::matrix sto_sigmakin(3,3);
	ModuleBase::matrix sigmaloc(3,3);
	ModuleBase::matrix sigmanl(3,3);
	ModuleBase::matrix sto_sigmanl(3,3);
	ModuleBase::matrix sigmaewa(3,3);
	ModuleBase::matrix sigmaxcc(3,3);

	//kinetic contribution
	if(GlobalV::NBANDS > 0 && GlobalV::MY_STOGROUP == 0) stress_kin(sigmakin);
	sto_stress_kin(sto_sigmakin, stowf);
	sigmakin = sigmakin + sto_sigmakin;
	
	//hartree contribution
	stress_har(sigmahar, 1);

    //ewald contribution
    stress_ewa(sigmaewa, 1);

    //xc contribution: add gradient corrections(non diagonal)
    for(int i=0;i<3;i++)
	{
       sigmaxc(i,i) = - (GlobalC::en.etxc - GlobalC::en.vtxc) / GlobalC::ucell.omega;
    }
    stress_gga(sigmaxc);

    //local contribution
    stress_loc(sigmaloc, 1);
    
    //nlcc
    stress_cc(sigmaxcc, 1);
   
    //nonlocal
	if(GlobalV::NBANDS > 0 && GlobalV::MY_STOGROUP == 0) stress_nl(sigmanl);
	sto_stress_nl(sto_sigmanl, stowf);
	sigmanl = sigmanl + sto_sigmanl;


    for(int ipol=0;ipol<3;ipol++)
	{
        for(int jpol=0;jpol<3;jpol++)
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

void Sto_Stress_PW::sto_stress_kin(ModuleBase::matrix& sigma, Stochastic_WF& stowf)
{
	ModuleBase::TITLE("Sto_Stress_PW","cal_stress");
	ModuleBase::timer::tick("Sto_Stress_PW","cal_stress");
	double **gk;
	gk=new double* [3];
	double tbsp,gk2,arg;
	int ik,l,m,i,j,ibnd,is;
	int npw;
	double s_kin[3][3];
	for(l=0;l<3;l++)
	{
		for(m=0;m<3;m++)
		{
			s_kin[l][m]=0.0;
		}
	}
		
	tbsp=2.0/sqrt(ModuleBase::PI);
		
	int npwx=0;
	int qtot = 0;
	for(int ik=0; ik<GlobalC::kv.nks; ik++)
	{
		for(int ig=0; ig<GlobalC::kv.ngk[ik]; ig++)
		{
			qtot += GlobalC::kv.ngk[ik];
		}
		if(npwx<GlobalC::kv.ngk[ik]) npwx=GlobalC::kv.ngk[ik];
	}	
	gk[0]= new double[npwx]; 
	gk[1]= new double[npwx];
	gk[2]= new double[npwx];
	double factor=ModuleBase::TWO_PI/GlobalC::ucell.lat0;
	//    if(nks>1){
	//       iunigh.clear();
	//       iunigh.seekg(0,ios::beg);
	//    }//go back to the beginning of file

	for(ik=0;ik<GlobalC::kv.nks;ik++)
	{

		npw = GlobalC::kv.ngk[ik];
	//       if(GlobalC::kv.nks>1){
	//          iunigk>>igk;
	//          get_buffer(evc,nwordwfc,iunwfc,ik);
	//       }
		for(i=0;i<npw;i++)
		{
			gk[0][i]=(GlobalC::kv.kvec_c[ik].x+GlobalC::pw.gcar[GlobalC::wf.igk(ik, i)].x)*factor;
			gk[1][i]=(GlobalC::kv.kvec_c[ik].y+GlobalC::pw.gcar[GlobalC::wf.igk(ik, i)].y)*factor;
			gk[2][i]=(GlobalC::kv.kvec_c[ik].z+GlobalC::pw.gcar[GlobalC::wf.igk(ik, i)].z)*factor;
		}

		//kinetic contribution

		for(l=0;l<3;l++)
		{
			for(m=0;m<l+1;m++)
			{
				for(ibnd=0;ibnd<stowf.nchip[ik];ibnd++)
				{
					for(i=0;i<npw;i++)
					{
						s_kin[l][m] += GlobalC::kv.wk[ik]*gk[l][i]*gk[m][i]
						*(conj(stowf.shchi[ik](ibnd, i))*stowf.shchi[ik](ibnd, i)).real();
					}
				}
			}
		}
		   
	}
		

	for(l=0;l<3;l++)
	{
		for(m=0;m<l;m++)
		{
			s_kin[m][l]=s_kin[l][m];
		}
	}

	if(INPUT.gamma_only)
	{
		for(l=0;l<3;l++)
		{
			for(m=0;m<3;m++)
			{
				s_kin[l][m] *= 2.0*ModuleBase::e2/GlobalC::ucell.omega;
			}
		}
	}
	else 
	{
		for(l=0;l<3;l++)
		{
			for(m=0;m<3;m++)
			{
				s_kin[l][m] *= ModuleBase::e2/GlobalC::ucell.omega;
			}
		}
	}

	for(l=0;l<3;l++)
	{
		for(m=0;m<3;m++)
		{
			Parallel_Reduce::reduce_double_all( s_kin[l][m] );
		}
	}


	for(l=0;l<3;l++)
	{
		for(m=0;m<3;m++)
		{
			sigma(l,m) = s_kin[l][m];
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

void Sto_Stress_PW::sto_stress_nl(ModuleBase::matrix& sigma, Stochastic_WF& stowf){
	ModuleBase::TITLE("Sto_Stress_Func","stres_nl");
	ModuleBase::timer::tick("Sto_Stress_Func","stres_nl");
	
	const int nkb = GlobalC::ppcell.nkb;
	if(nkb == 0) return;

	double sigmanlc[3][3];
	for(int l=0;l<3;l++)
	{
		for(int m=0;m<3;m++)
		{
			sigmanlc[l][m]=0.0;
		}
	}
	int* nchip = stowf.nchip;
	// dbecp: conj( -iG * <Beta(nkb,npw)|psi(nbnd,npw)> )
	

	// vkb1: |Beta(nkb,npw)><Beta(nkb,npw)|psi(nbnd,npw)>
	ModuleBase::ComplexMatrix vkb1( nkb, GlobalC::wf.npwx );
	ModuleBase::ComplexMatrix vkb0[3];
	for(int i=0;i<3;i++){
		vkb0[i].create(nkb, GlobalC::wf.npwx);
	}
	ModuleBase::ComplexMatrix vkb2( nkb, GlobalC::wf.npwx );
    for (int ik = 0;ik < GlobalC::kv.nks;ik++)
    {
		ModuleBase::ComplexMatrix dbecp( nkb, nchip[ik]);
		ModuleBase::ComplexMatrix becp( nkb, nchip[ik]);
		for(int i=0;i<3;i++){
			vkb0[i].zero_out();
		}
		vkb2.zero_out();      
		  
		if (GlobalV::NSPIN==2) GlobalV::CURRENT_SPIN = GlobalC::kv.isk[ik];
		GlobalC::wf.npw = GlobalC::kv.ngk[ik];
		// generate vkb
		if (GlobalC::ppcell.nkb > 0)
		{
			GlobalC::ppcell.getvnl(ik);
		}

		// get becp according to wave functions and vkb
		// important here ! becp must set zero!!
		// vkb: Beta(nkb,npw)
		// becp(nkb,nbnd): <Beta(nkb,npw)|psi(nbnd,npw)>
		becp.zero_out();
		for (int ib=0; ib<nchip[ik]; ib++)
		{
			for (int i=0;i<nkb;i++)
			{
				for (int ig=0; ig<GlobalC::wf.npw; ig++)
				{
					becp(i,ib) += stowf.shchi[ik](ib,ig)* conj( GlobalC::ppcell.vkb(i,ig) );
				}
			}
		}
		Parallel_Reduce::reduce_complex_double_pool( becp.c, becp.size);
			
		for(int i=0;i<3;i++) {
			get_dvnl1(vkb0[i],ik,i);
		}
				
		get_dvnl2(vkb2,ik);

		ModuleBase::Vector3<double> qvec;
		double qvec0[3];
				
		for (int ipol = 0; ipol<3; ipol++)
		{
			for(int jpol = 0; jpol < ipol+1; jpol++)
			{
				dbecp.zero_out();    
				vkb1.zero_out();
				for (int i = 0;i < nkb;i++)
				{
					for (int ig=0; ig<GlobalC::wf.npw; ig++)  
					{
						qvec = GlobalC::wf.get_1qvec_cartesian(ik,ig) ;
						qvec0[0] = qvec.x;
						qvec0[1] = qvec.y;
						qvec0[2] = qvec.z;
					  
						vkb1(i, ig) += 0.5 * qvec0[ipol] * vkb0[jpol](i,ig)
								   + 0.5 * qvec0[jpol] * vkb0[ipol](i,ig) ;
					}//end i	  
				}//end nkb
 

				for (int ib=0; ib<nchip[ik]; ib++)
				{
					for (int i=0; i<nkb; i++)
					{
						for (int ig=0; ig<GlobalC::wf.npw; ig++)
						{
							//first term
							dbecp(i,ib) = dbecp(i,ib) - 2.0 * stowf.shchi[ik](ib,ig) * conj( vkb1(i,ig) ) ;
							//second termi
							if(ipol == jpol)
								 dbecp(i,ib) += -1.0 * stowf.shchi[ik](ib,ig)* conj( GlobalC::ppcell.vkb(i,ig) );
							//third term
							qvec = GlobalC::wf.get_1qvec_cartesian(ik,ig);
							qvec0[0] = qvec.x;
							qvec0[1] = qvec.y;
							qvec0[2] = qvec.z;
							double qm1; 
							if(qvec.norm() > 1e-8) qm1 = 1.0 / qvec.norm();
							else qm1 = 0;
							dbecp(i,ib) +=  -2.0 * stowf.shchi[ik](ib,ig) * conj(vkb2(i,ig)) * qvec0[ipol] * qvec0[jpol] * qm1 * GlobalC::ucell.tpiba;
						}//end ig
					}//end i
				}//end ib

//              don't need to reduce here, keep dbecp different in each processor,
//              and at last sum up all the forces.
//              Parallel_Reduce::reduce_complex_double_pool( dbecp.ptr, dbecp.ndata);

//              double *cf = new double[GlobalC::ucell.nat*3];
//              ZEROS(cf, GlobalC::ucell.nat);
				for (int ib=0; ib<nchip[ik]; ib++)
				{
					double fac = 1.0 * GlobalC::kv.wk[ik];
					int iat = 0;
					int sum = 0;
					for (int it=0; it<GlobalC::ucell.ntype; it++)
					{
						const int Nprojs = GlobalC::ucell.atoms[it].nh;
						for (int ia=0; ia<GlobalC::ucell.atoms[it].na; ia++)
						{
							for (int ip=0; ip<Nprojs; ip++)
							{
								double ps = GlobalC::ppcell.deeq(GlobalV::CURRENT_SPIN, iat, ip, ip) ;
								const int inkb = sum + ip;
								//out<<"\n ps = "<<ps;

							 
								const double dbb = ( conj( dbecp( inkb, ib) ) * becp( inkb, ib) ).real();
								sigmanlc[ipol][ jpol] -= ps * fac * dbb;
							 
							}//end ip
							++iat;        
							sum+=Nprojs;
						}//ia
					} //end it
				} //end band
			}//end jpol
		}//end ipol
	}// end ik

	// sum up forcenl from all processors
	for(int l=0;l<3;l++){
		for(int m=0;m<3;m++){
			if(m>l) sigmanlc[l][m] = sigmanlc[m][l];
				Parallel_Reduce::reduce_double_all( sigmanlc[l][m] );
		}
	}

//        Parallel_Reduce::reduce_double_all(sigmanl.c, sigmanl.nr * sigmanl.nc);
        
	for (int ipol = 0; ipol<3; ipol++)
	{
		for(int jpol = 0; jpol < 3; jpol++)
		{
			sigmanlc[ipol][jpol] *= 1.0 / GlobalC::ucell.omega;
		}
	}
	
	for (int ipol = 0; ipol<3; ipol++)
	{
		for(int jpol = 0; jpol < 3; jpol++)
		{
			sigma(ipol,jpol) = sigmanlc[ipol][jpol] ;
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
