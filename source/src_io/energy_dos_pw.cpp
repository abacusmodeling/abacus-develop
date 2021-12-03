#include "dos.h"
#include "../src_pw/global.h"
#include "../src_pw/energy.h"

void energy::perform_dos_pw(void)
{
	ModuleBase::TITLE("energy","perform_dos_pw");

	if(out_dos !=0 || out_band !=0)
    {
        GlobalV::ofs_running << "\n\n\n\n";
        GlobalV::ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
        GlobalV::ofs_running << " |                                                                    |" << std::endl;
        GlobalV::ofs_running << " | Post-processing of data:                                           |" << std::endl;
        GlobalV::ofs_running << " | DOS (density of states) and bands will be output here.             |" << std::endl;
        GlobalV::ofs_running << " | If atomic orbitals are used, Mulliken charge analysis can be done. |" << std::endl;
        GlobalV::ofs_running << " | Also the .bxsf file containing fermi surface information can be    |" << std::endl;
        GlobalV::ofs_running << " | done here.                                                         |" << std::endl;
        GlobalV::ofs_running << " |                                                                    |" << std::endl;
        GlobalV::ofs_running << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
        GlobalV::ofs_running << "\n\n\n\n";
    }

	//qianrui modify 2020-10-18
	if(GlobalV::CALCULATION=="scf" || GlobalV::CALCULATION=="md" || GlobalV::CALCULATION=="relax")
	{
		std::stringstream ss;
		ss << GlobalV::global_out_dir << "istate.info" ;
		if(GlobalV::MY_RANK==0)
		{
			std::ofstream ofsi( ss.str().c_str() ); // clear istate.info
			ofsi.close();
		}
#ifdef __MPI
		for(int ip=0; ip<GlobalV::NPOOL; ip++)
		{
			MPI_Barrier(MPI_COMM_WORLD);
			if( GlobalV::MY_POOL == ip )
			{
				if( GlobalV::RANK_IN_POOL != 0 ) continue;
#endif
				std::ofstream ofsi2( ss.str().c_str(), ios::app );
				if(GlobalV::NSPIN == 1||GlobalV::NSPIN == 4)
				{
					for (int ik = 0;ik < GlobalC::kv.nks;ik++)
					{
						ofsi2<<"BAND"
						<<std::setw(25)<<"Energy(ev)"
						<<std::setw(25)<<"Occupation"
#ifdef __MPI
						<<std::setw(25)<<"Kpoint = "<<GlobalC::Pkpoints.startk_pool[ip]+ik+1
#else
						<<std::setw(25)<<"Kpoint = "<<ik+1
#endif
						<<std::setw(25)<<"("<<GlobalC::kv.kvec_d[ik].x<<" "<<GlobalC::kv.kvec_d[ik].y<<" "<<GlobalC::kv.kvec_d[ik].z<<")"<<std::endl;
						for(int ib=0;ib<GlobalV::NBANDS;ib++)
						{
							ofsi2<<std::setw(6)<<ib+1<<std::setw(25)<<GlobalC::wf.ekb[ik][ib]* ModuleBase::Ry_to_eV<<std::setw(25)<<GlobalC::wf.wg(ik,ib)<<std::endl;
						}
						ofsi2 <<std::endl;
						ofsi2 <<std::endl;
					}
				}
				else
				{
					for (int ik = 0;ik < GlobalC::kv.nks/2;ik++)
					{
						ofsi2<<"BAND"
						<<std::setw(25)<<"Spin up Energy(ev)"
						<<std::setw(25)<<"Occupation"
						<<std::setw(25)<<"Spin down Energy(ev)"
						<<std::setw(25)<<"Occupation"
#ifdef __MPI
						<<std::setw(25)<<"Kpoint = "<<GlobalC::Pkpoints.startk_pool[ip]+ik+1
#else
						<<std::setw(25)<<"Kpoint = "<<ik+1
#endif
						<<std::setw(25)<<"("<<GlobalC::kv.kvec_d[ik].x<<" "<<GlobalC::kv.kvec_d[ik].y<<" "<<GlobalC::kv.kvec_d[ik].z<<")"<<std::endl;

						for(int ib=0;ib<GlobalV::NBANDS;ib++)
						{
							ofsi2<<std::setw(6)<<ib+1
							<<std::setw(25)<<GlobalC::wf.ekb[ik][ib]* ModuleBase::Ry_to_eV
							<<std::setw(25)<<GlobalC::wf.wg(ik,ib)
							<<std::setw(25)<<GlobalC::wf.ekb[(ik+GlobalC::kv.nks/2)][ib]* ModuleBase::Ry_to_eV
							<<std::setw(25)<<GlobalC::wf.wg(ik+GlobalC::kv.nks/2,ib)<<std::endl;
						}
						ofsi2 <<std::endl;
						ofsi2 <<std::endl;

					}
				}

				ofsi2.close();
#ifdef __MPI
			}
		}
#endif
	}
	
	int nspin0=1;
	if(GlobalV::NSPIN==2) nspin0=2;

	if(this->out_dos)
	{
//find energy range
		double emax = GlobalC::wf.ekb[0][0];
		double emin = GlobalC::wf.ekb[0][0];
		for(int ik=0; ik<GlobalC::kv.nks; ++ik)
		{
			for(int ib=0; ib<GlobalV::NBANDS; ++ib)
			{
				emax = std::max( emax, GlobalC::wf.ekb[ik][ib] );
				emin = std::min( emin, GlobalC::wf.ekb[ik][ib] );
			}
		}

#ifdef __MPI
		Parallel_Reduce::gather_max_double_all(emax);
		Parallel_Reduce::gather_min_double_all(emin);
#endif

		emax *= ModuleBase::Ry_to_eV;
		emin *= ModuleBase::Ry_to_eV;

//scale up a little bit so the end peaks are displaced better
		const double scl=this->dos_scale;
		double delta=(emax-emin)*scl;
		std::cout << scl;
		emax=emax+delta/2.0;
		emin=emin-delta/2.0;

				ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"minimal energy is (eV)", emin);
				ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"maximal energy is (eV)", emax);
		// 		atom_arrange::set_sr_NL();
		//		atom_arrange::search( GlobalV::SEARCH_RADIUS );//qifeng-2019-01-21
		
//determine #. energy points	
		const double de_ev = this->dos_edelta_ev;
		std::cout << de_ev;

		const int npoints = static_cast<int>(std::floor ( ( emax - emin ) / de_ev ));
		const int np=npoints;

	 	for(int is=0; is<nspin0; ++is)
	 	{
//DOS_ispin contains not smoothed dos
			 std::stringstream ss;
			 ss << GlobalV::global_out_dir << "DOS" << is+1;

			 Dos::calculate_dos(
					 is,
					 GlobalC::kv.isk,
					 ss.str(), 
					 this->dos_edelta_ev, 
					 emax, 
					 emin, 
					 GlobalC::kv.nks, GlobalC::kv.nkstot, GlobalC::kv.wk, GlobalC::wf.wg, GlobalV::NBANDS, GlobalC::wf.ekb );
			 std::ifstream in(ss.str().c_str());
			 if(!in)
			 {
				       //std::cout<<"\n Can't find file : "<< name << std::endl;
				       //return 0;
			 }

			 //----------------------------------------------------------
			 // FOUND LOCAL VARIABLES :
			 // NAME : number(number of DOS points)
			 // NAME : nk(number of k point used)
			 // NAME : energy(energy range,from emin_ev to emax_ev)
			 // NAME : dos(old,count k points in the energy range)
			 // NAME : dos2(new,count k points in the energy range)
			 //----------------------------------------------------------
			 int number=0;
			 int nk=0;
			 in >> number;
			 in >> nk;
			 double *energy = new double[number];
			 double *dos = new double[number];
			 double *dos2 = new double[number];
			 for(int i=0 ;i<number; i++)
			 {
				 energy[i] = 0.0;
				 dos[i] = 0.0;
				 dos2[i] =0.0;
			 }

			 for(int i=0;i<number;i++)
			 {
				 in >> energy[i] >> dos[i];
			 }
			 if(!in.eof())
			 {
				 //std::cout<<"\n Read Over!"<<std::endl;
			 }
			 in.close();

//now use Gaussian smearing to smooth the dos and write to DOS_is_smearing

			 //----------------------------------------------------------
			 // EXPLAIN : b is an empirical value.
			 // please DIY b!!
			 //----------------------------------------------------------

			 //double b = INPUT.b_coef;
			 double b = sqrt(2.0)*bcoeff;
			 for(int i=0;i<number;i++)
			 {
				 double Gauss=0.0;
	
				 for(int j=0;j<number;j++)
				 {
					 double de = energy[j] - energy[i];
					 double de2 = de * de;
					 //----------------------------------------------------------
					 // EXPLAIN : if en
					 //----------------------------------------------------------
					 Gauss = exp(-de2/b/b)/sqrt(3.1415926)/b;
					 dos2[j] += dos[i]*Gauss;
				 }
			 }

		 //----------------------------------------------------------
		 // EXPLAIN : output DOS2.txt
		 //----------------------------------------------------------
		 std::stringstream sss;
		 sss << GlobalV::global_out_dir << "DOS" << is+1 << "_smearing" << ".dat" ;
		 std::ofstream out(sss.str().c_str());
		 double sum2=0.0;
		 for(int i=0;i<number;i++)
		 {
			 sum2 += dos2[i];
			 //            if(dos2[i]<1e-5)
			 //            {
			 //                    dos2[i] = 0.00;
			 //            }
			 out <<std::setw(20)<<energy[i]
				 <<std::setw(20)<<dos2[i]
				 <<std::setw(20)<<sum2<<"\n";
		 }
		 out.close();

		 //----------------------------------------------------------
		 // DELETE
		 //----------------------------------------------------------
		 delete[] dos;
		 delete[] dos2;
		 delete[] energy;

		 //std::cout<<" broden spectrum over, success : ) "<<std::endl;

	 }

	}//out_dos=1
	if(this->out_band) //pengfei 2014-10-13
	{
		int nks=0;
		if(nspin0==1) 
		{
			nks = GlobalC::kv.nkstot;
		}
		else if(nspin0==2) 
		{
			nks = GlobalC::kv.nkstot/2;
		}

		for(int is=0; is<nspin0; is++)
		{
			std::stringstream ss2;
			ss2 << GlobalV::global_out_dir << "BANDS_" << is+1 << ".dat";
			GlobalV::ofs_running << "\n Output bands in file: " << ss2.str() << std::endl;
			Dos::nscf_band(is, ss2.str(), nks, GlobalV::NBANDS, this->ef*0, GlobalC::wf.ekb);
		}

	}
}

