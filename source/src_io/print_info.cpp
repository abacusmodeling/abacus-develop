#include "print_info.h"
#include "../module_base/global_variable.h"
#include "../src_pw/global.h"

Print_Info::Print_Info(){}

Print_Info::~Print_Info(){}


void Print_Info::setup_parameters(void)
{
	TITLE("Print_Info","setup_parameters");
	
    if(GlobalV::CALCULATION=="scf" || GlobalV::CALCULATION=="relax" || GlobalV::CALCULATION=="cell-relax" || GlobalV::CALCULATION=="nscf"
	        || GlobalV::CALCULATION=="istate" || GlobalV::CALCULATION=="ienvelope" || GlobalV::CALCULATION=="md")
	{
		std::cout << " ---------------------------------------------------------" << std::endl;
		if(GlobalV::CALCULATION=="scf")
		{
			std::cout << " Self-consistent calculations for electrons" << std::endl;
		}
		else if(GlobalV::CALCULATION=="test")
		{
			std::cout << " Test run" << std::endl;
		}
		if(GlobalV::CALCULATION=="relax")
		{
            std::cout << " Ion relaxation calculations" << std::endl;
		}
        if(GlobalV::CALCULATION=="cell-relax")
        {
            std::cout << " Cell relaxation calculations" << std::endl;
        }
		if(GlobalV::CALCULATION=="md")
		{
			std::cout << " Molecular Dynamics simulations" << std::endl;

			std::cout << " ---------------------------------------------------------" << std::endl;

			if(INPUT.mdp.mdtype ==1 || INPUT.mdp.mdtype==2)
			{
				std::cout << " ENSEMBLE                 : " << "NVT" << std::endl;
				std::cout << " Qmass for NVT(a.u.)      : " << INPUT.mdp.Qmass/6.02/9.109*1e5 << std::endl;
			}
			else if(INPUT.mdp.mdtype==0)
			{
				std::cout << " ENSEMBLE                 : " << "NVE" << std::endl;
			}
			
			std::cout << " Time interval(fs)        : " << INPUT.mdp.dt << std::endl;
		}
		std::cout << " ---------------------------------------------------------" << std::endl;


		std::cout << " " << setw(8) << "SPIN"
		     << setw(16) << "KPOINTS"
		     << setw(12) << "PROCESSORS";

		if(GlobalV::BASIS_TYPE=="lcao" || GlobalV::BASIS_TYPE=="lcao_in_pw")
		{
			std::cout << setw(12) << "NBASE";
		}

		std::cout << std::endl;



		std::cout << " " << setw(8) << GlobalV::NSPIN;

		if(GlobalV::GAMMA_ONLY_LOCAL)
		{
			if(GlobalV::COLOUR && GlobalV::MY_RANK==0)
			{
				// red
				printf( "\e[31m%-16s\e[0m", "Gamma");
				//printf( "[31m%-16s[0m", "Gamma");
			}
			else
			{
				std::cout << setw(16) << "Gamma";
			}
		}
		else
		{
			if(GlobalV::COLOUR && GlobalV::MY_RANK==0)
			{
				// zi
				printf( "\e[35m%-16d\e[0m", GlobalC::kv.nkstot);
				//printf( "[35m%-16d[0m", kv.nkstot);
			}
			else
			{
				std::cout << setw(16) << GlobalC::kv.nkstot;
			}
		}

		std::cout << setw(12) << GlobalV::NPROC;

		if(GlobalV::BASIS_TYPE=="lcao" || GlobalV::BASIS_TYPE=="lcao_in_pw")
		{
			std::cout << setw(12) << GlobalV::NLOCAL;
		}

		std::cout << std::endl;




		std::cout << " ---------------------------------------------------------" << std::endl;
		if(GlobalV::BASIS_TYPE=="lcao") 
		{
			if(GlobalV::COLOUR && GlobalV::MY_RANK==0)
			{
				std::string a = "Use Systematically Improvable Atomic bases";
				printf( " \e[36m%-45s\e[0m\n", a.c_str());
				//printf( " [36m%-45s[0m\n", a.c_str());
			}
			else
			{
				std::cout << " Use Systematically Improvable Atomic bases" << std::endl;
			}
		}
		else if(GlobalV::BASIS_TYPE=="lcao_in_pw")
		{
			std::cout << " Expand Atomic bases into plane waves" << std::endl;
		}
		else if(GlobalV::BASIS_TYPE=="pw")
		{
			std::cout << " Use plane wave basis" << std::endl;
		}
		std::cout << " ---------------------------------------------------------" << std::endl;



		//----------------------------------
		// second part
		//----------------------------------

		std::cout << " " << setw(8) << "ELEMENT";

		if(GlobalV::BASIS_TYPE=="lcao" || GlobalV::BASIS_TYPE=="lcao_in_pw")
		{
			std::cout << setw(16) << "ORBITALS";
			std::cout << setw(12) << "NBASE";
		}
		std::cout << setw(12) << "NATOM";

		std::cout << setw(12) << "XC";
		std::cout << std::endl;



		for(int it=0; it<GlobalC::ucell.ntype; ++it)
		{
			if(GlobalV::COLOUR && GlobalV::MY_RANK==0)
			{
				printf( "\e[36m%-8s\e[0m", GlobalC::ucell.atoms[it].label.c_str());
			}
			else
			{
				std::cout << " " << setw(8) << GlobalC::ucell.atoms[it].label;
			}

			if(GlobalV::BASIS_TYPE=="lcao" || GlobalV::BASIS_TYPE=="lcao_in_pw")
			{
				std::stringstream orb;

				int norb = 0;
				/*for(int L=0; L<=GlobalC::ORB.Phi[it].getLmax(); ++L)
				{
					norb += (2*L+1)*GlobalC::ORB.Phi[it].getNchi(L);
					orb << GlobalC::ORB.Phi[it].getNchi(L);
					if(L==0) orb << "s";
					else if(L==1) orb << "p";
					else if(L==2) orb << "d";
					else if(L==3) orb << "f";
					else if(L==4) orb << "g";
					else if(L==5) orb << "h";
					else if(L==6) orb << "i";
				}
				orb << "-" << GlobalC::ORB.Phi[it].getRcut() << "au";*/


				for(int L=0; L<=GlobalC::ucell.atoms[it].nwl; ++L)        // pengfei Li 16-2-29
				{
					norb += (2*L+1)* GlobalC::ucell.atoms[it].l_nchi[L];
					orb << GlobalC::ucell.atoms[it].l_nchi[L];
					if(L==0) orb << "s";
					else if(L==1) orb << "p";
					else if(L==2) orb << "d";
					else if(L==3) orb << "f";
					else if(L==4) orb << "g";
					else if(L==5) orb << "h";
					else if(L==6) orb << "i";
				}
				orb << "-" << GlobalC::ucell.atoms[it].Rcut << "au";

				if(GlobalV::COLOUR && GlobalV::MY_RANK==0)
				{
					printf( "\e[36m%-16s\e[0m", orb.str().c_str());
					printf( "\e[36m%-12d\e[0m", norb);
					//printf( "[36m%-16s[0m", orb.str().c_str());
					//printf( "[36m%-12d[0m", norb);
				}
				else
				{
					std::cout << setw(16) << orb.str();
					std::cout << setw(12) << norb;
				}
			}


			std::cout << setw(12) << GlobalC::ucell.atoms[it].na;

//				if(GlobalC::ucell.atoms[it].dft[1]=="PZ")    // pengfei Li added 2015-1-31 cancelled by zws
//				{
//					//std::cout << setw(12) << "PZ-LDA";
//
//				}
//				else
//				{
//					//std::cout << setw(12) << GlobalC::ucell.atoms[it].dft[0];
//                                        std::cout << setw(12) << "PBE";
//				}
			GlobalC::xcf.ostreamdft(std::cout); // zws add 20150108
			//std::cout << " ( "  << setw(3) << GlobalC::xcf.iexch << setw(3) << GlobalC::xcf.icorr << setw(3) << GlobalC::xcf.igcx << setw(3) << GlobalC::xcf.igcc << ")";
			std::cout << std::endl;
		}

		std::cout << " ---------------------------------------------------------" << std::endl;
		std::cout << " Initial plane wave basis and FFT box" << std::endl;
		std::cout << " ---------------------------------------------------------" << std::endl;

	}

	return;
}
