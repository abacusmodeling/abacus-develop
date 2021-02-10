#include "print_info.h"
#include "../src_global/global_variable.h"
#include "../src_pw/global.h"

Print_Info::Print_Info(){}

Print_Info::~Print_Info(){}


void Print_Info::screen_output(void)
{
	
	// the following printing information should be moved to somewhere else -- mohan 2021-01-30
    if(CALCULATION=="scf" || CALCULATION=="relax" || CALCULATION=="cell-relax" || CALCULATION=="nscf"
	        || CALCULATION=="istate" || CALCULATION=="ienvelope"
	        || CALCULATION=="md") //pengfei add 2014-10-13
	{
		//LM.divide_HS_in_frag(); //move it above 2015-09-06, xiaohui
		cout << " ---------------------------------------------------------" << endl;
		if(CALCULATION=="scf")
		{
			cout << " This calculation is self-consistent" << endl;
		}
		else if(CALCULATION=="test")
		{
			cout << " This calculation is for test" << endl;
		}
		if(CALCULATION=="relax") //add 4 lines 2015-09-06, xiaohui
		{
            //cout << " This calculation is structure relaxation" << endl;
            cout << " This calculation is ion relaxation" << endl;
		}
        if(CALCULATION=="cell-relax")
        {
            cout << " This calculation is cell relaxation" << endl;
        }
		if(CALCULATION=="md") //add 4 lines 2015-09-06, xiaohui
		{
			cout << " This calculation is molecular dynamics" << endl;

			//xiaohui add 2015-09-15
			cout << " ---------------------------------------------------------" << endl;

			if(INPUT.md_mdtype ==1 || INPUT.md_mdtype==2)
			{
				cout << " ENSEMBLE                 : " << "NVT" << endl;
			}
			else if(INPUT.md_mdtype==0)
			{
				cout << " ENSEMBLE                 : " << "NVE" << endl;
			}
			cout << " Qmass for NVT(a.u.)      : " << INPUT.md_qmass/6.02/9.109*1e5 << endl;
			cout << " Time interval(fs)        : " << INPUT.md_dt << endl;
		}
		cout << " ---------------------------------------------------------" << endl;


		// TITLE
		cout << " " << setw(8) << "SPIN"
		     << setw(16) << "KPOINTS"
		     << setw(12) << "PROCESSORS";

		//if(LOCAL_BASIS==4) xiaohui modify 2013-09-01
		if(BASIS_TYPE=="lcao" || BASIS_TYPE=="lcao_in_pw") //xiaohui add 2013-09-01
		{
			cout << setw(12) << "NBASE";
			cout << setw(12) << "VNA";
		}

		cout << endl;



		// data
		cout << " " << setw(8) << NSPIN;

		if(GAMMA_ONLY_LOCAL)
		{
			if(COLOUR && MY_RANK==0)
			{
				// red
				//printf( "\e[31m%-16s\e[0m", "Gamma");
				printf( "[31m%-16s[0m", "Gamma");
			}
			else
			{
				cout << setw(16) << "Gamma";
			}
		}
		else
		{
			if(COLOUR && MY_RANK==0)
			{
				// zi
				//printf( "\e[35m%-16d\e[0m", kv.nkstot);
				printf( "[35m%-16d[0m", kv.nkstot);
			}
			else
			{
				cout << setw(16) << kv.nkstot;
			}
		}

		cout << setw(12) << NPROC;

		//if(LOCAL_BASIS==4) xiaohui modify 2013-09-01
		if(BASIS_TYPE=="lcao" || BASIS_TYPE=="lcao_in_pw") //xiaohui add 2013-09-01
		{
			cout << setw(12) << NLOCAL;
			// print VNA: no, should delete in future -- mohan 2021-02-09
			cout << setw(12) << "No";
		}

		cout << endl;




		cout << " ---------------------------------------------------------" << endl;
		//if(LOCAL_BASIS==4 && LINEAR_SCALING==1) xiaohui modify 2013-09-01
		if(BASIS_TYPE=="lcao") //xiaohui add 2013-09-01
		{
			if(COLOUR && MY_RANK==0)
			{
				string a = "Use Systematically Improvable Atomic bases";
				//printf( " \e[36m%-45s\e[0m\n", a.c_str());
				printf( " [36m%-45s[0m\n", a.c_str());
			}
			else
			{
				cout << " Use Systematically Improvable Atomic bases" << endl;
			}
		}
		//else if(LOCAL_BASIS==4 && LINEAR_SCALING==0) xiaohui modify 2013-09-01
		else if(BASIS_TYPE=="lcao_in_pw") //xiaohui add 2013-09-01
		{
			//cout << " Expand Systematically Improvable Atomic bases into plane waves" << endl;
			cout << " Expand Atomic bases into plane waves" << endl;
		}
		//else if(LOCAL_BASIS==0 && LINEAR_SCALING==0) xiaohui modify 2013-09-01
		else if(BASIS_TYPE=="pw") //xiaohui add 2013-09-01
		{
			cout << " Use plane wave basis" << endl;
		}
		cout << " ---------------------------------------------------------" << endl;



		//----------------------------------
		// second part
		//----------------------------------

		cout << " " << setw(8) << "ELEMENT";

		//if(LOCAL_BASIS==4) xiaohui modify 2013-09-01
		if(BASIS_TYPE=="lcao" || BASIS_TYPE=="lcao_in_pw") //xiaohui add 2013-09-01
		{
			cout << setw(16) << "ORBITALS";
			cout << setw(12) << "NBASE";
			//cout << setw(12) << "NATOM"; //move it below 2015-09-06, xiaohui
		}
		cout << setw(12) << "NATOM"; //add 2015-09-06, xiaohui

		cout << setw(12) << "XC";
		cout << endl;



		for(int it=0; it<ucell.ntype; ++it)
		{
			if(COLOUR && MY_RANK==0)
			{
				//printf( " \e[36m%-8s\e[0m", ucell.atoms[it].label.c_str());
				printf( " [36m%-8s[0m", ucell.atoms[it].label.c_str());
			}
			else
			{
				cout << " " << setw(8) << ucell.atoms[it].label;
			}

			//if(LOCAL_BASIS==4) xiaohui modify 2013-09-01
			if(BASIS_TYPE=="lcao" || BASIS_TYPE=="lcao_in_pw") //xiaohui add 2013-09-01
			{
				stringstream orb;

				int norb = 0;
				/*for(int L=0; L<=ORB.Phi[it].getLmax(); ++L)
				{
					norb += (2*L+1)*ORB.Phi[it].getNchi(L);
					orb << ORB.Phi[it].getNchi(L);
					if(L==0) orb << "s";
					else if(L==1) orb << "p";
					else if(L==2) orb << "d";
					else if(L==3) orb << "f";
					else if(L==4) orb << "g";
					else if(L==5) orb << "h";
					else if(L==6) orb << "i";
				}
				orb << "-" << ORB.Phi[it].getRcut() << "au";*/


				for(int L=0; L<=ucell.atoms[it].nwl; ++L)        // pengfei Li 16-2-29
				{
					norb += (2*L+1)* ucell.atoms[it].l_nchi[L];
					orb << ucell.atoms[it].l_nchi[L];
					if(L==0) orb << "s";
					else if(L==1) orb << "p";
					else if(L==2) orb << "d";
					else if(L==3) orb << "f";
					else if(L==4) orb << "g";
					else if(L==5) orb << "h";
					else if(L==6) orb << "i";
				}
				orb << "-" << ucell.atoms[it].Rcut << "au";

				if(COLOUR && MY_RANK==0)
				{
					//printf( "\e[36m%-16s\e[0m", orb.str().c_str());
					//printf( "\e[36m%-12d\e[0m", norb);
					printf( "[36m%-16s[0m", orb.str().c_str());
					printf( "[36m%-12d[0m", norb);
				}
				else
				{
					cout << setw(16) << orb.str();
					cout << setw(12) << norb;
				}
			}


			cout << setw(12) << ucell.atoms[it].na;

//				if(ucell.atoms[it].dft[1]=="PZ")    // pengfei Li added 2015-1-31 cancelled by zws
//				{
//					//cout << setw(12) << "PZ-LDA";
//
//				}
//				else
//				{
//					//cout << setw(12) << ucell.atoms[it].dft[0];
//                                        cout << setw(12) << "PBE";
//				}
			xcf.ostreamdft(cout); // zws add 20150108
			//cout << " ( "  << setw(3) << xcf.iexch << setw(3) << xcf.icorr << setw(3) << xcf.igcx << setw(3) << xcf.igcc << ")";
			cout << endl;
		}

		cout << " ---------------------------------------------------------" << endl;
		cout << " Initial plane wave basis and FFT box" << endl;
		cout << " ---------------------------------------------------------" << endl;

//			cout << " GRID_SPEED           : " << GRID_SPEED << endl;
	}

	return;
}
