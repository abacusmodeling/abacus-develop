#include "dc_driv.h"

#ifdef __EPM
#include "../src_epm/run_epm.h"
#else
#include "run_frag.h"
#endif

#include "input.h"
#include "input_conv.h"
#include "src_lcao/global_fp.h"
#include "src_pw/global.h"

DC_Driv::DC_Driv()
{}

DC_Driv::~DC_Driv()
{}

void DC_Driv::init()
{
	TITLE("DC_Driv","init");

	time_t  time_start = std::time(NULL);

	timer::start();

	// read the parameters.
	this->reading();

#ifdef __FP
	// divide the system into fragments.
	this->divide_frag();

	// setup the information for each fragment.
	this->setup_frag();
#endif

	// solve each fragment.
	this->solve_eachf();

	time_t	time_finish= std::time(NULL);

	cout << "\n START  Time  : " << ctime(&time_start);
	cout << " FINISH Time  : " << ctime(&time_finish);
	cout << " TOTAL  Time  : " << difftime(time_finish, time_start) << endl;
	cout << " SEE INFORMATION IN : "<<global_out_dir<<endl;

	ofs_running << "\n Start  Time  : " << ctime(&time_start);
	ofs_running << " Finish Time  : " << ctime(&time_finish);

	double total_time = difftime(time_finish, time_start);
	int hour = total_time / 3600;
	int mins = ( total_time - 3600 * hour ) / 60;
	int secs = total_time - 3600 * hour - 60 * mins ;
	ofs_running << " Total  Time  : " << hour << " h "
	            << mins << " mins "
	            << secs << " secs "<< endl;

	INPUT.close_log();

	return;
}

void DC_Driv::reading(void)
{
	timer::tick("DC_Driv","reading",'A');
	INPUT.Init( global_in_card );
	Input_Conv::Convert();
#ifdef __EPM
	Input_Conv::Convert_EPM();
#else
	Input_Conv::Convert_FP();
#endif

#ifdef __FP

	Parallel_Global::split_diag_world(DIAGO_PROC);
	Parallel_Global::split_grid_world(DIAGO_PROC);
	OUT(ofs_running,"DRANK",DRANK+1);
	OUT(ofs_running,"DSIZE",DSIZE);
	OUT(ofs_running,"DCOLOR",DCOLOR+1);
	OUT(ofs_running,"GRANK",GRANK+1);
	OUT(ofs_running,"GSIZE",GSIZE);

	Run_Frag::frag_init();
	//if(LOCAL_BASIS==4 && LINEAR_SCALING) xiaohui modify 2013-09-01

	if(BASIS_TYPE=="lcao") //xiaohui add 2013-09-01
	{

		// read orbital information.
		// init overlap matrix table.
		// init kinetical matrix element table.
		// init non-local pseudopotential matrix element table.
		hm.hon.set_orb_tables();

		LM.divide_HS_in_frag(); //add 2015-09-06, xiaohui
	} //add 2015-09-06, xiaohui



    if(CALCULATION=="scf" || CALCULATION=="relax" || CALCULATION=="cell-relax" || CALCULATION=="nscf"
	        || CALCULATION=="istate" || CALCULATION=="ienvelope"
	        || CALCULATION=="md") //pengfei 2014-10-13
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

			if(VNA==0)
				cout << setw(12) << "No";
			else if(VNA>0)
				cout << setw(12) << "Yes";
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


	//} //delete 2015-09-06, xiaohui
	/*
		else if(BASIS_TYPE=="pw") //2015-09-06, xiaohui
		{
			cout << " " << setw(8) << "ELEMENT";
			cout << setw(12) << "NATOM";
			cout << setw(12) << "XC";
			cout << endl;
			for(int it=0; it<ucell.ntype; ++it)
			{
				cout << " " << setw(8) << ucell.atoms[it].label;
				cout << setw(12) << ucell.atoms[it].na;
				if(ucell.atoms[it].dft[1]=="PZ")    // pengfei Li added 2015-1-31
				{
					//cout << " XC FUNCTIONAL      : " << "PZ-LDA" << endl;
					cout << setw(12) << "PZ-LDA";
				}
				else
				{
					//cout << " XC FUNCTIONAL      : " << "PBE" << endl;
					cout << setw(12) << "PZ-LDA";
				}
			}
			cout<<endl;
		}
	*/
#endif



	timer::tick("DC_Driv","reading",'A');
	return;
}

#include "src_pw/cal_test.h"
#include "src_pw/cal_test0.h"
//#include "../src_develop/src_dc/dc_info.h"
void DC_Driv::divide_frag(void)
{
	TITLE("DC_Driv","divide_frag");
	timer::tick("DC_Driv","divide_frag",'A');

	// (1) Init the plane wave.
	pw.gen_pw(ofs_running, ucell, kv);
	DONE(ofs_running,"INIT PLANEWAVE");
	cout << " UNIFORM GRID DIM     : " << pw.nx <<" * " << pw.ny <<" * "<< pw.nz << endl;
	cout << " UNIFORM GRID DIM(BIG): " << pw.nbx <<" * " << pw.nby <<" * "<< pw.nbz << endl;

	//bool test_dc = false;
	//if(test_dc)
	//{
	//	DC_Info::divide_fragments(ucell);
	//}

	// mohan add 2010-10-10, just to test the symmetry of a variety
	// of systems.
	//xiaohui modified 2013-03-23,adding "/*"
	if(CALCULATION == "test")
	{
		Cal_Test::adjacent_atoms();
		//Cal_Test::sparsity();
		Cal_Test::test_memory();
		QUIT();
	}

	// mohan add 2010-09-13
	// init the grid, then the charge
	// on grid can be distributed.
	Pgrid.init(pw.ncx, pw.ncy, pw.ncz, pw.nczp, pw.nrxx, pw.nbz, pw.bz); // mohan add 2010-07-22, update 2011-05-04


	timer::tick("DC_Driv","divide_frag",'A');
	return;
}

void DC_Driv::setup_frag(void)
{
	TITLE("DC_Driv","setup_frag");

}

void DC_Driv::solve_eachf(void)
{
	TITLE("DC_Driv","solve_eachf");
	timer::tick("DC_Driv","solve_eachf",'A');

#ifdef __EPM

	Run_EPM re;
	re.epm_line();

#else

	Run_Frag RF;


#ifdef __FP

	//xiaohui modify 2013-09-01
	//if(LINEAR_SCALING==-1)
	//{
	//    RF.frag_test();
	//}
	//else if (LINEAR_SCALING==0)
	//{
	//    RF.frag_pw_line();
	//}
	//else if (LINEAR_SCALING==1)
	//{
	//    RF.frag_LCAO_line();
	//}
	//else if (LINEAR_SCALING==2)
	//{
	//    RF.frag_linear_scaling_line();
	//}
	//xiaohui add 2013-09-01
	if(BASIS_TYPE=="pw" || BASIS_TYPE=="lcao_in_pw")
	{
		RF.frag_pw_line();
	}
	else if(BASIS_TYPE=="lcao")
	{
		RF.frag_LCAO_line();
	}

#else

	// if use plane wave basis,
	// the plane wave basis is generated in pw_line.
	RF.pw_line();

#endif

#endif

	timer::tick("DC_Driv","solve_eachf",'A');
	timer::finish( ofs_running );

	Memory::print_all( ofs_running ) ;

	return;
}

