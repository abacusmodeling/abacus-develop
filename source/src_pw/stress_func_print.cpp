#include"stress_func.h"

static double output_acc = 1.0e-8;

//print target stress term 
void Stress_Func::print_stress(const string &name, const matrix& f, const bool screen, bool ry)const
{
	ofs_running << " --------------------------- " << name << " ----------------------------" << endl;
	
	
	double fac = 1.0;
	
	if(!ry)
	{
	 //	fac = Ry_to_eV / 0.529177;
	}

	cout << setprecision(5);
	cout << setiosflags(ios::showpos);

	if(screen)
	{
		cout << " ------------------- " << name << " --------------------" << endl;
	
	}

	for (int i=0;i<3;i++)
	{
		ofs_running << setw(15)<< " ";
		if( abs(f(i,0)) >output_acc) ofs_running << setw(15) << f(i,0) * fac;
		else ofs_running << setw(15) << "0";
		if( abs(f(i,1)) >output_acc) ofs_running << setw(15) << f(i,1) * fac;
		else ofs_running << setw(15) << "0";
		if( abs(f(i,2)) >output_acc) ofs_running << setw(15) << f(i,2) * fac;
		else ofs_running << setw(15) << "0";
		ofs_running << endl;

		if(screen)
		{
			if( abs(f(i,0)) >output_acc) cout << setw(15) << f(i,0)*fac;
			else cout << setw(15) << "0";
			if( abs(f(i,1)) >output_acc) cout << setw(15) << f(i,1)*fac;
			else cout << setw(15) << "0";
			if( abs(f(i,2)) >output_acc) cout << setw(15) << f(i,2)*fac;
			else cout << setw(15) << "0";
			cout << endl;
		}	
	}


	cout << resetiosflags(ios::showpos);

    return;
}

//print total stress
void Stress_Func::printstress_total(const matrix& scs, bool ry)
{
// zhengdy update 2016-10-08
	double unit_transform = 1;

	if(!ry)
	{
		unit_transform = RYDBERG_SI / pow(BOHR_RADIUS_SI,3) * 1.0e-8;
	}
//	cout.setf(ios::fixed);


	//ofs_running << setiosflags(ios::right);
 	ofs_running << setprecision(6) << setiosflags(ios::showpos) << setiosflags(ios::fixed) << endl;
	NEW_PART("TOTAL-STRESS (KBAR)");//Ryd/(a.u.)^3
    cout << " ><><><><><><><><><><><><><><><><><><><><><><" << endl;
    cout << " TOTAL-STRESS (KBAR):" << endl;
    cout << " ><><><><><><><><><><><><><><><><><><><><><><" << endl;

//        if(INPUT.stress_set == 1)
//        int TEST_STRESS = 1;

 	if(TEST_STRESS) 
	{
		cout << setiosflags(ios::fixed) << setprecision(6);
		cout << setiosflags(ios::showpos);
		cout << " ------------------- TOTAL      STRESS --------------------" << endl;
    	cout << " " << setw(8) << "STRESS" << endl;
    	ofs_running << " " << setw(12) << "STRESS" << endl;
	}

    
	for (int i=0; i<3; i++)
	{

		//if(TEST_STRESS)
		cout << " " << setw(15) << scs(i,0)*unit_transform << setw(15)
			<< scs(i,1)*unit_transform << setw(15) << scs(i,2)*unit_transform << endl;

		ofs_running << " " << setw(15) << scs(i,0)*unit_transform << setw(15)
			<< scs(i,1)*unit_transform << setw(15) << scs(i,2)*unit_transform << endl;

	}
	ofs_running << setiosflags(ios::left);
	cout << resetiosflags(ios::showpos);

    return;
}
