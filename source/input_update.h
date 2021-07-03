//==========================================================
// Author: mohan
// DATE : 2011-06-06
//==========================================================
#ifndef UPDATE_INPUT_H
#define UPDATE_INPUT_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <stdio.h>
#include <string.h>

using namespace std;

class Update_input
{
public:
    
	Update_input();
    ~Update_input();
    void init(const string &fn);

private:
	
    int force;
    double force_thr;		// threshold of force in unit (Ry/Bohr)
	double force_thr_ev2;	// invalid force threshold, mohan add 2011-04-17
    double dr2;				// \sum |rhog_out - rhog_in |^2
    int niter;				// number of max elec iter
    int nstep;				// number of max ionic iter
    double mixing_beta;		// 0 : no_mixing
	int printe;
	//int extra_pot;	// mohan add 2011-03-13, xiaohui modify 2015-02-01
	string charge_extrap;	//xiaohui add 2015-02-01
    int out_charge;		// output charge density.
	int out_dm; // output density matrix.
	int out_dos;			// dos calculation. mohan add 20090909
	int	out_lowf;			// output the wave functions in local basis.
    
	bool Read(const string &fn);
#ifdef __MPI
    void Bcast(void);
#endif

	template <class T>
	static void change(ofstream &ofs, const string &name, T &v1, T &v2)
	{
		ofs << " " << name << " change from " << v1 << " to " << v2 << endl;
	}

    template <class T>
    static void read_value(ifstream &ifs, T &var)
    {
        ifs >> var;
        ifs.ignore(150, '\n');
        return;
    }
};

#endif //UPDATE_INPUT
