#ifndef MD_H
#define MD_H

#include<iostream>
#include<string>
#include<fstream>
#include<stdlib.h>
#include<time.h>
#include<cmath>
#include"tools.h"
#include"global.h"
#include"../src_lcao/FORCE_STRESS.h"
#include"unitcell_pseudo.h"

using namespace std;

//const double fundamentalTime = 2.418884326505e-17; 

//-------------------------------------------------------------------
// mohan reconstruction note: 2021-02-07
// 1) explain every variable and every function
// 2) if too many, divide the class into more files 
// 3) leave a space between different functions
// 4) write the comments above the variables or functions
//-------------------------------------------------------------------

class md
{
	public:

	md(int n=1);

	void md_allocate();
	void initMD();
	bool RestartMD();
	void mstout(int step);
	void RemoveMovementOfCenterOfMass();
	double GetAtomKE();
	void MakeIntCoeff();
	void InitVelocity();
	void MonitorMeanSquareDisplacement(int step);
	void OutputIonVelocityAndStats(int iter);
	void OutputMDHeader();
	void OutputMDGeom(int iter);
	void ReadNewTemp(int step);
	string intTurnTostring(int iter,string path);
	void connection0();
	void connection1();
	void connection2();
	void callforce();
	void moveatoms(int step);
    void printpos(string file,int iter);
    void md_release();
    void scalevel();
    void RDF();
    void printRDF(int step);
    void PDF(int step);

	// variables:

	int mdtype; // 1: NVT:Nose Hoover, 2:NVT:Nose Hoover add velocity scaling 3: NPT, 0: NVE 
	double xLogS; // thermostat's position
	double vLogS;// thermostat's speed
	Vector3<double> *vel;//   velocity of each atom, unit is a.u.
	Matrix3 vBoxG;        // barostat's velocity
	Matrix3 ionlatvec;
	double  tauthermo;  // The thermostat fluctuation period, for Al fcc at 300K ~ 1000fs
	double taubaro; // The barostat   fluctuation period, for Al fcc at 300K ~ 500 fs
	double Qmass;  // Inertia of extended system variable
	double dt ; // Time increment (hbar/E_hartree)
	double temperature;  //Temperature (in Hartree, 1 Hartree ~ 3E5 K)
	bool doMSD;   // compute the mean square displacement.
	bool doMSDAtom;        // compute the mean square displacement for each atom.
	double msdstartTime;   //beyond this time, msd is going to be computed
	double *msd    ;
	double diffuCoeff ;
	int numIon;
	int ntype;
	int* na;
	int nfrozen;
	double *allmass;
	Vector3<double> *force;
	Vector3<int> *ionmbl;
	Vector3<double> *taudirac;           //dirac coord of atoms
	Vector3<double> *cartNoWrap;        //cartensian coord of atoms, *not* wrapped
	Vector3<double> *msdcoords0;        
	Vector3<double> *msdcoordst;
	int nResn;   // Coarseness of Trotter factorization of Nose-Hoover thermostat
	int nYosh;  // Order of thermostat scheme corresponding with "intCoeff"
	double *intCoeff;
	int outputstressperiod;// The period to output stress
	int dumpmdfreq;       // The period to dump MD information for monitoring and restarting MD
	int fixTemperature;       // The period to read in new temperature during MD run.
	string mdoutputpath;         // output directory of md files: .ion .vel
	int rstMD;// 1 : restart MD, vel.restart and ion.restart files will be read
	// 0 : not restart from ion and vel files

	double ediff;           //parameter for constrain 
	double ediffg;          //parameter for constrain
	int startStep;// the first step number of MD
	bool velRescale;// whether to rescale velocities at the first step
	ofstream out;
	double mdenergy;
	const static double fundamentalTime;
	double rdf[100];
	int step_rst; 

	matrix stress_lcao;
};

#endif
