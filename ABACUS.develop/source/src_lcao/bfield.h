#ifndef BFIELD_H
#define BFIELD_H

#include <complex>
#include "../src_global/vector3.h"
using namespace std;

class Bfield
{
	public:
	Bfield();
	~Bfield();

	double tesla_x;
	double tesla_y;
	double tesla_z;
	double Rydberg_x;
	double Rydberg_y;
	double Rydberg_z;  //sun zhiyuan added at 2011-12-26//
	double Gauge_Origin_x;
	double Gauge_Origin_y;
	double Gauge_Origin_z; //sun zhiyuan added at 2011-12-26//
    double fac_of_phase;   //sun zhiyuan added at 2011-12-31//
    double bohr_mag;   //sun zhiyuan added at 2012-02-17//
    bool Ready_A_of_Atom;  //judge whether A_of_Atom has been done, 2012-01-05//

    double fac;  //transform from IU to Rydberg//
	double c;    //speed of light in Rydberg unit//
     
    int nlocal;  //nlocal is the number of local basis this processor includes//
	int *trace;  //trace[NLOCAL]  this array gives the index of a basis among the bases in this processor according to its index among the whole bases//

	double** A_of_Atom;   //nAtom*3, Zhiyuan add at 2011-12-24, former version:RxB// 
    complex<double>** Tab;

    void check();   //Zhiyuan add 2011-12-27 in order to check the input//	
	void convert(); //Convert bfield to Rydberg unit, sun zhiyuan added at 2011-12-26//
	void cal_A_of_Atom(void);
	void make_table(void);
	bool allocate_tab;

    complex<double> cal_phase(double A[3], double r[3], int sign);  //Zhiyuan add 2011-12-31//
	void cal_A(double A[3],double r[3]);  //sun zhiyuan add at 2011-08-07//
	static void a_x_b(const double* a, const double* b, double*c);
	static double a_dot_b(const double* a, const double* b); 
	static double a_dot_b(const Vector3<double> &a, const double *b);

	void calculate_NL_B(void);

	void snap_psibeta(complex<double> &nlm, const int& job,
			const Vector3<double> &R1,
			const int &T1,
			const int &iw1_all,
			const Vector3<double> &R2,
			const int &T2,
			const int &iw2_all,
			const Vector3<double> &R0,// The projector.
			const int &T0,
			const int &I0
			) const;
	void get_nlocal(void); //sun zhiyuan add at 2011-08-22, for the dimension of <beta|psi> table//

	void add_zeeman();

	

};

extern Bfield bfid;

#endif
