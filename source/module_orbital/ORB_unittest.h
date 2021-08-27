#include "ORB_control.h"
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace std;

class test_orb
{
public:

	test_orb();
	~test_orb();

	LCAO_Orbitals ORB;
	ORB_gen_tables OGT;
	ORB_control ooo;

	std::ofstream ofs_running;

	void count_ntype(); //from STRU, count types of elements
	void set_files();   //from STRU, read names of LCAO files
	void set_ekcut();	//from LCAO files, read and set ekcut
	void read_files();	//interface to Read_PAO

	bool force_flag = 0;
	int my_rank = 0;
	int ntype;

	double lcao_ecut = 0; // (Ry)
	double lcao_dk = 0.01;
	double lcao_dr = 0.01;
	double lcao_rmax = 30; // (a.u.)

	int out_descriptor = 0;
	int out_r_matrix = 0;

	int lmax=1;
	double lat0=1.0;
};
