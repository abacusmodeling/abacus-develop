#include "ORB_control.h"

class test_orb
{
public:

	test_orb();
	~test_orb();

	LCAO_Orbitals ORB;

	void count_ntype();
	void set_files();

	double lcao_ecut = 0; // (Ry)
	double lcao_dk = 0.01;
	double lcao_dr = 0.01;
	double lcao_rmax = 30; // (a.u.)

	int out_descriptor = 0;
	int out_r_matrix = 0;
	bool force_flag = 0;

	int lmax;
	int my_rank = 0;

	int ntype;

};
