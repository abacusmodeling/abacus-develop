#ifndef STRESS_H
#define SSTRESS_H

#include "tools.h"

class Stress 
{
public:
    Stress(){};
    ~Stress(){};

	void init();
	
private:
	
	Matrix3 skin;
	Matrix3 sloc;
	Matrix3 shar;
	Matrix3 sxc;
	Matrix3 snlcc;
	Matrix3 sion;

	void s_loc();
	void s_hartree();
	void s_nlcc();
	void s_ewald();	
	void s_t_vnl();
	void s_exx();

};

#endif
