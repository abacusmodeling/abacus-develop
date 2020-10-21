//=========================================================
//AUTHOR : mohan
//DATE : 2008-11-29
//=========================================================
#ifndef ORBITAL_INFORMATION_0_H
#define ORBITAL_INFORMATION_0_H

#include "../src_pw/tools.h"

class Orbital_Information_0
{
	public:
	
	Orbital_Information_0();
	~Orbital_Information_0();

	void Read(const string &fn,ofstream &ofs);

	const int& getNtype(void) const {return ntype;}
	const string& getLabel(const int it) const {return label[it];}
	const string& getEnv(const int it)const {return enviroment[it]; }
	const double& getTran(const int it)const {return transparence[it];}
	const int* getNchi(const int it) const { return nchi[it]; }
	const int& getNchi_l(const int it,const int il){ return nchi[it][il]; }
	const int& getLmax(const int it) const { return lmax[it]; }

	private:

	int ntype;
	string *label;
	string *enviroment;
	double *transparence;
	int* lmax;
	int** nchi;

	bool init_done;


};

#endif
