//=========================================================
//AUTHOR : mohan
//DATE : 2010-11-09
//=========================================================
#ifndef GINT_K_INIT_H
#define GINT_K_INIT_H
#include "../src_pw/tools.h"
class Gint_k_init
{
	friend class Gint_k;
	
public:

	void init(
		const int &nbx_in,
		const int &nby_in,
		const int &nbz_in,
		const int &nbz_start_in,
		const int &ncxyz_in);

protected:

	Gint_k_init();
	~Gint_k_init();

	int nbx;
	int nby;
	int nbz;
	int ncxyz;
	int nbz_start;
	int* nnn;

};

#endif
