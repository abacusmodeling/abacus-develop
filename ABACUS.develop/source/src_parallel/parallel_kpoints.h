#ifndef PARALLEL_KPOINTS_H
#define PARALLEL_KPOINTS_H

#include "../src_pw/tools.h"

class Parallel_Kpoints
{
	public:
	
	Parallel_Kpoints();
	~Parallel_Kpoints();	

	void init();
	void kinfo(int &nkstot);
	
	// collect value from each pool to wk.
	void pool_collection(double &value, const double *wk, const int &ik);

	// collect value from each pool to overlap.
	void pool_collection(double *valuea, double *valueb, const realArray &a, const realArray &b, const int &ik);

	// information about pool
	int *nproc_pool;
	int *startpro_pool;

	// inforamation about kpoints
	int* nks_pool;
	int* startk_pool;
	int* whichpool;

	private:

#ifdef __MPI
	void divide_pools(void);
	void get_nks_pool(const int &nkstot);
	void get_startk_pool(const int &nkstot);
	void get_whichpool(const int &nkstot);
#endif


};

#endif
