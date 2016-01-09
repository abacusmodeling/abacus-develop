//==========================================================
// Author: mohan
// Last Update : 2009-07-06
//==========================================================
#ifndef MPI_DATA_H
#define MPI_DATA_H

#include "../src_pw/tools.h" 
#include "../src_pw/unitcell.h"
#include "parallel_global.h"

class Parallel_PW//: public PW_Basis
{
public:
	
	Parallel_PW();
	~Parallel_PW();

	void init(
		const double &gcut_in,
		const int &n1_in,
		const int &n2_in,
		const int &n3_in,
		const int &bz,
		const int &nproc_in,
		const int &rank_in);

	void columns_map();
	void restore_st();
	void columns_and_pw_distribution();


	void max_pw_column(int &pw,int &max_i,int &max_j);
	void fft_dlay_set();
	void fft_map(int *ig2fft,const int ngm, const int &ngmc_g_in);
	void print_data(ofstream &print)const;

	int *isind;
	int *ismap;
	int *st_start;
	int *nst_per;
	int *npw_per;
	int *npps;
	int *st_i;
	int *st_j;
	int *st_n;

	int nst;//number of total columns

	bool *index;
	int *index_ip;
	int n1;
	int n2;
	int n3;
	int **st;
	int nproc_use;
	int rank_use;

private:
	double gcut;
	int *ig_l2g;
	int ngm_i;
	int ngm_i_record;
	bool allocate_igl2g;

};

#endif
