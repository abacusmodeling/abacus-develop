//=========================================================
//AUTHOR : mohan
//DATE : 2012-03-29
//=========================================================
#ifndef GINT_SPEED_H
#define GINT_SPEED_H

#include "../src_pw/tools.h"
#include "grid_base_beta.h"
#include "grid_technique.h"
#include "lcao_matrix.h"

class Gint_Speed : public Grid_Base_Beta
{
private:

	int grid_index;
	int max_size;
	
	void gamma_vlocal(void);  
	double gamma_charge(void);
	void save_phi(void);

	int* nr; //[nat], number of phi value in real space
	double*** phiylm; //[nat,nw,nr[iat]], save phi in real space
	bool allocate_phiylm;

	int** grid_label;//[nat,nr], record the grid information for each phi(r)
	int* vindex;//index for local potential
	int nov;//max number of overlap points in real space.

public:

	Gint_Speed();
	~Gint_Speed();

	double cal_rho(void);
	void cal_vlocal( const double* vlocal_in);

};

#endif
