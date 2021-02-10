#ifndef LCAO_NNR_H
#define LCAO_NNR_H

#include "grid_technique.h"

class LCAO_nnr
{
	public:

	LCAO_nnr();
	~LCAO_nnr();	
	
	//---------------------------------------
	// nnrg: number of matrix elements on
	// each processor's real space grid.
	// use: GridT.in_this_processor
	//---------------------------------------
	int nnrg;
	int *nlocdimg;
	int *nlocstartg;
	int *nad; // number of adjacent atoms for each atom.
	int **find_R2;
	int **find_R2st;
	bool allocate_find_R2;

	//---------------------------------------
	// number of elements in H or S matrix,
	// nnr -> 2D block distribution;
	// use: ParaO.trace_loc_row
	// use: ParaO.trace_loc_col
	//---------------------------------------
	int nnr; 
	int *nlocdim;
	int *nlocstart;

	void cal_nnr();	
	void cal_nnrg(const Grid_Technique &GT);


	// folding the fixed Hamiltonian (T+Vnl) if
	// k-point algorithm is used.
	void folding_fixedH(const int &ik);

	int cal_RindexAtom(const int &u1, const int &u2, const int &u3, const int &iat2);

	private:

	void cal_max_box_index(void);

	int maxB1, maxB2, maxB3;
	int minB1, minB2, minB3;
	int nB1, nB2, nB3;
	int nbox;

	
};

extern LCAO_nnr LNNR;

#endif 
