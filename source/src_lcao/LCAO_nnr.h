#ifndef LCAO_NNR_H
#define LCAO_NNR_H

#include "grid_technique.h"
#include "src_lcao/LCAO_matrix.h"

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
	void folding_fixedH(const int &ik, LCAO_Matrix &lm);

	int cal_RindexAtom(const int &u1, const int &u2, const int &u3, const int &iat2);

	private:

	void cal_max_box_index(void);

	int maxB1;
	int maxB2;
	int maxB3;

	int minB1;
	int minB2;
	int minB3;

	int nB1;
	int nB2;
	int nB3;

	int nbox;
};

namespace GlobalC
{
extern LCAO_nnr LNNR;
}

#endif 
