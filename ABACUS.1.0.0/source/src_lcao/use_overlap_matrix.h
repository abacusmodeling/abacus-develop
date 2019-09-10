/************************************
//LiaoChen Modify on 2010-3-22
***********************************/
#ifndef USE_OVERLAP_MATRIX_H
#define USE_OVERLAP_MATRIX_H

#include "../src_pw/tools.h"
#include "use_overlap_table.h"
#include "sltk_grid_driver.h"

class Use_Overlap_Matrix
{
	friend class Force_LCAO_gamma;
	friend class Force_LCAO_k;
	friend class energy;//qifeng 2019/9/10
        friend class  Mulliken_Charge;//qifeng  2019/9/10
	public:
	Use_Overlap_Matrix();
	~Use_Overlap_Matrix();

	void calculate_NL_no(void);
	void calculate_T_no(void);
	void calculate_S_no(void);

	private:

	void build_ST_new(const char& dtype, const bool& cal_deri);	
	
	// can used in gamma algorithm.
	void build_Nonlocal_beta (const bool& calc_deri);

	// used if k point is used,
	// also can used in gamma algorithm.
	void build_Nonlocal_mu (const bool &calc_deri); 
	void test_Nonlocal(void);

};

#endif
