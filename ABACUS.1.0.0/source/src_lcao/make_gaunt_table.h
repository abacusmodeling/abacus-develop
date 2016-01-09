//=========================================================
//AUTHOR : Mohan 
//DATE : 2009-04-23
//=========================================================
#ifndef MAKE_GAUNT_TABLE_H
#define MAKE_GAUNT_TABLE_H

#include "../src_pw/tools.h"

class Make_Gaunt_Table
{
	public:

	Make_Gaunt_Table();
	~Make_Gaunt_Table();

	/***************************************
	 * Method 2
	 * using WIgner 3j symbols 
	 * Y(l1,m1), Y(l2,m2), Y(L,M)
	 * *************************************/
	
	void init_Gaunt_CH(const int& Lmax);
	double Get_Gaunt_CH(
							const int& l1,
							const int& m1,
							const int& l2,
							const int& m2,
							const int& l3,
							const int& m3	);

	//M defined here are restricted within 0 to 2l+1
	//should be transformed first
	double Get_Gaunt_SH(
							const int& l1,
							const int& mm1,
							const int& l2,
							const int& mm2,
							const int& l3,
							const int& mm3	);

	double Calc_Gaunt_CH(
							const int& l1,
							const int& m1,
							const int& l2,
							const int& m2,
							const int& l3,
							const int& m3	);
					

	/***************************************
	 * MEthod 2
	 * Directly Calculate integral of 
	 * S(l1,m1), S(l2,m2), S(L,M)
	 * ************************************/
	realArray Gaunt_Coefficients;

	//============================================================
	// (1) Make Ylm_Gaunt Table.
	//============================================================
	void init_Ylm_Gaunt(const int &lmax, const double &s1,const double &e1,
		const double &s2,const double &e2);

	//============================================================
	// (2) Use Ylm_Gaunt to calculate Gaunt Coefficinets element
	//============================================================
	double Cal_Gaunt_single(
	   	const int &l1, const int &m1, 
	   	const int &l2, const int &m2, 
	   	const int &l, const int &m, 
	   	const double &s1, const double &e1,    
	   	const double &s2, const double &e2);

	//============================================================
	// (3) Make the whole Gaunt Coefficients table
	//============================================================
	void init_Gaunt(const int &lmax);


	//========================================================
	// Small function
	//========================================================
	static int get_lm_index(const int l, const int m);

	private:
	
	//Index Function
	//Yu's mehtod
	void ZEROS();
	
	int P_EL(const int& L);
	int EP_EL(const int& L);
	int index_func(
					const int& l1,
					const int& l2,
					const int& l3,
					const int& m3	);
	
	double Fact(const int& n);
	void Swap(
				int& l1,
				int& m1,
				int& l2,
				int& m2	);
	
	int Index_M(const int& m);
	//2*Lmax+1
	//Lmax <= 6
	//dim1 <= 5000
	//dim2 <= 30
	double Gaunt_CH[5000][30];
	
	//direct integral
	matrix Ylm_Gaunt;
};
#endif
