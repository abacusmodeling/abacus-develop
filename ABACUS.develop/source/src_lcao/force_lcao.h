#ifndef FORCE_LCAO_H
#define FORCE_LCAO_H

#include "../src_pw/tools.h"
#include "force_lcao_k.h"
#include "stress_lcao.h"
#include "../input_conv.h"

class Force_LCAO : private Force_LCAO_k
{
	// mohan add 2021-02-09
	friend class md;
	friend void Input_Conv::Convert();
	friend class Update_input;
	friend class Local_Orbital_Ions;

	public :
	
	Force_LCAO ();
	~Force_LCAO ();

	private:

	void allocate (void);
	void destroy (void);

	void cal_force_loc (void);
	void cal_force_ew (void);
	void cal_force_scc (void);
	void cal_force_cc (void);
	
	void print_force(const string &name, double** f, const bool screen, bool ry)const;
	void printforce_total (bool ry);
	
	void start_force(void);

	void cal_stress(matrix &stress);
	
	// total force
	matrix fcs; 
	static double force_invalid_threshold_ev; // mohan add 2011-04-17

	matrix stress_vdw;//zhengdy added for vdw-stress-part

	//each part of force
	double** fvl_dvl;
	double** fewalds; 
	double** fcc;
	double** fscc;
	
	bool allocate_flag;
	double output_acc; // control the accuracy
};
#endif
