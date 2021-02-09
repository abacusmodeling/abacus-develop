#ifndef FORCE_LCAO_K_H
#define FORCE_LCAO_K_H

#include "../src_pw/tools.h"
#include "lcao_matrix.h" 
#include "force_lcao_gamma.h"

class Force_LCAO_k : public Force_LCAO_gamma
{
	public :
	
	Force_LCAO_k ();
	~Force_LCAO_k ();

	protected:
	
	//orthonormal force + contribution from T and VNL
	void ftable_k (void);

	private:

	// get the ds, dt, dvnl.
	void allocate_k(void);


	void finish_k(void);


	void set_EDM_k(double **dmR, const bool with_energy);


	// mohan add 2012-01-09
	complex<double> set_EDM_k_element(
		const complex<double> &phase,
		const bool with_energy,
		complex<double> &coef1, complex<double> &coef2,
		const double &ekb);

	
	// calculate the force due to < dphi | beta > < beta | phi >
	void cal_ftvnl_dphi_k(double **dmR);


	// calculate the force due to < phi | Vlocal | dphi >
	void cal_fvl_dphi_k(double **dmR);


	// calculate the overlap force
	void cal_foverlap_k(void);


	// calculate the force due to < phi | dbeta > < beta | phi >
	void cal_fvnl_dbeta_k(double** dmR);


	void test(double* mm, const string &name);
	
};
#endif
