#ifndef FORCE_LCAO_K_H
#define FORCE_LCAO_K_H

#include "../src_pw/tools.h"
#include "LCAO_matrix.h" 
#include "FORCE_gamma.h"

class Force_LCAO_k : public Force_LCAO_gamma
{
	public :
	
	friend class Force_Stress_LCAO;

	Force_LCAO_k ();
	~Force_LCAO_k ();

	private:
	
	//orthonormal force + contribution from T and VNL
	void ftable_k (
		const bool isforce,
		const bool isstress,
		matrix& foverlap,
		matrix& ftvnl_dphi,
		matrix& fvnl_dbeta,	
		matrix& fvl_dphi,
		matrix& soverlap,
		matrix& stvnl_dphi,
		matrix& svnl_dbeta,
		matrix& svl_dphi
		);

	// get the ds, dt, dvnl.
	void allocate_k(void);

	void finish_k(void);

	void set_EDM_k(double** dm2d, const bool with_energy);

	// mohan add 2012-01-09
	complex<double> set_EDM_k_element(
		const complex<double> &phase,
		const bool with_energy,
		complex<double> &coef1, complex<double> &coef2,
		const double &ekb);
	
	// calculate the force due to < dphi | beta > < beta | phi >
	void cal_ftvnl_dphi_k(double** dm2d, const bool isforce, const bool isstress, matrix& ftvnl_dphi, matrix& stvnl_dphi);

	// calculate the overlap force
	void cal_foverlap_k(const bool isforce, const bool isstress, matrix& foverlap, matrix& soverlap);

	// calculate the force due to < phi | Vlocal | dphi >
	void cal_fvl_dphi_k(double** dm2d, const bool isforce, const bool isstress, matrix& fvl_dphi, matrix& svl_dphi);

	// calculate the force due to < phi | dbeta > < beta | phi >
	void cal_fvnl_dbeta_k(double** dm2d, const bool isforce, const bool isstress, matrix& fvnl_dbeta, matrix& svnl_dbeta);


	void test(double* mm, const string &name);

};
#endif
