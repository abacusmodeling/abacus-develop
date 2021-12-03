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
		ModuleBase::matrix& foverlap,
		ModuleBase::matrix& ftvnl_dphi,
		ModuleBase::matrix& fvnl_dbeta,	
		ModuleBase::matrix& fvl_dphi,
		ModuleBase::matrix& soverlap,
		ModuleBase::matrix& stvnl_dphi,
		ModuleBase::matrix& svnl_dbeta,
		ModuleBase::matrix& svl_dphi
		);

	// get the ds, dt, dvnl.
	void allocate_k(void);

	void finish_k(void);

	void set_EDM_k(double** dm2d, const bool with_energy);

	// mohan add 2012-01-09
	std::complex<double> set_EDM_k_element(
		const std::complex<double> &phase,
		const bool with_energy,
		std::complex<double> &coef1, std::complex<double> &coef2,
		const double &ekb);
	
	// calculate the force due to < dphi | beta > < beta | phi >
	void cal_ftvnl_dphi_k(double** dm2d, const bool isforce, const bool isstress, ModuleBase::matrix& ftvnl_dphi, ModuleBase::matrix& stvnl_dphi);

	// calculate the overlap force
	void cal_foverlap_k(const bool isforce, const bool isstress, ModuleBase::matrix& foverlap, ModuleBase::matrix& soverlap);

	// calculate the force due to < phi | Vlocal | dphi >
	void cal_fvl_dphi_k(double** dm2d, const bool isforce, const bool isstress, ModuleBase::matrix& fvl_dphi, ModuleBase::matrix& svl_dphi);

	// old method to calculate the force due to < phi | dbeta > < beta | phi >
	void cal_fvnl_dbeta_k(double** dm2d, const bool isforce, const bool isstress, ModuleBase::matrix& fvnl_dbeta, ModuleBase::matrix& svnl_dbeta);
	// new method to calculate the force due to < phi | dbeta > < beta | phi > , developed by wenfei-li
	void cal_fvnl_dbeta_k_new(double** dm2d, const bool isforce, const bool isstress, ModuleBase::matrix& fvnl_dbeta, ModuleBase::matrix& svnl_dbeta);

	void test(double* mm, const std::string &name);

	// calculate the force due to < phi | dbeta > < beta | phi >
	void calFvnlDbeta(
		double** dm2d, 
		const bool &isforce, 
		const bool &isstress, 
		ModuleBase::matrix& fvnl_dbeta, 
		ModuleBase::matrix& svnl_dbeta,
		const int &vnl_method);

};
#endif
