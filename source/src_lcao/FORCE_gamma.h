#ifndef FORCE_LCAO_GAMMA_H
#define FORCE_LCAO_GAMMA_H

#include "../src_pw/tools.h"
#include "LCAO_matrix.h" 

class Force_LCAO_gamma
{
	public :

	friend class Force_Stress_LCAO;

	Force_LCAO_gamma ();
	~Force_LCAO_gamma ();

	private:
	
	//orthonormal force + contribution from T and VNL
	void ftable_gamma (
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
	void allocate_gamma(void);

	void finish_ftable_gamma(void);

	void average_force(double* fm);

	void test_gamma(double* mm, const std::string &name);


	//-------------------------------------------------------------
	// forces reated to overlap matrix
	// forces related to energy density matrix 
	//-------------------------------------------------------------
	void set_EDM_gamma(ModuleBase::matrix& dm, bool with_energy);

	double set_EDM_element(const int &ii, const int &jj, const bool with_energy, 
	 double*** coef1, double*** coef2, const int &is);

	void cal_foverlap(
		const bool isforce, 
		const bool isstress, 
		ModuleBase::matrix& foverlap, 
		ModuleBase::matrix& soverlap);	

	//-------------------------------------------------------------
	// forces related to kinetic and non-local pseudopotentials
	//--------------------------------------------------------------
	void cal_ftvnl_dphi(
		ModuleBase::matrix& dm2d, 
		const bool isforce, 
		const bool isstress, 
		ModuleBase::matrix& ftvnl_dphi, 
		ModuleBase::matrix& stvnl_dphi);

	void cal_ftvnl_dphi(
		const std::vector<ModuleBase::matrix> &dm2d, 
		const bool isforce, 
		const bool isstress, 
		ModuleBase::matrix& ftvnl_dphi, 
		ModuleBase::matrix& stvnl_dphi);

	void cal_fvnl_dbeta(
		ModuleBase::matrix& dm2d, 
		const bool isforce, 
		const bool isstress, 
		ModuleBase::matrix& fvnl_dbeta, 
		ModuleBase::matrix& svnl_dbeta);

	void cal_fvnl_dbeta(
		const std::vector<ModuleBase::matrix> &dm2d, 
		const bool isforce, 
		const bool isstress, 
		ModuleBase::matrix& fvnl_dbeta, 
		ModuleBase::matrix& svnl_dbeta);

	void cal_fvnl_dbeta_new(
		ModuleBase::matrix& dm2d, 
		const bool isforce, 
		const bool isstress, 
		ModuleBase::matrix& fvnl_dbeta, 
		ModuleBase::matrix& svnl_dbeta);

	void cal_fvnl_dbeta_new(
		const std::vector<ModuleBase::matrix> &dm2d, 
		const bool isforce, 
		const bool isstress, 
		ModuleBase::matrix& fvnl_dbeta, 
		ModuleBase::matrix& svnl_dbeta);

	//-------------------------------------------
	// forces related to local pseudopotentials
	//-------------------------------------------
	void cal_fvl_dphi(
		ModuleBase::matrix& dm2d, 
		const bool isforce, 
		const bool isstress, 
		ModuleBase::matrix& fvl_dphi, 
		ModuleBase::matrix& svl_dphi);

	void cal_fvl_dphi(
		const std::vector<ModuleBase::matrix> &dm2d, 
		const bool isforce, 
		const bool isstress, 
		ModuleBase::matrix& fvl_dphi, 
		ModuleBase::matrix& svl_dphi);

	void calFvnlDbeta(
		const std::vector<ModuleBase::matrix> &dm2d, 
		const bool &isforce, 
		const bool &isstress, 
		ModuleBase::matrix& fvnl_dbeta, 
		ModuleBase::matrix& svnl_dbeta,
		const int &vnl_method);

	public:
	// calculate dVnl=<phi|dVnl|dphi> in LCAO
	void NonlocalDphi(const int& nspin, const int& vnl_method, const bool& cal_deri);
};
#endif
