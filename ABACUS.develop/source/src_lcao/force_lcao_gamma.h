#ifndef FORCE_LCAO_GAMMA_H
#define FORCE_LCAO_GAMMA_H

#include "../src_pw/tools.h"
#include "LCAO_matrix.h" 

class Force_LCAO_gamma
{
	public :
	
	Force_LCAO_gamma ();
	~Force_LCAO_gamma ();

	protected:
	
	//orthonormal force + contribution from T and VNL
	void ftable_gamma (void);

	//each part of force
	double** foverlap;
	double** ftvnl_dphi;
	double** fvnl_dbeta;	
	double** fvl_dphi;

	public:
	//each part of stress
	double soverlap[3][3];
	double stvnl_dphi[3][3];
	double svnl_dbeta[3][3];
	double svl_dphi[3][3];	

	private:

	// get the ds, dt, dvnl.
	void allocate_gamma(void);
	void finish_ftable_gamma(void);
	void set_EDM_gamma(double* dm, bool with_energy, const int &is);

	// mohan fix bug 2011-06-15
	double set_EDM_element(const int &ii, const int &jj, const bool with_energy, 
	 double*** coef1, double*** coef2, const int &is);

	void cal_foverlap(void);
	void cal_ftvnl_dphi(double** dm2d);
	void cal_fvnl_dbeta(double** dm2d);
	void cal_fvl_dphi(double** dm2d);

	void cal_ftvnl_dphi(const std::vector<matrix> &dm2d);
	void cal_fvnl_dbeta(const std::vector<matrix> &dm2d);
	void cal_fvl_dphi(const std::vector<matrix> &dm2d);

	void DerivS_PW (void);
	void DerivT_PW (void);

	void average_force(double* fm);

	void test_gamma(double* mm, const string &name);
	
};
#endif
