#ifndef NUMERICAL_VNA_LM_H
#define NUMERICAL_VNA_LM_H

#include "../src_pw/tools.h"
#include "../src_parallel/fftw.h"

//=========================================================
//CLASS Numerical_Vna_Lm
//Note : contain information about each Vna projector
//=========================================================

class Numerical_Vna_Lm
{
	friend class Neutral_Pot;

public:

	Numerical_Vna_Lm();
	~Numerical_Vna_Lm();

	const int& getL() const { return angular_momentum_l; }
	const int& getType() const { return index_atom_type; }
	const int& getChi() const {return index_chi; }

	const double* getRadial() const { return r_radial; }
	const double& getRadial(const int ir) const { return r_radial[ir]; }

	const double* getPsi() const { return psi;}
	const double* getPsi_r() const { return psir; }
	const double& getPsi_r(const int ir) const { return psir[ir]; }

	const double& getDk()const { return dk; }
	const double* getKpoint()const { return k_radial; }
	const double& getKpoint(const int ik) const { return k_radial[ik]; }
	const double* getPsi_k() const { return psik; }
	const double& getPsi_k(const int ik) const { return psik[ik]; }
	
//==========================================================
// EXPLAIN : set information about Numerical_Orbital_Lm
// MEMBER FUNCTION :
//==========================================================
	void set_vna_proj_info
	(
 		const string &label_in,
	 	const int &index_atom_type_in,
		const int &angular_momentum_l_in,
		const int &index_chi_in,

	    const int &nr_in,
		const double *rab_in,
		const double *r_radial_in,
		const double *psi_in,

		const int &nk_in,
		const double &dk_in,

		const double &lat0
	);

private:

	void cal_kradial(void);
	void norm_test()const;
	void plot()const;

	string label;
	int index_atom_type;
	int angular_momentum_l;
	int index_chi;

	int nr;
	int nk;

	double rcut;
	double kcut;
	double dk;

	double* r_radial; //points of r
	double* k_radial;

	double* rab;

	double* psi;//psi(r)
	double* psir; //psi(r) * r
	double* psik;

};

#endif

