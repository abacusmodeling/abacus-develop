//=========================================================
//AUTHOR : mohan
//DATE : 2009-04-23
//Last Update : 2009-04-23
//=========================================================
#ifndef LCAO_ORBITALS_H
#define LCAO_ORBITALS_H

#include "../src_pw/tools.h"
#include "ORB_atomic.h"
#include "ORB_atomic_lm.h"
#include "ORB_nonlocal.h"

class LCAO_Orbitals
{
	public:

	LCAO_Orbitals();
	~LCAO_Orbitals();

	void Read_Orbitals(void);

	void Read_PAO(const int& it);

	// in order to get rid of the .NONLOCAL file.
	void Set_NonLocal(const int &it, int &n_projectors);

	// read in the NONLOCAL projector from file.
	void Read_NonLocal(const int& it, int &n_projectors);

	void set_nl_index(void);

#ifdef __MPI
	void bcast_files(void);
#endif

	const double& get_ecutwfc(void) const {return ecutwfc;}
	const int& get_kmesh(void) const{return kmesh;}
	const double& get_dk(void) const {return dk;}
	const double& get_dR(void) const {return dR;}
	const double& get_Rmax(void) const {return Rmax;}
	const int& get_lmax(void) const {return lmax;}
	const int& get_nchimax(void) const {return nchimax;}
	const int& get_ntype(void) const {return ntype;}
	const double& get_dr_uniform(void) const {return dr_uniform;}

	Numerical_Orbital* Phi;
	Numerical_Nonlocal* Beta;
	
	// init in input.cpp
	double ecutwfc;
	double dk;
	double dR;
	double Rmax;
	int *nproj; //mohan add 2010-12-19
	int nprojmax; // mohan add 2010-03-07
	int nkb; // total number of projectors.
	IntArray itiaib2ib_all;
	IntArray ib2_ylm;
	
	double dr_uniform;

	// init in unitcell_pseudo
	// assume ntype < 20.
	bool read_in_flag;
	std::vector<string> orbital_file;
	std::vector<string> nonlocal_file;

	private:

	int kmesh;
	int lmax;
	int nchimax;
	int ntype;

};

extern LCAO_Orbitals ORB;
#endif
