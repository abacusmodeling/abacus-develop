#ifndef MD_BASIC_H
#define MD_BASIC_H

#include "MD_func.h"
#include "MD_thermo.h"
#include "MD_parameters.h"
#include "MD_fire.h"
#include "../module_cell/unitcell_pseudo.h"

//using namespace std;
class MD_basic
{
	public:	

	MD_basic(MD_parameters& MD_para_in, UnitCell_pseudo &unit_in);
	~MD_basic();

	void runNVT(int step1, double potential, ModuleBase::Vector3<double> *force, const ModuleBase::matrix &stress);//NVT ensemble MD
	void runNVE(int step1, double potential, ModuleBase::Vector3<double> *force, const ModuleBase::matrix &stress); //NVE ensemble MD
	bool runFIRE(int step1, double potential, ModuleBase::Vector3<double> *force, const ModuleBase::matrix &stress); //relax method FIRE
	int getRealStep();

	private:
    MD_func mdf;
    MD_thermo mdt;
	MD_parameters &mdp;
	UnitCell_pseudo &ucell;
	MD_fire fire;

	double temperature_;
	const static double fundamentalTime;
	int outputstressperiod_;
	int step_rst_;
	int step_;
	double energy_;
	int frozen_freedom_;
	double oldEtot_;

	//allocable variaints
	ModuleBase::Vector3<double> *vel;//   velocity of each atom, unit is a.u.
	double *allmass;     //mass of each atom
	//ModuleBase::Vector3<double> *force;  //force of each atom
	//matrix stress;           //stress for this lattice
	ModuleBase::Vector3<int> *ionmbl;    //if frozen of each atom
	ModuleBase::Vector3<double> *tauDirectChange;           //change of dirac coord of atoms
	ModuleBase::Vector3<double> *cart_change;        //cartensian coord of atoms, *not* wrapped

	//repete code
	void update_half_velocity(ModuleBase::Vector3<double> *force);
	void update_half_direct(const bool is_restart);
	void save_output_position();
	void outStressMD(const ModuleBase::matrix& stress, const double& twiceKE);
	void getTaudUpdate();
};

#endif
