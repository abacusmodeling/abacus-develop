#ifndef MD_FUNC_H
#define MD_FUNC_H

#include "MD_parameters.h"
#include "../module_cell/unitcell_pseudo.h"

class MD_func
{
    public:

    MD_func(){};
    ~MD_func(){};

	static double gaussrand();

	static bool RestartMD(const int& numIon, ModuleBase::Vector3<double>* vel, int& step_rst);
    static void mdRestartOut(const int& step, const int& recordFreq, const int& numIon, ModuleBase::Vector3<double>* vel);
	
	static void InitPos( const UnitCell_pseudo &unit_in, ModuleBase::Vector3<double>* pos);

	static void InitVel(
		const UnitCell_pseudo &unit_in, 
		const double& temperature, 
		double* allmass,
		int& frozen_freedom,
		ModuleBase::Vector3<int>* ionmbl,
		ModuleBase::Vector3<double>* vel);

	static void ReadVel(const UnitCell_pseudo &unit_in, ModuleBase::Vector3<double>* vel);

	static void RandomVel(
		const int& numIon, 
		const double& temperature,
		const double* allmass,
		const int& frozen_freedom,
		const ModuleBase::Vector3<int> frozen,
		const ModuleBase::Vector3<int>* ionmbl,
		ModuleBase::Vector3<double>* vel);

	static void force_virial(
		const MD_parameters &mdp,
		const UnitCell_pseudo &unit_in,
		double &potential,
		ModuleBase::Vector3<double> *force,
		ModuleBase::matrix &stress);

	static double GetAtomKE(
		const int &numIon,
		const ModuleBase::Vector3<double> *vel, 
		const double *allmass);

	static void kinetic_stress(
		const UnitCell_pseudo &unit_in,
		const ModuleBase::Vector3<double> *vel, 
		const double *allmass, 
		double &kinetic,
		ModuleBase::matrix &stress);

	static void outStress(const ModuleBase::matrix &virial, const ModuleBase::matrix &stress);

	static void MDdump(
		const int &step, 
		const int &natom,
		const ModuleBase::matrix &virial, 
		const ModuleBase::Vector3<double> *force);

	static std::string intTurnTostring(long int iter,std::string path);
	static void getMassMbl(const UnitCell_pseudo &unit_in, 
		double* allmass, 
		ModuleBase::Vector3<int> &frozen,
		ModuleBase::Vector3<int>* ionmbl);
    static void printpos(const std::string& file, const int& iter, const int& recordFreq, const UnitCell_pseudo& unit_in);
    static void scalevel(
		const int& numIon,
		const int& nfrozen,
		const double& temperature,
		ModuleBase::Vector3<double>* vel,
		const double* allmass);
	static double MAXVALF(const int numIon, const ModuleBase::Vector3<double>* force);
	static double Conserved(const double KE, const double PE, const int number);
};
#endif
