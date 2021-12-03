#ifndef MD_FUNC_H
#define MD_FUNC_H

#include "MD_parameters.h"
#include "../module_cell/unitcell_pseudo.h"

class MD_func
{
    public:

    MD_func(){};
    ~MD_func(){};

	bool RestartMD(const int& numIon, ModuleBase::Vector3<double>* vel, int& step_rst);
    void mdRestartOut(const int& step, const int& recordFreq, const int& numIon, ModuleBase::Vector3<double>* vel);
	double GetAtomKE(const int& numIon, const ModuleBase::Vector3<double>* vel, const double* allmass);
	void InitVelocity(
		const int& numIon, 
		const double& temperature, 
		const double& fundamentalTime, 
		const double* allmass,
		ModuleBase::Vector3<double>* vel);
	
	void InitVel(
		const UnitCell_pseudo &unit_in, 
		const double& temperature, 
		double* allmass,
		int& frozen_freedom,
		ModuleBase::Vector3<int>* ionmbl,
		ModuleBase::Vector3<double>* vel);

	void ReadVel(
		const UnitCell_pseudo &unit_in, 
		ModuleBase::Vector3<double>* vel);

	void RandomVel(
		const int& numIon, 
		const double& temperature,
		const double* allmass,
		const int& frozen_freedom,
		const ModuleBase::Vector3<int>* ionmbl,
		ModuleBase::Vector3<double>* vel);

//	void ReadNewTemp(int step);
	std::string intTurnTostring(long int iter,std::string path);
	int getMassMbl(const UnitCell_pseudo &unit_in, 
		const MD_parameters &mdp,
		double* allmass, 
		ModuleBase::Vector3<int>* ionmbl);
    void printpos(const std::string& file, const int& iter, const int& recordFreq, const UnitCell_pseudo& unit_in);
    void scalevel(
		const int& numIon,
		const int& nfrozen,
		const double& temperature,
		ModuleBase::Vector3<double>* vel,
		const double* allmass);
	double MAXVALF(const int numIon, const ModuleBase::Vector3<double>* force);
	double Conserved(const double KE, const double PE, const int number);
};
#endif
