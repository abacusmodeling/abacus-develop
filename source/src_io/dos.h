#ifndef DOS_H
#define DOS_H
#include "../src_pw/tools.h"

namespace Dos
{
	bool calculate_dos(
		const int &is,
		const std::vector<int> &isk,
		const std::string &fn,// file address.
		const double &de_ev, // delta energy in ev.
		const double &emax_ev,// maximal energy in ev.
		const double &emin_ev,// minimal energy in ev.
		const int &nks,//number of k points
		const int &nkstot,
		const std::vector<double> &wk,//weight of k points
		const ModuleBase::matrix &wg,//weight of (kpoint,bands)
		const int &nbands,// number of bands
		double **ekb);//store energy for each k point and each band

	void calculate_Mulliken(const std::string &fn);

	void nscf_band(
		const int &is,
		const std::string &out_band_dir, 
		const int &nks, 
		const int &nband, 
		const double &fermie,
		double **ekb);

	void nscf_fermi_surface(const std::string &out_band_dir,
		const int &nks,
		const int &nband,
		double **ekb);
	
}

#endif 
