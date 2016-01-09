#ifndef DOS_H
#define DOS_H
#include "../src_pw/tools.h"

namespace Dos
{
	bool calculate_dos(
		const int &is,
		const int *isk,
		const string &fn,// file address.
		const double &de_ev, // delta energy in ev.
		const double &emax_ev,// maximal energy in ev.
		const double &emin_ev,// minimal energy in ev.
		const int &nks,//number of k points
		const int &nkstot,
		const double *wk,//weight of k points
		const matrix &wg,//weight of (kpoint,bands)
		const int &nbands,// number of bands
		double **ekb);//store energy for each k point and each band

	void calculate_Mulliken(const string &fn);

	void nscf_band(
		const int &is,
		const string &out_band_dir, 
		const int &nks, 
		const int &nband, 
		const double &fermie,
		double **ekb);

	void nscf_fermi_surface(const string &out_band_dir,
		const int &nks,
		const int &nband,
		double **ekb);
	
}

#endif 
