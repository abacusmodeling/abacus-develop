#ifndef WF_LOCAL_H
#define WF_LOCAL_H

#include "../src_pw/tools.h"

namespace WF_Local
{
	// mohan add 2010-09-09
	void write_lowf(const string &name, double** ctot);
	void write_lowf_complex(const string &name, complex<double>** ctot, const int &ik);

	// mohan add 2010-09-10
	void distri_lowf(double** ctot, double **c);
	void distri_lowf_complex(complex<double>** ctot, complex<double> **cc);

	// mohan add 2010-09-24
	void distri_lowf_aug(double** ctot, double **c_aug);
	// mohan add 2012-01-09
	void distri_lowf_aug_complex(complex<double>** ctot, complex<double> **c_aug);

	// mohan add 2010-09-10
	int read_lowf(double **c);

	// mohan add 2012-03-04
	int read_lowf_complex(complex<double> **c, const int &ik);
}

#endif
