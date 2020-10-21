#ifndef SYMM_OTHER_H
#define SYMM_OTHER_H

#include "tools.h"

namespace Symm_Other
{
	void print1(const int &ibrav, const double *cel_const);
	bool right_hand_sense(Vector3<double> &v1,Vector3<double> &v2,Vector3<double> &v3);
	double celvol(const Vector3<double> &a, 
		const Vector3<double> &b, const Vector3<double> &c);


}

#endif
