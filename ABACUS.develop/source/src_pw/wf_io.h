#ifndef WF_IO_H
#define WF_IO_H

#include "../src_pw/tools.h"
#include "../src_global/rwstream.h"

namespace WF_io
{
    void write_wfc(const string &fn, const ComplexMatrix *psi);
	
	// mohan add 2011-02-21
	void write_wfc2(const string &fn, const ComplexMatrix *psi,const Vector3<double> *gkk);
    void read_wfc(const string &fn, const ComplexMatrix *psi);
}

#endif
