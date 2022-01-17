#ifndef WF_IO_H
#define WF_IO_H

#include "../src_pw/tools.h"
#include "../module_base/complexmatrix.h"
#include "rwstream.h"

namespace WF_io
{
    void write_wfc(const std::string &fn, const ModuleBase::ComplexMatrix *psi);
	
	// mohan add 2011-02-21
	void write_wfc2(const std::string &fn, const ModuleBase::ComplexMatrix *psi,const ModuleBase::Vector3<double> *gkk);
    void read_wfc(const std::string &fn, const ModuleBase::ComplexMatrix *psi);
}

#endif
