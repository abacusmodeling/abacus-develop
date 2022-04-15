#ifndef UNKOVERLAP_PW
#define UNKOVERLAP_PW

#include <complex>
#include <cmath>
#include <fstream>
#include <string>
#include "../module_base/complexmatrix.h"
#include "../module_base/lapack_connector.h"
#include "../module_base/global_variable.h"
#include "../src_parallel/parallel_reduce.h"
#include "../module_base/vector3.h"

using namespace std;

class unkOverlap_pw
{
public:

	unkOverlap_pw();
	~unkOverlap_pw();	
	std::complex<double> unkdotp_G(const int ik_L, const int ik_R, const int iband_L, const int iband_R, const ModuleBase::ComplexMatrix *evc);
	std::complex<double> unkdotp_G0(const int ik_L, const int ik_R, const int iband_L, const int iband_R, const ModuleBase::ComplexMatrix *evc, const ModuleBase::Vector3<double> G);
	std::complex<double> unkdotp_soc_G(const int ik_L, const int ik_R, const int iband_L, const int iband_R, const ModuleBase::ComplexMatrix *evc);
	std::complex<double> unkdotp_soc_G0(const int ik_L, const int ik_R, const int iband_L, const int iband_R, const ModuleBase::ComplexMatrix *evc, const ModuleBase::Vector3<double> G);
	
	// this function just for test the class unkOverlap_pw that is works successful.
	void test_for_unkOverlap_pw();
	
	//std::complex<double> unkdotp_R(int ik_L, int ik_R, int iband_L, int iband_R, ModuleBase::ComplexMatrix *evc);
	//std::complex<double> g00(int ik_R, int ik_L, int ib_L, int ib_R, double x, double y, double z, ModuleBase::ComplexMatrix *evc);

};

#endif
