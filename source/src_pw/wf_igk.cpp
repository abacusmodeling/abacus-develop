#include "wf_igk.h"

#include "../src_pw/global.h"
#include "../module_base/memory.h"
#include "../module_base/timer.h"

WF_igk::WF_igk()
{
}

WF_igk::~WF_igk()
{
	if (GlobalV::test_deconstructor)
	{
		std::cout << " ~WF_igk()" << std::endl;
	}
}


//--------------------------------------------------------
// Compute kinetic energy for each k-point
//--------------------------------------------------------

ModuleBase::Vector3<double> WF_igk::get_1qvec_cartesian(const int ik, const int ig) const
{
	ModuleBase::Vector3<double> qvec = GlobalC::wfcpw->getgpluskcar(ik,ig);

	/*
	if(igk(ik,ig)==0)
	{
		std::cout << " g add = " << GlobalC::sf.gcar << std::endl;
		std::cout << "\n igk = " << igk(ik,ig);
		std::cout << " k=" << GlobalC::kv.kvec_c[ik].x << " " << GlobalC::kv.kvec_c[ik].y << " " <<
	GlobalC::kv.kvec_c[ik].z; std::cout << " g=" << GlobalC::sf.gcar[ this->igk(ik, ig) ].x << " " << GlobalC::sf.gcar[
	this->igk(ik, ig) ].y << " " << GlobalC::sf.gcar[ this->igk(ik, ig) ].z;
	}
	*/

	return qvec;
}

double *WF_igk::get_qvec_cartesian(const int &ik)
{
	double *qmod = new double[GlobalC::kv.ngk[ik]];
	for (int ig = 0; ig < GlobalC::kv.ngk[ik]; ig++)
	{
		// cartesian coordinate
		// modulus, in GlobalC::ucell.tpiba unit.
		const double q2 = GlobalC::wfcpw->getgk2(ik,ig);
		qmod[ig] = GlobalC::ucell.tpiba * sqrt(q2); // sqrt(q2) * GlobalC::ucell.tpiba;
	}
	return qmod;
}

std::complex<double> *WF_igk::get_sk(const int ik, const int it, const int ia, ModulePW::PW_Basis_K* wfc_basis) const
{
	ModuleBase::timer::tick("WF_igk", "get_sk");
	const double arg = (GlobalC::kv.kvec_c[ik] * GlobalC::ucell.atoms[it].tau[ia]) * ModuleBase::TWO_PI;
	const std::complex<double> kphase = std::complex<double>(cos(arg), -sin(arg));
	std::complex<double> *sk = new std::complex<double>[GlobalC::kv.ngk[ik]];
	const int nx = wfc_basis->nx, ny = wfc_basis->ny, nz = wfc_basis->nz;
	for(int igl = 0; igl < GlobalC::kv.ngk[ik]; ++igl)
	{
		const int isz = wfc_basis->getigl2isz(ik,igl);
    	int iz = isz % nz;
    	const int is = isz / nz;
		const int ixy = wfc_basis->is2fftixy[is];
    	int ix = ixy / wfc_basis->fftny;
    	int iy = ixy % wfc_basis->fftny;
		if (ix < int(nx/2) + 1) ix += nx;
        if (iy < int(ny/2) + 1) iy += ny;
        if (iz < int(nz/2) + 1) iz += nz;
		const int iat = GlobalC::ucell.itia2iat(it, ia);
		sk[igl] = kphase * GlobalC::sf.eigts1(iat, ix) * GlobalC::sf.eigts2(iat, iy)
				 * GlobalC::sf.eigts3(iat, iz);
	}
	ModuleBase::timer::tick("WF_igk", "get_sk");
	return sk;
}

std::complex<double> *WF_igk::get_skq(int ik,
									  const int it,
									  const int ia,
									  ModuleBase::Vector3<double> q) // pengfei 2016-11-23
{
	std::complex<double> *skq = new std::complex<double>[GlobalC::kv.ngk[ik]];

	for (int ig = 0; ig < GlobalC::kv.ngk[ik]; ig++)
	{
		ModuleBase::Vector3<double> qkq = GlobalC::wfcpw->getgpluskcar(ik,ig) + q;
		double arg = (qkq * GlobalC::ucell.atoms[it].tau[ia]) * ModuleBase::TWO_PI;
		skq[ig] = std::complex<double>(cos(arg), -sin(arg));
	}

	return skq;
}
