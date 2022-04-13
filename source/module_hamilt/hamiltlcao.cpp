#include "hamiltlcao.h"
#include "module_base/global_variable.h"
#include "module_base/timer.h"

namespace ModuleHamilt
{

template<typename T, typename T1>
void HamiltLCAO<T,T1>::constructFixedReal()
{
    //update the flag
    if(HamiltLCAO<T,T1>::isFixedDone)return;
    else HamiltLCAO<T,T1>::isFixedDone=true;

    this->fixedRealM.resize(1, T1(0));

    ModuleBase::TITLE("HamiltLCAO","constructFixedReal"); 
    ModuleBase::timer::tick("HamiltLCAO","constructFixedReal"); 

	// (9) compute S, T, Vnl, Vna matrix.
    this->genH.calculate_S_no(this->fixedRealM.getS());

    if(GlobalV::VNL_IN_H)
    {
        this->genH.calculate_NL_no(this->fixedRealM.getH());
    }

    if(GlobalV::T_IN_H)
    {
        this->genH.calculate_T_no(this->fixedRealM.getH());
    }

#ifdef __DEEPKS
    //for each ionic step, the overlap <psi|alpha> must be rebuilt
    //since it depends on ionic positions
    if (GlobalV::deepks_setorb)
    {
        const Parallel_Orbitals* pv = this->UHM->LM->ParaV;
        //build and save <psi(0)|alpha(R)> at beginning
        GlobalC::ld.build_psialpha(GlobalV::CAL_FORCE,
			GlobalC::ucell,
			GlobalC::ORB,
			GlobalC::GridD,
			pv->trace_loc_row,
			pv->trace_loc_col,
			GlobalC::UOT);

		if(GlobalV::deepks_out_unittest)
		{
			GlobalC::ld.check_psialpha(GlobalV::CAL_FORCE,
					GlobalC::ucell,
					GlobalC::ORB,
					GlobalC::GridD,
					pv->trace_loc_row,
					pv->trace_loc_col,
					GlobalC::UOT);
		}
    }
#endif

    ModuleBase::timer::tick("HamiltLCAO","constructFixedReal"); 
	return;

    return;
}

template<typename T, typename T1>
void HamiltLCAO<T,T1>::constructUpdateReal()
{
    return;
}

template<typename T, typename T1>
void HamiltLCAO<T,T1>::hk_fixed_mock(const int ik)
{
    this->kM.resize(1, T(0));
	//folding_fixedH() should be refactored to there, but now deepks code in it
	return;
}

template<typename T, typename T1>
void HamiltLCAO<T, T1>::hk_update_mock(const int ik)
{
    //update_Hk and update_Hgamma should be refactored to there
    return;
}

template<typename T, typename T1>
void HamiltLCAO<T, T1>::getMatrix(MatrixBlock<T> hk_in, MatrixBlock<T> sk_in)const
{
    hk_in.p = (T*)kM.hloc.data();
    sk_in.p = (T*)kM.sloc.data();
}


// case for nspin==4
/*void HamiltLCAO<std::complex<double>, std::complex<double>>::matrix(MatrixBlock<std::complex<double>> hk_in, MatrixBlock<std::complex<double>> sk_in)const
{
    this->getMatrix(hk_in, sk_in);
}*/
// case for nspin<4, multi-k-points 
void HamiltLCAO<std::complex<double>, double>::matrix(MatrixBlock<std::complex<double>> hk_in, MatrixBlock<std::complex<double>> sk_in)const
{
    this->getMatrix(hk_in, sk_in);
}

// case for nspin<4, gamma_only
void HamiltLCAO<double, double>::matrix(MatrixBlock<double> hk_in, MatrixBlock<double> sk_in)const
{
    this->getMatrix(hk_in, sk_in);
}

// nspin==4 case not supported yet
/*void HamiltLCAO<std::complex<double>, std::complex<double>>::constructHamilt(const int iter, const MatrixBlock<double> rho)
{
    this->constructFixedReal();
    this->constructUpdateReal();
}*/
void HamiltLCAO<std::complex<double>, double>::constructHamilt(const int iter, const MatrixBlock<double> rho)
{
    this->constructFixedReal();
    this->constructUpdateReal();
}
void HamiltLCAO<double, double>::constructHamilt(const int iter, const MatrixBlock<double> rho)
{
    this->constructFixedReal();
    this->constructUpdateReal();
}

/*void HamiltLCAO<std::complex<double>, std::complex<double>>::updateHk(const int ik)
{
    this->hk_fixed_mock(ik);
    this->hk_update_mock(ik);
};*/
void HamiltLCAO<double, std::complex<double>>::updateHk(const int ik)
{
    this->hk_fixed_mock(ik);
    this->hk_update_mock(ik);
};
void HamiltLCAO<double, double>::updateHk(const int ik)
{
    this->hk_fixed_mock(ik);
    this->hk_update_mock(ik);
};

}