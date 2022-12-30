#include "stress_func.h"
#include "../module_base/math_polyint.h"
#include "../module_base/math_ylmreal.h"
#include "../module_base/timer.h"
#include "global.h"
#include "module_psi/include/device.h"

//calculate the nonlocal pseudopotential stress in PW
template <typename FPTYPE, typename Device>
void Stress_Func<FPTYPE, Device>::stress_nl(ModuleBase::matrix& sigma, const ModuleBase::matrix& wg, const psi::Psi<complex<FPTYPE>, Device>* psi_in)
{
	ModuleBase::TITLE("Stress_Func","stres_nl");
	ModuleBase::timer::tick("Stress_Func","stres_nl");
	
	const int nkb = GlobalC::ppcell.nkb;
	if(nkb == 0) 
	{
		ModuleBase::timer::tick("Stress_Func","stres_nl");
		return;
	}

    this->device = psi::device::get_device_type<Device>(this->ctx);

    // FPTYPE sigmanlc[3][3];
	// for(int l=0;l<3;l++)
	// {
	// 	for(int m=0;m<3;m++)
	// 	{
	// 		sigmanlc[l][m]=0.0;
	// 	}
	// }

	// dbecp: conj( -iG * <Beta(nkb,npw)|psi(nbnd,npw)> )
	// ModuleBase::ComplexMatrix dbecp( GlobalV::NBANDS, nkb );
	// ModuleBase::ComplexMatrix becp( GlobalV::NBANDS, nkb );

	// vkb1: |Beta(nkb,npw)><Beta(nkb,npw)|psi(nbnd,npw)>
	// ModuleBase::ComplexMatrix vkb1( nkb, GlobalC::wf.npwx );
	ModuleBase::ComplexMatrix vkb0[3];
	for(int i=0;i<3;i++){
		vkb0[i].create(nkb, GlobalC::wf.npwx);
	}
	ModuleBase::ComplexMatrix vkb2( nkb, GlobalC::wf.npwx );
    std::complex<FPTYPE> * dbecp = nullptr, * becp = nullptr, * dbecp_noevc = nullptr, * vkb = nullptr, * pvkb0 = nullptr, * vkb1 = nullptr, * pvkb2 = nullptr;
    std::complex<FPTYPE> * _vkb0[3] = {nullptr, nullptr, nullptr};
    resmem_complex_op()(this->ctx, becp, GlobalV::NBANDS * nkb);
    resmem_complex_op()(this->ctx, dbecp, GlobalV::NBANDS * nkb);
    resmem_complex_op()(this->ctx, dbecp_noevc, nkb * GlobalC::wf.npwx);
    resmem_complex_op()(this->ctx, vkb1, nkb * GlobalC::wf.npwx);

    int wg_nc = wg.nc;
    int * atom_nh = nullptr, * atom_na = nullptr, * h_atom_nh = new int[GlobalC::ucell.ntype], * h_atom_na = new int[GlobalC::ucell.ntype];
    for (int ii = 0; ii < GlobalC::ucell.ntype; ii++) {
        h_atom_nh[ii] = GlobalC::ucell.atoms[ii].ncpp.nh;
        h_atom_na[ii] = GlobalC::ucell.atoms[ii].na;
    }
    FPTYPE * stress = nullptr, * sigmanlc = nullptr, * d_wg = nullptr, * deeq = nullptr, * gcar = nullptr, * kvec_c = nullptr, * qvec = nullptr;
    resmem_var_op()(this->ctx, qvec, 3);
    resmem_var_op()(this->ctx, stress, 9);
    setmem_var_op()(this->ctx, stress, 0, 9);
    resmem_var_h_op()(this->cpu_ctx, sigmanlc, 9);
    if (this->device == psi::GpuDevice) {
        deeq = GlobalC::ppcell.d_deeq;
        resmem_var_op()(this->ctx, d_wg, wg.nr * wg.nc);
        resmem_var_op()(this->ctx, gcar, 3 * GlobalC::kv.nks * GlobalC::wfcpw->npwk_max);
        resmem_var_op()(this->ctx, kvec_c, 3 * GlobalC::kv.nks);
        syncmem_var_h2d_op()(this->ctx, this->cpu_ctx, d_wg, wg.c, wg.nr * wg.nc);
        syncmem_var_h2d_op()(this->ctx, this->cpu_ctx, gcar, &GlobalC::wfcpw->gcar[0][0], 3 * GlobalC::kv.nks * GlobalC::wfcpw->npwk_max);
        syncmem_var_h2d_op()(this->ctx, this->cpu_ctx, kvec_c, &GlobalC::wfcpw->kvec_c[0][0], 3 * GlobalC::kv.nks);
        resmem_complex_op()(this->ctx, pvkb2, nkb * GlobalC::wf.npwx);
        resmem_complex_op()(this->ctx, pvkb0, 3 * nkb * GlobalC::wf.npwx);
        for (int ii = 0; ii < 3; ii++) {
            _vkb0[ii] = pvkb0 + ii * nkb * GlobalC::wf.npwx;
        }
        resmem_int_op()(this->ctx, atom_nh, GlobalC::ucell.ntype);
        resmem_int_op()(this->ctx, atom_na, GlobalC::ucell.ntype);
        syncmem_int_h2d_op()(this->ctx, this->cpu_ctx, atom_nh, h_atom_nh, GlobalC::ucell.ntype);
        syncmem_int_h2d_op()(this->ctx, this->cpu_ctx, atom_na, h_atom_na, GlobalC::ucell.ntype);
    }
    else {
        deeq = GlobalC::ppcell.deeq.ptr;
        d_wg = wg.c;
        gcar = &GlobalC::wfcpw->gcar[0][0];
        kvec_c = &GlobalC::wfcpw->kvec_c[0][0];
        atom_nh = h_atom_nh;
        atom_na = h_atom_na;
        for (int ii = 0; ii < 3; ii++) {
            _vkb0[ii] = vkb0[ii].c;
        }
    }

    for (int ik = 0;ik < GlobalC::kv.nks;ik++)
    {   	  
		if (GlobalV::NSPIN==2) GlobalV::CURRENT_SPIN = GlobalC::kv.isk[ik];
		const int npw = GlobalC::kv.ngk[ik];
		// generate vkb
		if (GlobalC::ppcell.nkb > 0) {
            if (this->device == psi::GpuDevice) {
                vkb = GlobalC::ppcell.d_vkb;
            }
            else {
                vkb = GlobalC::ppcell.vkb.c;
            }
			GlobalC::ppcell.getvnl(ctx, ik, vkb);
		}

		// get becp according to wave functions and vkb
		// important here ! becp must set zero!!
		// vkb: Beta(nkb,npw)
		// becp(nkb,nbnd): <Beta(nkb,npw)|psi(nbnd,npw)>
        // becp.zero_out();
		const std::complex<FPTYPE>* ppsi = nullptr;
		if(psi_in!=nullptr)
		{
			ppsi = &(psi_in[0](ik, 0, 0));
		}
		else
		{
			ppsi = &(GlobalC::wf.evc[ik](0, 0));
		}
		char transa = 'C';
        char transb = 'N';
        ///
        ///only occupied band should be calculated.
        ///
        int nbands_occ = GlobalV::NBANDS;
        while(wg(ik, nbands_occ - 1) < ModuleBase::threshold_wg) {
            nbands_occ--;
        }
        int npm = GlobalV::NPOL * nbands_occ;
        gemm_op()(
            this->ctx,
            transa,
            transb,
            nkb,
            npm,
            npw,
            &ModuleBase::ONE,
            vkb,
            GlobalC::wf.npwx,
            ppsi,
            GlobalC::wf.npwx,
            &ModuleBase::ZERO,
            becp,
            nkb);
		//becp calculate is over , now we should broadcast this data.
        if (this->device == psi::GpuDevice) {
            std::complex<FPTYPE> * h_becp = nullptr;
            resmem_complex_h_op()(this->cpu_ctx, h_becp, GlobalV::NBANDS * nkb);
            syncmem_complex_d2h_op()(this->cpu_ctx, this->ctx, h_becp, becp, GlobalV::NBANDS * nkb);
            Parallel_Reduce::reduce_complex_double_pool(becp, GlobalV::NBANDS * nkb);
            syncmem_complex_h2d_op()(this->ctx, this->cpu_ctx, becp, h_becp, GlobalV::NBANDS * nkb);
            delmem_complex_h_op()(this->cpu_ctx, h_becp);
        }
        else {
            Parallel_Reduce::reduce_complex_double_pool(becp, GlobalV::NBANDS * nkb);
        }
		for (int i = 0; i < 3; i++) 
		{
			get_dvnl1(vkb0[i], ik, i);
            if (this->device == psi::GpuDevice) {
                syncmem_complex_h2d_op()(this->ctx, this->cpu_ctx, _vkb0[i], vkb0[i].c, nkb * GlobalC::wf.npwx);
            }
        }
        get_dvnl2(vkb2, ik);
        if (this->device == psi::GpuDevice) {
            syncmem_complex_h2d_op()(this->ctx, this->cpu_ctx, pvkb2, vkb2.c, nkb * GlobalC::wf.npwx);
        }
        else {
            pvkb2 = vkb2.c;
        }

        for (int ipol = 0; ipol < 3; ipol++) 
		{
            for (int jpol = 0; jpol < ipol + 1; jpol++) 
			{
                setmem_complex_op()(this->ctx, vkb1, 0, nkb * GlobalC::wf.npwx);
                setmem_complex_op()(this->ctx, dbecp_noevc, 0, nkb * GlobalC::wf.npwx);
                cal_dbecp_noevc_nl_op()(
                    this->ctx,
                    ipol,
                    jpol,
                    nkb,
                    npw,
                    GlobalC::wf.npwx,
                    ik,
                    GlobalC::ucell.tpiba,
                    gcar,
                    kvec_c,
                    _vkb0[ipol],
                    _vkb0[jpol],
                    vkb,
                    vkb1,
                    pvkb2,
                    dbecp_noevc);
                gemm_op()(
                    this->ctx,
                    transa,
                    transb,
                    nkb,
                    npm,
                    npw,
                    &ModuleBase::ONE,
                    dbecp_noevc,
                    GlobalC::wf.npwx,
                    ppsi,
                    GlobalC::wf.npwx,
                    &ModuleBase::ZERO,
                    dbecp,
                    nkb);
				//              don't need to reduce here, keep
				//              dbecp different in each
				//              processor, and at last sum up
				//              all the forces.
				//              Parallel_Reduce::reduce_complex_double_pool(
				//              dbecp.ptr, dbecp.ndata);

                //              FPTYPE *cf = new
                //              FPTYPE[GlobalC::ucell.nat*3];
                //              ModuleBase::GlobalFunc::ZEROS(cf,
                //              GlobalC::ucell.nat);
                cal_stress_nl_op()(
                    this->ctx,
                    GlobalC::ppcell.multi_proj,
                    ipol,
                    jpol,
                    nkb,
                    nbands_occ,
                    GlobalC::ucell.ntype,
                    GlobalV::CURRENT_SPIN,
                    wg_nc,
                    ik,
                    GlobalC::ppcell.deeq.getBound2(),
                    GlobalC::ppcell.deeq.getBound3(),
                    GlobalC::ppcell.deeq.getBound4(),
                    atom_nh,
                    atom_na,
                    d_wg,
                    deeq,
                    becp,
                    dbecp,
                    stress);

            }//end jpol
		}//end ipol
	}// end ik

    syncmem_var_d2h_op()(this->cpu_ctx, this->ctx, sigmanlc, stress, 9);
	// sum up forcenl from all processors
	for(int l=0;l<3;l++)
	{
		for(int m=0;m<3;m++)
		{
			if(m>l) 
			{
				sigmanlc[l * 3 + m] = sigmanlc[m * 3 + l];
			}
			Parallel_Reduce::reduce_double_all( sigmanlc[l * 3 + m] ); //qianrui fix a bug for kpar > 1
		}
	}

//        Parallel_Reduce::reduce_double_all(sigmanl.c, sigmanl.nr * sigmanl.nc);
        
	for (int ipol = 0; ipol<3; ipol++)
	{
		for(int jpol = 0; jpol < 3; jpol++)
		{
			sigmanlc[ipol * 3 + jpol] *= 1.0 / GlobalC::ucell.omega;
		}
	}
	
	for (int ipol = 0; ipol<3; ipol++)
	{
		for(int jpol = 0; jpol < 3; jpol++)
		{
			sigma(ipol,jpol) = sigmanlc[ipol * 3 + jpol] ;
		}
	}
	//do symmetry
	if(ModuleSymmetry::Symmetry::symm_flag == 1)
	{
		GlobalC::symm.stress_symmetry(sigma, GlobalC::ucell);
	}//end symmetry

    delete [] h_atom_nh;
    delete [] h_atom_na;
    delmem_var_op()(this->ctx, qvec);
    delmem_var_op()(this->ctx, stress);
    delmem_complex_op()(this->ctx, becp);
    delmem_complex_op()(this->ctx, vkb1);
    delmem_complex_op()(this->ctx, pvkb0);
    delmem_complex_op()(this->ctx, dbecp);
    delmem_complex_op()(this->ctx, dbecp_noevc);
    if (this->device == psi::GpuDevice) {
        delmem_var_op()(this->ctx, d_wg);
        delmem_var_op()(this->ctx, gcar);
        delmem_var_op()(this->ctx, kvec_c);
        delmem_var_h_op()(this->cpu_ctx, sigmanlc);
        delmem_int_op()(this->ctx, atom_nh);
        delmem_int_op()(this->ctx, atom_na);
        delmem_complex_op()(this->ctx, pvkb2);
    }
	//  this->print(GlobalV::ofs_running, "nonlocal stress", stresnl);
	ModuleBase::timer::tick("Stress_Func","stres_nl");
}

template <typename FPTYPE, typename Device>
void Stress_Func<FPTYPE, Device>::get_dvnl1
(
	ModuleBase::ComplexMatrix &vkb,
	const int ik,
	const int ipol
)
{
	if(GlobalV::test_pp) ModuleBase::TITLE("Stress_Func","get_dvnl1");

	const int lmaxkb = GlobalC::ppcell.lmaxkb;
	if(lmaxkb < 0)
	{
		return;
	}

	const int npw = GlobalC::kv.ngk[ik];
	const int nhm = GlobalC::ppcell.nhm;
	int ig, ia, nb, ih;
	ModuleBase::matrix vkb1(nhm, npw);
	vkb1.zero_out();
	FPTYPE *vq = new FPTYPE[npw];
	const int x1= (lmaxkb + 1)*(lmaxkb + 1);

	ModuleBase::matrix dylm(x1, npw);
	ModuleBase::Vector3<FPTYPE> *gk = new ModuleBase::Vector3<FPTYPE>[npw];
	for (ig = 0;ig < npw;ig++)
	{
		gk[ig] = GlobalC::wf.get_1qvec_cartesian(ik, ig);
	}
			   
	dylmr2(x1, npw, gk, dylm, ipol);

	int jkb = 0;
	for(int it = 0;it < GlobalC::ucell.ntype;it++)
	{
		if(GlobalV::test_pp>1) ModuleBase::GlobalFunc::OUT("it",it);
		// calculate beta in G-space using an interpolation table
		const int nbeta = GlobalC::ucell.atoms[it].ncpp.nbeta;
		const int nh = GlobalC::ucell.atoms[it].ncpp.nh;

		if(GlobalV::test_pp>1) ModuleBase::GlobalFunc::OUT("nbeta",nbeta);

		for (nb = 0;nb < nbeta;nb++)
		{
			if(GlobalV::test_pp>1) ModuleBase::GlobalFunc::OUT("ib",nb);
			for (ig = 0;ig < npw;ig++)
			{
				const FPTYPE gnorm = gk[ig].norm() * GlobalC::ucell.tpiba;

				//cout << "\n gk[ig] = " << gk[ig].x << " " << gk[ig].y << " " << gk[ig].z;
				//cout << "\n gk.norm = " << gnorm;

				vq [ig] = ModuleBase::PolyInt::Polynomial_Interpolation(
						GlobalC::ppcell.tab, it, nb, GlobalV::NQX, GlobalV::DQ, gnorm );

			} // enddo

			// add spherical harmonic part
			for (ih = 0;ih < nh;ih++)
			{
				if (nb == GlobalC::ppcell.indv(it, ih))
				{
					const int lm = static_cast<int>( GlobalC::ppcell.nhtolm(it, ih) );
					for (ig = 0;ig < npw;ig++)
					{
						vkb1(ih, ig) = dylm(lm, ig) * vq [ig];
					}

				}

			} // end ih

		} // end nbeta

		// vkb1 contains all betas including angular part for type nt
		// now add the structure factor and factor (-i)^l
		for (ia=0; ia<GlobalC::ucell.atoms[it].na; ia++)
		{
			std::complex<FPTYPE> *sk = GlobalC::wf.get_sk(ik, it, ia,GlobalC::wfcpw);
			for (ih = 0;ih < nh;ih++)
			{
				std::complex<FPTYPE> pref = pow( ModuleBase::NEG_IMAG_UNIT, GlobalC::ppcell.nhtol(it, ih));      //?
				for (ig = 0;ig < npw;ig++)
				{
					vkb(jkb, ig) = vkb1(ih, ig) * sk [ig] * pref;
				}
				++jkb;
			} // end ih
		delete [] sk;
		} // end ia
	} // enddo
	delete [] gk;
	delete [] vq;
	return;
}//end get_dvnl1

template <typename FPTYPE, typename Device>
void Stress_Func<FPTYPE, Device>::get_dvnl2(ModuleBase::ComplexMatrix &vkb,
		const int ik)
{
	if(GlobalV::test_pp) ModuleBase::TITLE("Stress","get_dvnl2");
//	ModuleBase::timer::tick("Stress","get_dvnl2");

	const int lmaxkb = GlobalC::ppcell.lmaxkb;
	if(lmaxkb < 0)
	{
		return;
	}

	const int npw = GlobalC::kv.ngk[ik];
	const int nhm = GlobalC::ppcell.nhm;
	int ig, ia, nb, ih;
	ModuleBase::matrix vkb1(nhm, npw);
	FPTYPE *vq = new FPTYPE[npw];
	const int x1= (lmaxkb + 1)*(lmaxkb + 1);

	ModuleBase::matrix ylm(x1, npw);
	ModuleBase::Vector3<FPTYPE> *gk = new ModuleBase::Vector3<FPTYPE>[npw];
	for (ig = 0;ig < npw;ig++)
	{
		gk[ig] = GlobalC::wf.get_1qvec_cartesian(ik, ig);
	}
	ModuleBase::YlmReal::Ylm_Real(x1, npw, gk, ylm);

	int jkb = 0;
	for(int it = 0;it < GlobalC::ucell.ntype;it++)
	{
		if(GlobalV::test_pp>1) ModuleBase::GlobalFunc::OUT("it",it);
		// calculate beta in G-space using an interpolation table
		const int nbeta = GlobalC::ucell.atoms[it].ncpp.nbeta;
		const int nh = GlobalC::ucell.atoms[it].ncpp.nh;

		if(GlobalV::test_pp>1) ModuleBase::GlobalFunc::OUT("nbeta",nbeta);

		for (nb = 0;nb < nbeta;nb++)
		{
			if(GlobalV::test_pp>1) ModuleBase::GlobalFunc::OUT("ib",nb);
			for (ig = 0;ig < npw;ig++)
			{
				const FPTYPE gnorm = gk[ig].norm() * GlobalC::ucell.tpiba;
	//cout << "\n gk[ig] = " << gk[ig].x << " " << gk[ig].y << " " << gk[ig].z;
	//cout << "\n gk.norm = " << gnorm;
				vq [ig] = Polynomial_Interpolation_nl(
						GlobalC::ppcell.tab, it, nb, GlobalV::DQ, gnorm );

			} // enddo

							// add spherical harmonic part
			for (ih = 0;ih < nh;ih++)
			{
				if (nb == GlobalC::ppcell.indv(it, ih))
				{
					const int lm = static_cast<int>( GlobalC::ppcell.nhtolm(it, ih) );
					for (ig = 0;ig < npw;ig++)
					{
						vkb1(ih, ig) = ylm(lm, ig) * vq [ig];
					}
				}
			} // end ih
		} // end nbeta

		// vkb1 contains all betas including angular part for type nt
		// now add the structure factor and factor (-i)^l
		for (ia=0; ia<GlobalC::ucell.atoms[it].na; ia++)
		{
			std::complex<FPTYPE> *sk = GlobalC::wf.get_sk(ik, it, ia,GlobalC::wfcpw);
			for (ih = 0;ih < nh;ih++)
			{
				std::complex<FPTYPE> pref = pow( ModuleBase::NEG_IMAG_UNIT, GlobalC::ppcell.nhtol(it, ih));      //?
				for (ig = 0;ig < npw;ig++)
				{
					vkb(jkb, ig) = vkb1(ih, ig) * sk [ig] * pref;
				}
			++jkb;
			} // end ih
			delete [] sk;
		} // end ia
	}	 // enddo

	delete [] gk;
	delete [] vq;
//	ModuleBase::timer::tick("Stress","get_dvnl2");

	return;
}


template <typename FPTYPE, typename Device>
FPTYPE Stress_Func<FPTYPE, Device>::Polynomial_Interpolation_nl
(
    const ModuleBase::realArray &table,
    const int &dim1,
    const int &dim2,
    const FPTYPE &table_interval,
    const FPTYPE &x                             // input value
)
{

	assert(table_interval>0.0);
	const FPTYPE position = x  / table_interval;
	const int iq = static_cast<int>(position);

	const FPTYPE x0 = position - static_cast<FPTYPE>(iq);
	const FPTYPE x1 = 1.0 - x0;
	const FPTYPE x2 = 2.0 - x0;
	const FPTYPE x3 = 3.0 - x0;
	const FPTYPE y=
			( table(dim1, dim2, iq)   * (-x2*x3-x1*x3-x1*x2) / 6.0 +
			table(dim1, dim2, iq+1) * (+x2*x3-x0*x3-x0*x2) / 2.0 -
			table(dim1, dim2, iq+2) * (+x1*x3-x0*x3-x0*x1) / 2.0 +
			table(dim1, dim2, iq+3) * (+x1*x2-x0*x2-x0*x1) / 6.0 )/table_interval ;


	return y;
}

template <typename FPTYPE, typename Device>
void Stress_Func<FPTYPE, Device>::dylmr2 (
	const int nylm,
	const int ngy,
	ModuleBase::Vector3<FPTYPE> *gk,
	ModuleBase::matrix &dylm,
	const int ipol)
{
  //-----------------------------------------------------------------------
  //
  //     compute \partial Y_lm(G) \over \partial (G)_ipol
  //     using simple numerical derivation (SdG)
  //     The spherical harmonics are calculated in ylmr2
  //
  //int nylm, ngy, ipol;
  // number of spherical harmonics
  // the number of g vectors to compute
  // desired polarization
  //FPTYPE g (3, ngy), gg (ngy), dylm (ngy, nylm)
  // the coordinates of g vectors
  // the moduli of g vectors
  // the spherical harmonics derivatives
  //
	int ig, lm;
	// counter on g vectors
	// counter on l,m component

	const FPTYPE delta = 1e-6;
	FPTYPE *dg, *dgi;

	ModuleBase::matrix ylmaux;
	// dg is the finite increment for numerical derivation:
	// dg = delta |G| = delta * sqrt(gg)
	// dgi= 1 /(delta * sqrt(gg))
	// gx = g +/- dg


	ModuleBase::Vector3<FPTYPE> *gx = new ModuleBase::Vector3<FPTYPE> [ngy];
	 

	dg = new FPTYPE [ngy];
	dgi = new FPTYPE [ngy];

	ylmaux.create (nylm, ngy);

	dylm.zero_out();
	ylmaux.zero_out();

	for( ig = 0;ig< ngy;ig++){
		gx[ig] = gk[ig];
	}
	//$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ig)
	for( ig = 0;ig< ngy;ig++){
		dg [ig] = delta * gx[ig].norm() ;
		if (gx[ig].norm2() > 1e-9) {
			dgi [ig] = 1.0 / dg [ig];
		}
		else{
			dgi [ig] = 0.0;
		}
	}
	//$OMP END PARALLEL DO

	//$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ig)
	for( ig = 0;ig< ngy;ig++)
	{
		if(ipol==0)
			gx [ig].x = gk[ ig].x + dg [ig];
		else if(ipol==1)
			gx [ig].y = gk [ ig].y + dg [ig];
		else if(ipol==2)
			gx [ig].z = gk [ ig].z + dg [ig];
	}
	//$OMP END PARALLEL DO

	ModuleBase::YlmReal::Ylm_Real(nylm, ngy, gx, dylm);
	//$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ig)
	for(ig = 0;ig< ngy;ig++)
	{
		if(ipol==0)
			gx [ig].x = gk [ ig].x - dg [ig];
		else if(ipol==1)
			gx [ig].y = gk [ ig].y - dg [ig];
		else if(ipol==2)
			gx [ig].z = gk [ ig].z - dg [ig];
	}
	//$OMP END PARALLEL DO

	ModuleBase::YlmReal::Ylm_Real(nylm, ngy, gx, ylmaux);


	//  zaxpy ( - 1.0, ylmaux, 1, dylm, 1);
	for( lm = 0;lm< nylm;lm++)
	{
		for(ig = 0;ig< ngy;ig++)
		{
			dylm (lm,ig) = dylm(lm,ig) - ylmaux(lm,ig);
		}
	}


	for( lm = 0;lm< nylm;lm++)
	{
		//$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ig)
		for(ig = 0;ig< ngy;ig++)
		{
			dylm (lm,ig) = dylm(lm,ig) * 0.5 * dgi [ig];
		}
		//$OMP END PARALLEL DO
	}
	delete[] gx;
	delete[] dg;
	delete[] dgi;

	return;
}

template class Stress_Func<double, psi::DEVICE_CPU>;
#if ((defined __CUDA) || (defined __ROCM))
template class Stress_Func<double, psi::DEVICE_GPU>;
#endif