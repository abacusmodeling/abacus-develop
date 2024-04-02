#include "paw_cell.h"
#include "module_base/tool_title.h"
#include "module_base/tool_quit.h"
#ifdef __MPI
#include "module_base/parallel_reduce.h"
#endif

namespace GlobalC
{
    Paw_Cell paw_cell;
}

void Paw_Cell::init_paw_cell(
    const double ecutwfc_in, const double cell_factor_in,
    const double omega_in,
    const int nat_in, const int ntyp_in,
    const int * atom_type_in, const double ** atom_coord_in,
    const std::vector<std::string> & filename_list_in)
{
    ModuleBase::TITLE("Paw_Element","init_paw_cell");

    this -> nat = nat_in;
    this -> ntyp = ntyp_in;
    this -> omega = omega_in;

    atom_coord.resize(nat);
    atom_type.resize(nat);
    for(int iat = 0; iat < nat; iat ++)
    {
        atom_coord[iat].resize(3);
        for(int i = 0; i < 3; i ++)
        {
            atom_coord[iat][i] = atom_coord_in[iat][i];
        }

        atom_type[iat] = atom_type_in[iat];
    }

    nat_type.resize(ntyp);
    atom_map.resize(nat);

    // Not sure if these will be used, but it doesn't
    // hurt to prepare them I suppose
    // Doesn't take much time, also not much space
    int ind = 0;
    for(int it = 0; it < ntyp; it ++)
    {
        nat_type[it] = 0;
        for(int iat = 0; iat < nat; iat ++)
        {
            if(atom_type[iat] == it)
            {
                nat_type[it] ++;
                atom_map[ind] = iat;
                ind ++;
            }
        }
    }

    paw_element_list.resize(ntyp);
    assert(filename_list_in.size() == ntyp);
    for(int ityp = 0; ityp < ntyp; ityp ++)
    {
        paw_element_list[ityp].init_paw_element(ecutwfc_in, cell_factor_in);
        paw_element_list[ityp].read_paw_xml(filename_list_in[ityp]);
    }

    this->map_paw_proj();

    paw_atom_list.resize(nat);
    for(int iat = 0; iat < nat; iat ++)
    {
        int it = atom_type[iat];
        int nproj = paw_element_list[it].get_mstates();
        paw_atom_list[iat].init_paw_atom(nproj);
    }

    this -> init_rhoij();
    //this -> init_mix_dij();
}

void Paw_Cell::init_mix_dij()
{
    first_iter = true;
    count = 0;

    if(GlobalV::RANK_IN_POOL == 0)
    {
        dij_save.resize(natom);
        for(int iat = 0; iat < natom; iat ++)
        {
            const int it = atom_type[iat];
            const int nproj = paw_element_list[it].get_mstates();
            const int size_dij = nproj * (nproj+1) / 2;
            dij_save[iat].resize(size_dij * nspden);
            for(int i = 0; i < size_dij * nspden; i ++)
            {
                dij_save[iat][i] = 0.0;
            }
        }
    }
}

void Paw_Cell::init_rhoij()
{
    ModuleBase::TITLE("Paw_Cell","init_rhoij");

    for(int iat = 0; iat < nat; iat ++)
    {
        const int it = atom_type[iat];
        const int nproj = paw_element_list[it].get_mstates();

        const int size_rhoij = nproj * (nproj + 1) / 2;

        std::vector<double> mstate_occ = paw_element_list[it].get_mstate_occ();

        std::vector<double> rhoij_in;
        rhoij_in.resize(size_rhoij);
        for(int i = 0; i < size_rhoij; i ++)
        {
            rhoij_in[i] = 0.0;
        }

        for(int iproj = 0; iproj < nproj; iproj ++)
        {
            int i0 = iproj * (iproj + 1) / 2;
            rhoij_in[i0 + iproj] = mstate_occ[iproj] / GlobalV::NSPIN;
        }

        paw_atom_list[iat].set_rhoij(rhoij_in);
    }
}

void Paw_Cell::set_eigts(const int nx_in, const int ny_in, const int nz_in,
        const std::complex<double> * eigts1_in,
        const std::complex<double> * eigts2_in,
        const std::complex<double> * eigts3_in)
{
    this -> nx = nx_in;
    this -> ny = ny_in;
    this -> nz = nz_in;

    eigts1.resize(nat);
    eigts2.resize(nat);
    eigts3.resize(nat);

    for(int iat = 0; iat < nat; iat ++)
    {
        eigts1[iat].resize(2*nx+1);
        eigts2[iat].resize(2*ny+1);
        eigts3[iat].resize(2*nz+1);

        for(int i = 0; i < 2*nx+1; i ++)
        {
            eigts1[iat][i] = eigts1_in[ iat * (2*nx+1) + i];
        }

        for(int i = 0; i < 2*ny+1; i ++)
        {
            eigts2[iat][i] = eigts2_in[ iat * (2*ny+1) + i];
        }

        for(int i = 0; i < 2*nz+1; i ++)
        {
            eigts3[iat][i] = eigts3_in[ iat * (2*nz+1) + i];
        }
    }
}

// exp(-i(k+G)R_I) = exp(-ikR_I) exp(-iG_xR_Ix) exp(-iG_yR_Iy) exp(-iG_zR_Iz)
void Paw_Cell::set_paw_k(
    const int npw_in, const int npwx_in, const double * kpt,
    const int * ig_to_ix, const int * ig_to_iy, const int * ig_to_iz,
    const double ** kpg, const double tpiba, const double ** gcar)
{
    ModuleBase::TITLE("Paw_Element","set_paw_k");

    const double pi = 3.141592653589793238462643383279502884197;
    const double twopi = 2.0 * pi;

    this -> npw = npw_in;
    this -> npwx = npwx_in;

    struc_fact.resize(nat);
    for(int iat = 0; iat < nat; iat ++)
    {
        double arg = 0.0;
        for(int i = 0; i < 3; i ++)
        {
            arg += atom_coord[iat][i] * kpt[i];
        }
        arg *= twopi;
        const std::complex<double> kphase = std::complex<double>(cos(arg), -sin(arg));

        struc_fact[iat].resize(npw);

        for(int ipw = 0; ipw < npw; ipw ++)
        {
            int ix = ig_to_ix[ipw];
            int iy = ig_to_iy[ipw];
            int iz = ig_to_iz[ipw];

            struc_fact[iat][ipw] = kphase * eigts1[iat][ix] * eigts2[iat][iy] * eigts3[iat][iz];
        }
    }

    this -> set_ylm(npw, kpg);
    
    gnorm.resize(npw);
    for(int ipw = 0; ipw < npw; ipw ++)
    {
        gnorm[ipw] = std::sqrt(kpg[ipw][0]*kpg[ipw][0] + kpg[ipw][1]*kpg[ipw][1] + kpg[ipw][2]*kpg[ipw][2]) * tpiba;
    }

    std::complex<double> i_cplx(0.0,1.0);
    // ig : i(G)
    if(GlobalV::CAL_FORCE || GlobalV::CAL_STRESS)
    {
        ig.resize(npw);
        for(int ipw = 0; ipw < npw; ipw ++)
        {
            ig[ipw].resize(3);
            for(int i = 0; i < 3; i ++)
            {
                ig[ipw][i] = gcar[ipw][i] * tpiba * i_cplx;
            }
        }
    }
}

void Paw_Cell::set_isk(const int nk, const int * isk_in)
{
    isk.resize(nk);
    for(int ik = 0; ik < nk; ik ++)
    {
        isk[ik] = isk_in[ik];
    }
}

void Paw_Cell::set_currentk(const int ik)
{
    current_k = ik;
    current_spin = isk[ik];
}

void Paw_Cell::map_paw_proj()
{
    ModuleBase::TITLE("Paw_Element","map_paw_proj");

    nproj_tot = 0;
    
    for(int ia = 0; ia < nat; ia ++)
    {
        int it = atom_type[ia];
        nproj_tot += paw_element_list[it].get_mstates();
    }

    // Not sure if all of them will be used, but it doesn't
    // hurt to prepare them I suppose
    // Doesn't take much time, also not much space

    iprj_to_ia.resize(nproj_tot);
    iprj_to_im.resize(nproj_tot);
    iprj_to_il.resize(nproj_tot);
    iprj_to_l.resize(nproj_tot);
    iprj_to_m.resize(nproj_tot);
    start_iprj_ia.resize(nat);

    int iproj = 0;
    for(int ia = 0; ia < nat; ia ++)
    {
        start_iprj_ia[ia] = iproj;
        int it = atom_type[ia];
        int mstates = paw_element_list[it].get_mstates();
        std::vector<int> im_to_istate = paw_element_list[it].get_im_to_istate();
        std::vector<int> lstate = paw_element_list[it].get_lstate();
        std::vector<int> mstate = paw_element_list[it].get_mstate();

        for(int im = 0; im < mstates; im ++)
        {
            iprj_to_ia[iproj] = ia;
            iprj_to_im[iproj] = im;
            iprj_to_il[iproj] = im_to_istate[im];
            iprj_to_l [iproj] = lstate[im_to_istate[im]];
            iprj_to_m [iproj] = mstate[im];
            iproj ++;
        }
    }

    lmax = 0;
    for(int it = 0; it < ntyp; it ++)
    {
        lmax = std::max(lmax, paw_element_list[it].get_lmax());
    }

    assert(iproj == nproj_tot);
}

void Paw_Cell::set_ylm(const int npw_in, const double ** kpg)
{
    ylm_k.resize(npw_in);
    for(int ipw = 0; ipw < npw_in; ipw ++)
    {
        ylm_k[ipw] = calc_ylm(lmax, kpg[ipw]);
    }
}

std::vector<double> Paw_Cell::calc_ylm(const int lmax, const double * r)
{
    // there are altogether (l+1)^2 spherical harmonics
    // up to angular momentum l
    int size_ylm = (lmax + 1) * (lmax + 1);

    std::vector<double> ylm;
    ylm.resize(size_ylm);

    double xx = r[0];
    double yy = r[1];
    double zz = r[2];
    double rr = std::sqrt(xx*xx+yy*yy+zz*zz);

    double tol = 1e-10;

    double cphi, sphi, ctheta, stheta;
    std::vector<std::complex<double>> phase;
    phase.resize(lmax+1);

    const double pi = 3.141592653589793238462643383279502884197;
    const double fourpi = 4.0 * pi;

    // set zero
    for(int i = 0; i < size_ylm; i++)
    {
        ylm[i] = 0;
    }

    // l = 0
    ylm[0] = 1.0/std::sqrt(fourpi);

    if(rr>tol)
    {
        cphi = 1.0;
        sphi = 0.0;
        ctheta = zz/rr;
        stheta = std::sqrt(std::abs((1.0-ctheta)*(1.0+ctheta)));

        if(stheta > tol)
        {
            cphi = xx/(rr*stheta); // cos (phi)
            sphi = yy/(rr*stheta); // sin (phi)
        }
        
        for(int m = 1; m < lmax + 1; m ++)
        {
            std::complex<double> tmp(cphi,sphi);
            phase[m] = std::pow(tmp,m); //exp(i m phi)
        }

        for(int l = 1; l < lmax + 1; l++)
        {
            const int l0 = l*l + l;
            double fact  = 1.0/double(l*(l+1));
            const double ylmcst = std::sqrt(double(2*l+1)/fourpi);
            
            // m = 0
            ylm[l0] = ylmcst*Paw_Cell::ass_leg_pol(l,0,ctheta);
            
            // m > 0
            double onem = 1.0;
            for(int m = 1; m < l + 1; m ++)
            {
                onem = - onem; // the factor (-1)^m
                double work1 = ylmcst*std::sqrt(fact)*onem
                    *Paw_Cell::ass_leg_pol(l,m,ctheta)*std::sqrt(2.0);
                ylm[l0+m] = work1 * phase[m].real();
                ylm[l0-m] = work1 * phase[m].imag();
                if(m != l) fact = fact/double((l+m+1)*(l-m));
            }
        }
    }

    return ylm;

}

double Paw_Cell::ass_leg_pol(const int l, const int m, const double arg)
{

    double x = arg;
    if(m < 0 || m > l || std::abs(x) > 1.0)
    {
        if(std::abs(x) > 1.0 + 1e-10)
        {
            ModuleBase::WARNING_QUIT("Paw_Cell","bad argument l, m, or x");
        }
        x = 1.0;
    }

    double polmm = 1.0;
    if(m > 0)
    {
        double sqrx = std::sqrt(std::abs((1.0-x) * (1.0+x)));
        for(int i = 0; i < m; i++)
        {
            polmm = polmm * (1.0 - 2.0 * double(i+1)) * sqrx;
        }
    }

    if(l == m)
    {
        return polmm;
    }
    else
    {
        double tmp1 = x * (2.0 * double(m) + 1.0) * polmm;
        if(l == m+1)
        {
            return tmp1;
        }
        else
        {

            if (m+1<l)
            {
                double pll = 0.0;
                for(int ll = m+2; ll < l+1; ll ++)
                {
                    pll = (x*(2.0*ll-1.0)*tmp1-(ll+m-1.0)*polmm)/double(ll-m);
                    polmm = tmp1;
                    tmp1 = pll;
                }
                return pll;
            } else {
                ModuleBase::WARNING_QUIT("Paw_Cell","l is too small");
                return 0.0;
            }
        }
    }
}

void Paw_Cell::get_vkb()
{
    ModuleBase::TITLE("Paw_Cell","get_vkb");

    // Not sure if the loop structure is the best
    // may need to improve the efficiency later

    // vkb = ptilde(G)
    // = [\int f(r)r^2j_l(r)dr] * [(-i)^l] * [ylm(\hat{G})] * [exp(-GR_I)]
    // = ptilde * fact * ylm * sk (in the code below)
    vkb.resize(nproj_tot);
    for(int iproj = 0; iproj < nproj_tot; iproj ++)
    {
        vkb[iproj].resize(npw);
    }

    const std::complex<double> neg_i = std::complex<double>(0.0,-1.0);

    for(int ipw = 0; ipw < npw; ipw ++)
    {
        int ind_at = 0;
        for(int it = 0; it < ntyp; it ++)
        {
            const int nstates = paw_element_list[it].get_nstates(); // number of l channels
            const int nproj = paw_element_list[it].get_mstates(); // number of (l,m) channels

            std::vector<double> ptilde;
            ptilde.resize(nproj);

            int iproj = 0;
            for(int istate = 0; istate < nstates; istate ++)
            {
                const double tmp = paw_element_list[it].get_ptilde(istate,gnorm[ipw],omega);
                const int l = paw_element_list[it].get_lstate()[istate];

                for(int im = 0; im < 2*l+1; im ++)
                {
                    ptilde[iproj] = tmp;
                    iproj ++;
                }
            }
            assert(iproj == nproj);

            const int na = nat_type[it];
            for(int ia = 0; ia < na; ia ++)
            {
                const int iat = atom_map[ind_at+ia];
                const int proj_start = start_iprj_ia[iat];

                for(int ip = 0; ip < nproj; ip ++)
                {
                    const int il = paw_element_list[it].get_im_to_istate()[ip];
                    const int l  = paw_element_list[it].get_lstate()[il];
                    const int m  = paw_element_list[it].get_mstate()[ip];

                    // I put these codes in comment here as a reminder of
                    // how this part is treated differently from nonlocal PP
                    // First of all, the order of m quantum number
                    // in nonlocal PP (calculated using ylm_real) is 0,1,-1,2,-2,...
                    // but -l,-l+1,-l+2,...,l here
                    // Secondly, there is a factor (-1)^m which is present here
                    // but not in ylm_real

                    /*
                    int m1 = 0;
                    if(l==1 && m == -1) m1 = 0;
                    if(l==1 && m == 0 ) m1 = 1;
                    if(l==1 && m == 1 ) m1 = -1;
                    */

                    //const int lm_ind = l*l + l + m1;
                    const int lm_ind = l*l + l + m;

                    const std::complex<double> fact = std::pow(neg_i,l);

                    //const double ylm = ylm_k[ipw][lm_ind] * std::pow(-1.0,m1);
                    const double ylm = ylm_k[ipw][lm_ind];

                    const std::complex<double> sk = struc_fact[iat][ipw];

                    vkb[ip+proj_start][ipw] = fact * ptilde[ip] * ylm * sk;
                }
            }
            
            ind_at += na;
        }
        assert(ind_at == nat);
    }
}

void Paw_Cell::reset_rhoij()
{
    for(int iat = 0; iat < nat; iat ++)
    {
        paw_atom_list[iat].reset_rhoij();
    }
}

void Paw_Cell::accumulate_rhoij(const std::complex<double> * psi, const double weight)
{
    ModuleBase::TITLE("Paw_Cell","accumulate_rhoij");

    for(int iat = 0; iat < nat; iat ++)
    {
        // ca : <ptilde(G)|psi(G)>
        // = \sum_G [\int f(r)r^2j_l(r)dr] * [(-i)^l] * [ylm(\hat{G})] * [exp(-GR_I)] *psi(G)
        // = \sum_ipw ptilde * fact * ylm * sk * psi (in the code below)
        // This is what is called 'becp' in nonlocal pp
        std::vector<std::complex<double>> ca;

        const int it = atom_type[iat];
        const int nproj = paw_element_list[it].get_mstates();
        const int proj_start = start_iprj_ia[iat];

        ca.resize(nproj);

        for(int iproj = 0; iproj < nproj; iproj ++)
        {
            ca[iproj] = 0.0;
            
            // consider use blas subroutine for this part later
            for(int ipw = 0; ipw < npw; ipw ++)
            {
                ca[iproj] += psi[ipw] * std::conj(vkb[iproj+proj_start][ipw]);
            }
        }

#ifdef __MPI
        Parallel_Reduce::reduce_pool(ca.data(), nproj);

        if(GlobalV::RANK_IN_POOL == 0)
        {
            paw_atom_list[iat].set_ca(ca, weight);
            paw_atom_list[iat].accumulate_rhoij(current_spin);
        }
#else
        paw_atom_list[iat].set_ca(ca, weight);
        paw_atom_list[iat].accumulate_rhoij(current_spin);
#endif
    }
}

/*
std::vector<std::vector<std::vector<double>>> Paw_Cell::get_rhoij()
{
    std::vector<std::vector<std::vector<double>>> rhoij_all;
    rhoij_all.resize(nat);
    for(int iat = 0; iat < nat; iat ++)
    {
        rhoij_all[iat] = paw_atom_list[iat].get_rhoij();
    }

    return rhoij_all;
}
*/

void Paw_Cell::get_rhoijp(std::vector<std::vector<double>> & rhoijp,
        std::vector<std::vector<int>> & rhoijselect, std::vector<int> & nrhoijsel)
{

    rhoijp.resize(nat);
    rhoijselect.resize(nat);
    nrhoijsel.resize(nat);

    for(int iat = 0; iat < nat; iat ++)
    {
        paw_atom_list[iat].convert_rhoij();
        rhoijp[iat] = paw_atom_list[iat].get_rhoijp();
        rhoijselect[iat] = paw_atom_list[iat].get_rhoijselect();
        nrhoijsel[iat] = paw_atom_list[iat].get_nrhoijsel();
    }
}

void Paw_Cell::paw_nl_psi(const int mode, const std::complex<double> * psi, std::complex<double> * vnlpsi)
{
    ModuleBase::TITLE("Paw_Cell","paw_nl_psi");

    if(mode == 1)
    {
        for(int ipw = 0; ipw < npw; ipw ++)
        {
            vnlpsi[ipw] = psi[ipw];
        }        
    }
    else
    {
        for(int ipw = 0; ipw < npw; ipw ++)
        {
            vnlpsi[ipw] = 0.0;
        }
    }

    for(int iat = 0; iat < nat; iat ++)
    {
        // ca : <ptilde(G)|psi(G)>
        // = \sum_G [\int f(r)r^2j_l(r)dr] * [(-i)^l] * [ylm(\hat{G})] * [exp(-GR_I)] *psi(G)
        // = \sum_ipw ptilde * fact * ylm * sk * psi (in the code below)
        // This is what is called 'becp' in nonlocal pp
        std::vector<std::complex<double>> ca;

        const int it = atom_type[iat];
        const int nproj = paw_element_list[it].get_mstates();
        const int proj_start = start_iprj_ia[iat];

        ca.resize(nproj);

        for(int iproj = 0; iproj < nproj; iproj ++)
        {
            ca[iproj] = 0.0;
            
            // consider use blas subroutine for this part later
            for(int ipw = 0; ipw < npw; ipw ++)
            {
                ca[iproj] += psi[ipw] * std::conj(vkb[iproj+proj_start][ipw]);
            }
        }

#ifdef __MPI
        Parallel_Reduce::reduce_pool(ca.data(), nproj);
#endif

        // sum_ij D_ij ca_j
        std::vector<std::complex<double>> v_ca;
        v_ca.resize(nproj);

        for(int iproj = 0; iproj < nproj; iproj ++)
        {
            v_ca[iproj] = 0.0;
            for(int jproj = 0; jproj < nproj; jproj ++)
            {
                if(mode == 0) // V_{NL}|psi>
                {
                    v_ca[iproj] += paw_atom_list[iat].get_dij()[current_spin][iproj*nproj+jproj] * ca[jproj];
                }
                else if(mode == 1) // (S+I)|psi>
                {
                    v_ca[iproj] += paw_atom_list[iat].get_sij()[iproj*nproj+jproj] * ca[jproj];
                }
            }
        }

        // vnl(ipw) = \sum_i ptilde_{iproj}(G) v_ca[iproj]
        for(int ipw = 0; ipw < npw; ipw ++)
        {
            for(int iproj = 0; iproj < nproj; iproj ++)
            {
                vnlpsi[ipw] += vkb[iproj+proj_start][ipw] * v_ca[iproj];
            }
        }
    }
}

void Paw_Cell::paw_nl_force(const std::complex<double> * psi, const double * epsilon, const double * weight, const int nbands , double * force)
{
    ModuleBase::TITLE("Paw_Cell","paw_nl_force");

    for(int i = 0; i < nat * 3; i ++)
    {
        force[i] = 0.0;
    }

    for(int iband = 0; iband < nbands; iband ++)
    {
        if(weight[iband] < 1e-8) continue;
        for(int iat = 0; iat < nat; iat ++)
        {
            // ca : <ptilde(G)|psi(G)>
            // = \sum_G [\int f(r)r^2j_l(r)dr] * [(-i)^l] * [ylm(\hat{G})] * [exp(-GR_I)] *psi(G)
            // = \sum_ipw ptilde * fact * ylm * sk * psi (in the code below)
            // This is what is called 'becp' in nonlocal pp
            // (but complex conjugate)
            std::vector<std::complex<double>> ca;
            std::vector<std::vector<std::complex<double>>> dca;

            const int it = atom_type[iat];
            const int nproj = paw_element_list[it].get_mstates();
            const int proj_start = start_iprj_ia[iat];

            ca.resize(nproj);
            dca.resize(3);
            for(int i = 0; i < 3; i ++)
            {
                dca[i].resize(nproj);
            }

            for(int iproj = 0; iproj < nproj; iproj ++)
            {
                ca[iproj] = 0.0;
                dca[0][iproj] = 0.0;
                dca[1][iproj] = 0.0;
                dca[2][iproj] = 0.0;
                
                // consider use blas subroutine for this part later
                for(int ipw = 0; ipw < npw; ipw ++)
                {
                    std::complex<double> overlp = psi[iband*npwx+ipw] * std::conj(vkb[iproj+proj_start][ipw]);
                    ca[iproj] += overlp;
                    dca[0][iproj] += overlp * ig[ipw][0];
                    dca[1][iproj] += overlp * ig[ipw][1];
                    dca[2][iproj] += overlp * ig[ipw][2];
                }
            }

#ifdef __MPI
            Parallel_Reduce::reduce_pool(ca.data(), nproj);
            Parallel_Reduce::reduce_pool(dca[0].data(), nproj);
            Parallel_Reduce::reduce_pool(dca[1].data(), nproj);
            Parallel_Reduce::reduce_pool(dca[2].data(), nproj);
#endif
            // sum_ij (D_ij - epsilon_n O_ij) ca_j
            std::vector<std::complex<double>> v_ca;
            v_ca.resize(nproj);

            for(int iproj = 0; iproj < nproj; iproj ++)
            {
                v_ca[iproj] = 0.0;
                for(int jproj = 0; jproj < nproj; jproj ++)
                {
                    double coeff = paw_atom_list[iat].get_dij()[current_spin][iproj*nproj+jproj] -
                        paw_atom_list[iat].get_sij()[iproj*nproj+jproj] * epsilon[iband];
                    v_ca[iproj] += coeff * ca[jproj];
                }
            }

            // force += conjg(v_ca[iproj]) * d_ca[iproj]
            // \sum_i ptilde_{iproj}(G) v_ca[iproj]
            for(int iproj = 0; iproj < nproj; iproj ++)
            {
                force[iat*3] -= (v_ca[iproj] * std::conj(dca[0][iproj])).real() * weight[iband];
                force[iat*3+1] -= (v_ca[iproj] * std::conj(dca[1][iproj])).real() * weight[iband];
                force[iat*3+2] -= (v_ca[iproj] * std::conj(dca[2][iproj])).real() * weight[iband];
            }
        }
    }

    for(int i = 0; i < nat*3; i ++) force[i] = force[i] * 2.0;
}