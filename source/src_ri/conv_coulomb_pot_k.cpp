#include "conv_coulomb_pot_k.h"
#include "../module_base/constants.h"
#include "../module_orbital/ORB_atomic_lm.h"
#include "../src_pw/global.h"

std::vector<double>
Conv_Coulomb_Pot_K::cal_psi_ccp(const std::vector<double> &psif) {
    std::vector<double> psik2_ccp(psif.size());
    for (size_t ik = 0; ik < psif.size(); ++ik)
        psik2_ccp[ik] = ModuleBase::FOUR_PI * psif[ik];
    return psik2_ccp;
}

std::vector<double>
Conv_Coulomb_Pot_K::cal_psi_hf(const std::vector<double> &psif,
                               const std::vector<double> &k_radial,
                               const double omega = 0) {
    // cout<<"cal_psi_hf"<<endl;

    // double hf_rcut = 0.867 * pow((kv.nks * ucell.omega) / NSPIN, 0.3334);
    // double dr = 0.01;
    // int nr = 100 * hf_rcut / dr;
    // std::vector<double> psik2_ccp(psif.size());

    // for (size_t ik = 0; ik < psif.size(); ++ik) {
    //     double sum_r = 0;
    //     for (size_t ir = 0; ir != nr; ++ir) {
    //         sum_r += 0.5 * std::erfc(omega * (dr * ir - hf_rcut)) *
    //                  std::sin(dr * ir * k_radial[ik]);
    //     }
    //     psik2_ccp[ik] = FOUR_PI * psif[ik] * dr * sum_r * k_radial[ik];
    // }

    // Sphere truction -- Spencer
    double Rc = pow(0.75 * GlobalC::ucell.omega / (ModuleBase::PI), 0.3333334);
    cout << "hf_Rc:  " << Rc << endl;
    std::vector<double> psik2_ccp(psif.size());
    for (size_t ik = 0; ik < psif.size(); ++ik)
        psik2_ccp[ik] =
            ModuleBase::FOUR_PI * psif[ik] * (1 - std::cos(k_radial[ik] * Rc));
    return psik2_ccp;
}

std::vector<double>
Conv_Coulomb_Pot_K::cal_psi_hse(const std::vector<double> &psif,
                                const std::vector<double> &k_radial,
                                const double omega) {
    std::vector<double> psik2_ccp(psif.size());
    for (size_t ik = 0; ik < psif.size(); ++ik)
        psik2_ccp[ik] = ModuleBase::FOUR_PI * psif[ik] *
                        (1 - std::exp(-(k_radial[ik] * k_radial[ik]) /
                                      (4 * omega * omega)));
    return psik2_ccp;
}

template <>
Numerical_Orbital_Lm Conv_Coulomb_Pot_K::cal_orbs_ccp<Numerical_Orbital_Lm>(
    const Numerical_Orbital_Lm &orbs, const Ccp_Type &ccp_type,
    const std::map<std::string, double> &parameter, const double rmesh_times) {
    std::vector<double> psik2_ccp;
    switch (ccp_type) {
    case Ccp_Type::Ccp:
        psik2_ccp = cal_psi_ccp(orbs.get_psif());
        break;
    case Ccp_Type::Hf:
        psik2_ccp = cal_psi_hf(orbs.get_psif(), orbs.get_k_radial());
        break;
    case Ccp_Type::Hse:
        psik2_ccp = cal_psi_hse(orbs.get_psif(), orbs.get_k_radial(),
                                parameter.at("hse_omega"));
        break;
    default:
        throw(ModuleBase::GlobalFunc::TO_STRING(__FILE__) + " line " +
              ModuleBase::GlobalFunc::TO_STRING(__LINE__));
        break;
    }

    const double dr = orbs.get_rab().back();
    const int Nr = (static_cast<int>(orbs.getNr() * rmesh_times)) | 1;
    std::vector<double> rab(Nr);
    for (size_t ir = 0; ir < std::min(orbs.getNr(), Nr); ++ir)
        rab[ir] = orbs.getRab(ir);
    for (size_t ir = orbs.getNr(); ir < Nr; ++ir)
        rab[ir] = dr;
    std::vector<double> r_radial(Nr);
    for (size_t ir = 0; ir < std::min(orbs.getNr(), Nr); ++ir)
        r_radial[ir] = orbs.getRadial(ir);
    for (size_t ir = orbs.getNr(); ir < Nr; ++ir)
        r_radial[ir] = orbs.get_r_radial().back() + (ir - orbs.getNr()) * dr;

    Numerical_Orbital_Lm orbs_ccp;
    orbs_ccp.set_orbital_info(
        orbs.getLabel(), orbs.getType(), orbs.getL(), orbs.getChi(), Nr,
        ModuleBase::GlobalFunc::VECTOR_TO_PTR(rab),
        ModuleBase::GlobalFunc::VECTOR_TO_PTR(r_radial),
        Numerical_Orbital_Lm::Psi_Type::Psik2,
        ModuleBase::GlobalFunc::VECTOR_TO_PTR(psik2_ccp), orbs.getNk(),
        orbs.getDk(), orbs.getDruniform(), false, true, GlobalV::FORCE);
    return orbs_ccp;
}

template <>
double
Conv_Coulomb_Pot_K::get_rmesh_proportion(const Numerical_Orbital_Lm &orbs,
                                         const double psi_threshold) {
    for (int ir = orbs.getNr() - 1; ir >= 0; --ir) {
        if (std::abs(orbs.getPsi(ir)) >= psi_threshold)
            return static_cast<double>(ir) / orbs.getNr();
    }
    return 0.0;
}
