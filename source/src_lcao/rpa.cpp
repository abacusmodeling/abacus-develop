#include "rpa.h"

#include "../module_base/constants.h"
#include "../module_base/global_variable.h"
#include "../module_base/matrix3.h"
#include "../src_pw/global.h"
#include "global_fp.h"

#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <mpi.h>

#include "../src_ri/conv_coulomb_pot-inl.h"
#include "../src_ri/conv_coulomb_pot.h"
#include "../src_ri/conv_coulomb_pot_k-template.h"
#include "../src_ri/conv_coulomb_pot_k.h"
#include "../src_ri/exx_abfs.h"
#include "../src_ri/exx_abfs-construct_orbs.h"
#include "../src_ri/exx_abfs-parallel-distribute-htime.h"
#include "../src_ri/exx_abfs-parallel-distribute-kmeans.h"
#include "../src_ri/exx_abfs-parallel-distribute-order.h"
#include "../src_ri/exx_abfs-parallel.h"
#include "../src_ri/exx_abfs-util.h"
#include "../src_ri/exx_abfs-io.h"
#include "../src_ri/exx_abfs-abfs_index.h"

namespace GlobalC
{
ModuleRPA::DFT_RPA_interface rpa;
//Exx_Lcao rpa_exx_lcao(GlobalC::exx_global.info);
}

namespace ModuleRPA
{
void DFT_RPA_interface::out_for_RPA(Local_Orbital_wfc &lowf,Local_Orbital_Charge &loc)
{
    ModuleBase::TITLE("DFT_RPA_interface", "out_for_RPA");
    this->out_bands();
    this->out_eigen_vector(lowf);
    this->out_struc();
    
    
    //void exx_inti();
    //GlobalC::rpa_exx_lcao=GlobalC::exx_lcao;
    if(GlobalV::DFT_FUNCTIONAL=="default")
    {
        this->exx_init();
        cout<<"rpa_pca_threshold: "<<GlobalC::rpa_exx_lcao.info.pca_threshold<<endl;
        cout<<"rpa_ccp_rmesh_times: "<<GlobalC::rpa_exx_lcao.info.ccp_rmesh_times<<endl;
        GlobalC::rpa_exx_lcao.cal_exx_ions(*lowf.ParaV);
        GlobalC::rpa_exx_lcao.cal_exx_elec(loc, lowf.wfc_k_grid);
    }
    cout<<"rpa_lcao_exx(Ha): "<< std::fixed << std::setprecision(15)<<GlobalC::rpa_exx_lcao.get_energy()/2.0<<endl;
    this->out_Cs();
    this->out_coulomb_k();
    //cout << "EXX_energy(Ha):  " << std::setprecision(6) << GlobalC::exx_lcao.get_energy()/2.0 << endl;
    cout<<"etxc(Ha):"<< std::fixed << std::setprecision(15)<<GlobalC::en.etxc/2.0<<endl;
    cout<<"etot(Ha):"<< std::fixed << std::setprecision(15)<<GlobalC::en.etot/2.0<<endl;
    cout<<"Etot_without_rpa(Ha):"<< std::fixed << std::setprecision(15)<<
        (GlobalC::en.etot-GlobalC::en.etxc+GlobalC::rpa_exx_lcao.get_energy())/2.0<<endl;
    //cout<<"Etot(Ha):"<<GlobalC::en.etot/2.0<<endl;
    //cout<<"etxcc(Ha):"<<GlobalC::en.etxcc<<endl;
    return;
}
void DFT_RPA_interface::out_eigen_vector(Local_Orbital_wfc &lowf)
{

    ModuleBase::TITLE("DFT_RPA_interface", "out_eigen_vector");

    const int nks_tot = GlobalV::NSPIN == 2 ? (int)GlobalC::kv.nks / 2 : GlobalC::kv.nks;
    const int npsin_tmp = GlobalV::NSPIN == 2 ? 2 : 1;
    const std::complex<double> zero(0.0, 0.0);
    // for (int ik = 0; ik != nks_tot; ik++)
    //     this->print_complex_matrix("wfc_k", wfc[ik]);
    // for (int ik = 0; ik != nks_tot; ik++)
    //     this->print_complex_matrix("dm_k", GlobalC::LOC.wfc_dm_2d.dm_k[ik]);

    // GlobalV::ofs_running << "GlobalV::NLOCAL:" << GlobalV::NLOCAL
    // << "  GlobalC::this->ParaV.nrow:" << GlobalC::this->ParaV.nrow
    // << "  nb:" << wfc[0].nr
    // << "  GlobalC::this->ParaV.ncol:" << GlobalC::this->ParaV.ncol
    // << "  iw:" << wfc[0].nc << std::endl;

    for (int ik = 0; ik < nks_tot; ik++)
    {
        std::stringstream ss;
        ss << "KS_eigenvector_" << ik << ".dat";

        std::ofstream ofs;
        if (GlobalV::MY_RANK == 0)
            ofs.open(ss.str().c_str(), std::ios::out);

        for (int is = 0; is < npsin_tmp; is++)
        {
            // ofs<<GlobalC::kv.isk[ik+nks_tot*is]<<endl;
            ofs << ik + 1 << endl;
            for (int ib_global = 0; ib_global < GlobalV::NBANDS; ++ib_global)
            {
                std::vector<std::complex<double>> wfc_iks(GlobalV::NLOCAL, zero);

                const int ib_local = lowf.ParaV->trace_loc_col[ib_global];

                if (ib_local >= 0)
                    for (int ir = 0; ir < lowf.wfc_k[ik + nks_tot * is].nc; ir++)
                        wfc_iks[lowf.ParaV->MatrixInfo.row_set[ir]] = lowf.wfc_k[ik + nks_tot * is](ib_local, ir);

                std::vector<std::complex<double>> tmp = wfc_iks;
#ifdef __MPI
                MPI_Allreduce(&tmp[0], &wfc_iks[0], GlobalV::NLOCAL, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
#endif
                for (int iw = 0; iw < GlobalV::NLOCAL; iw++)
                    ofs << std::setw(21) << std::fixed << std::setprecision(15) << wfc_iks[iw].real() << std::setw(21)
                        << std::fixed << std::setprecision(15) << wfc_iks[iw].imag() << std::endl;
            } // ib
        } // is
        ofs.close();
    } // ik
    return;
}

void DFT_RPA_interface::out_struc()
{
    if (GlobalV::MY_RANK != 0)
        return;
    ModuleBase::TITLE("DFT_RPA_interface", "out_struc");
    double TWOPI_Bohr2A = ModuleBase::TWO_PI * ModuleBase::BOHR_TO_A;
    const int nks_tot = GlobalV::NSPIN == 2 ? (int)GlobalC::kv.nks / 2 : GlobalC::kv.nks;
    ModuleBase::Matrix3 lat = GlobalC::ucell.latvec / ModuleBase::BOHR_TO_A;
    ModuleBase::Matrix3 G = GlobalC::ucell.G * TWOPI_Bohr2A;
    std::stringstream ss;
    ss << "stru_out";
    std::ofstream ofs;
    ofs.open(ss.str().c_str(), std::ios::out);
    ofs << lat.e11 << std::setw(15) << lat.e12 << std::setw(15) << lat.e13 << std::endl;
    ofs << lat.e21 << std::setw(15) << lat.e22 << std::setw(15) << lat.e23 << std::endl;
    ofs << lat.e31 << std::setw(15) << lat.e32 << std::setw(15) << lat.e33 << std::endl;

    ofs << G.e11 << std::setw(15) << G.e12 << std::setw(15) << G.e13 << std::endl;
    ofs << G.e21 << std::setw(15) << G.e22 << std::setw(15) << G.e23 << std::endl;
    ofs << G.e31 << std::setw(15) << G.e32 << std::setw(15) << G.e33 << std::endl;

    ofs << GlobalC::kv.nmp[0] << std::setw(6) << GlobalC::kv.nmp[1] << std::setw(6) << GlobalC::kv.nmp[2]
        << std::setw(6) << std::endl;

    for (int ik = 0; ik != nks_tot; ik++)
        ofs << std::setw(15) << std::fixed << std::setprecision(9) << GlobalC::kv.kvec_c[ik].x * TWOPI_Bohr2A
            << std::setw(15) << std::fixed << std::setprecision(9) << GlobalC::kv.kvec_c[ik].y * TWOPI_Bohr2A
            << std::setw(15) << std::fixed << std::setprecision(9) << GlobalC::kv.kvec_c[ik].z * TWOPI_Bohr2A
            << std::endl;
    ofs.close();
    return;
}

void DFT_RPA_interface::out_bands()
{
    ModuleBase::TITLE("DFT_RPA_interface", "out_bands");
    if (GlobalV::MY_RANK != 0)
        return;
    const int nks_tot = GlobalV::NSPIN == 2 ? (int)GlobalC::kv.nks / 2 : GlobalC::kv.nks;
    const int nspin_tmp = GlobalV::NSPIN == 2 ? 2 : 1;
    std::stringstream ss;
    ss << "band_out";
    std::ofstream ofs;
    ofs.open(ss.str().c_str(), std::ios::out);
    ofs << GlobalC::kv.nks << endl;
    ofs << GlobalV::NSPIN << endl;
    ofs << GlobalV::NBANDS << endl;
    ofs << GlobalV::NLOCAL << endl;
    ofs << (GlobalC::en.ef / 2.0) << endl;

    for (int ik = 0; ik != nks_tot; ik++)
    {
        for (int is = 0; is != nspin_tmp; is++)
        {
            ofs << std::setw(6) << ik + 1 << std::setw(6) << is + 1 << endl;
            for (int ib = 0; ib != GlobalV::NBANDS; ib++)
                ofs << std::setw(5) << ib + 1 << "   " << std::setw(8) << GlobalC::wf.wg(ik, ib) * nks_tot
                    << std::setw(18) << std::fixed << std::setprecision(8)
                    << GlobalC::wf.ekb[ik + is * nks_tot][ib] / 2.0 << std::setw(18) << std::fixed
                    << std::setprecision(8) << GlobalC::wf.ekb[ik + is * nks_tot][ib] * ModuleBase::Ry_to_eV << endl;
        }
    }
    ofs.close();
    return;
}

void DFT_RPA_interface::out_Cs()
{
    auto cal_atom_centres_core = [](const std::vector<std::pair<size_t, size_t>> &atom_pairs_core) -> std::set<size_t> {
        std::set<size_t> atom_centres_core;
        for (const std::pair<size_t, size_t> &atom_pair: atom_pairs_core)
        {
            atom_centres_core.insert(atom_pair.first);
            atom_centres_core.insert(atom_pair.second);
        }
        return atom_centres_core;
    };

    const std::set<size_t> atom_centres_core = cal_atom_centres_core(GlobalC::rpa_exx_lcao.atom_pairs_core);
    // map<size_t, map<size_t, map<Abfs::Vector3_Order<int>, shared_ptr<ModuleBase::matrix>>>> Cs_temp;
    // map<size_t, map<size_t, map<Abfs::Vector3_Order<int>, shared_ptr<ModuleBase::matrix>>>> Cs_m;
    // Cs_temp = Abfs::cal_Cs(atom_centres_core,
    //                        GlobalC::rpa_exx_lcao.m_abfs_abfs,
    //                        GlobalC::rpa_exx_lcao.m_abfslcaos_lcaos,
    //                        GlobalC::rpa_exx_lcao.index_abfs,
    //                        GlobalC::rpa_exx_lcao.index_lcaos,
    //                        GlobalC::rpa_exx_lcao.info.c_threshold,
    //                        GlobalC::rpa_exx_lcao.Cws,
    //                        GlobalC::rpa_exx_lcao.Vws);
    // Cs_m = Abfs::cal_mps(GlobalC::rpa_exx_lcao.Born_von_Karman_period, Cs_temp);

    std::stringstream ss;
    ss << "Cs_data.txt";
    std::ofstream ofs;
    ofs.open(ss.str().c_str(), std::ios::out);
    ofs << atom_centres_core.size() << "    " << 0 << endl;
    for (auto &Ip: GlobalC::rpa_exx_lcao.Cps)
    {
        size_t I = Ip.first;
        size_t i_num = GlobalC::ucell.atoms[GlobalC::ucell.iat2it[I]].nw;
        size_t mu_num = GlobalC::rpa_exx_lcao.index_abfs[GlobalC::ucell.iat2it[I]].count_size;
        for (auto &Jp: Ip.second)
        {
            size_t J = Jp.first;
            size_t j_num = GlobalC::ucell.atoms[GlobalC::ucell.iat2it[J]].nw;
            for (auto &Rp: Jp.second)
            {
                auto R = Rp.first;
                ofs << I + 1 << "   " << J + 1 << "   " << R.x << "   " << R.y << "   " << R.z << "   " << i_num
                    << endl;
                ofs << j_num << "   " << mu_num << endl;
                auto &tmp_Cs = *Rp.second;
                size_t ele_num = tmp_Cs.nr * tmp_Cs.nc;
                for (int it = 0; it != ele_num; it++)
                    ofs << std::setw(15) << std::fixed << std::setprecision(9) << tmp_Cs.c[it] << endl;
            }
        }
    }
    ofs.close();
    return;
}

void DFT_RPA_interface::out_coulomb_k()
{
    // auto cal_atom_centres_core = [](const std::vector<std::pair<size_t, size_t>> &atom_pairs_core) -> std::set<size_t> {
    //     std::set<size_t> atom_centres_core;
    //     for (const std::pair<size_t, size_t> &atom_pair: atom_pairs_core)
    //     {
    //         atom_centres_core.insert(atom_pair.first);
    //         atom_centres_core.insert(atom_pair.second);
    //     }
    //     return atom_centres_core;
    // };
    // map<size_t,
    //     map<size_t,
    //         map<Abfs::Vector3_Order<int>,
    //             shared_ptr<ModuleBase::matrix>>>>
    //     Vs_m; // Vs[iat1][iat2][box2]
    // map<size_t, map<size_t, map<Abfs::Vector3_Order<int>, shared_ptr<ModuleBase::matrix>>>> Vps_m;
    // Vs_m = Abfs::cal_Vs(GlobalC::rpa_exx_lcao.atom_pairs_core_origin,
    //                     GlobalC::rpa_exx_lcao.m_abfs_abfs,
    //                     GlobalC::rpa_exx_lcao.index_abfs,
    //                     GlobalC::rpa_exx_lcao.info.ccp_rmesh_times,
    //                     GlobalC::rpa_exx_lcao.info.v_threshold,
    //                     GlobalC::rpa_exx_lcao.Vws);
    // Vps_m = Abfs::cal_mps(GlobalC::rpa_exx_lcao.Born_von_Karman_period, Vs_m);

    int all_mu = 0;
    vector<int> mu_shift(GlobalC::rpa_exx_lcao.Vps.size());
    for (auto &Ip: GlobalC::rpa_exx_lcao.Vps)
    {
        auto I = Ip.first;
        mu_shift[I] = all_mu;
        all_mu += GlobalC::rpa_exx_lcao.index_abfs[GlobalC::ucell.iat2it[I]].count_size;
    }
    std::stringstream ss;
    ss << "coulomb_mat.txt";
    std::ofstream ofs;
    ofs.open(ss.str().c_str(), std::ios::out);

    ofs << GlobalC::kv.nks << endl;
    for (auto &Ip: GlobalC::rpa_exx_lcao.Vps)
    {
        auto I = Ip.first;
        size_t mu_num = GlobalC::rpa_exx_lcao.index_abfs[GlobalC::ucell.iat2it[I]].count_size;
        for (auto &Jp: Ip.second)
        {
            auto J = Jp.first;
            size_t nu_num = GlobalC::rpa_exx_lcao.index_abfs[GlobalC::ucell.iat2it[J]].count_size;

            for (int ik = 0; ik != GlobalC::kv.nks; ik++)
            {
                ModuleBase::ComplexMatrix tmp_Vk(mu_num, nu_num);
                tmp_Vk.zero_out();
                ofs << all_mu << "   " << mu_shift[I] + 1 << "   " << mu_shift[I] + mu_num << "  " << mu_shift[J] + 1
                    << "   " << mu_shift[J] + nu_num << endl;
                ofs << ik + 1 << "  " << GlobalC::kv.wk[ik] / 2.0 << endl;
                for (auto &Rp: Jp.second)
                {
                    auto R = Rp.first;
                    ModuleBase::ComplexMatrix tmp_VR(*Rp.second);
                    const double arg
                        = 1 * (GlobalC::kv.kvec_c[ik] * (R * GlobalC::ucell.latvec)) * ModuleBase::TWO_PI; // latvec
                    const complex<double> kphase = complex<double>(cos(arg), sin(arg));
                    tmp_Vk += tmp_VR * kphase;
                }
                for (int i = 0; i != tmp_Vk.size; i++)
                {
                    ofs << std::setw(21) << std::fixed << std::setprecision(12) << tmp_Vk.c[i].real() << std::setw(21)
                        << std::fixed << std::setprecision(12) << tmp_Vk.c[i].imag() << endl;
                }
            }
        }
    }
    ofs.close();
}

void DFT_RPA_interface::exx_init()
{
    cout<<"rpa_exx_init!!!"<<endl;
    #ifdef __MPI
    if (GlobalC::exx_global.info.separate_loop)
    {
        GlobalC::rpa_exx_lcao.Hexx_para.mixing_mode = Exx_Abfs::Parallel::Communicate::Hexx::Mixing_Mode::No;
        GlobalC::rpa_exx_lcao.Hexx_para.mixing_beta = 0;
    }
    else
    {
        if ("plain" == GlobalC::CHR.mixing_mode)
        {
            GlobalC::rpa_exx_lcao.Hexx_para.mixing_mode = Exx_Abfs::Parallel::Communicate::Hexx::Mixing_Mode::Plain;
        }
        else if ("pulay" == GlobalC::CHR.mixing_mode)
        {
            GlobalC::rpa_exx_lcao.Hexx_para.mixing_mode = Exx_Abfs::Parallel::Communicate::Hexx::Mixing_Mode::Pulay;
        }
        else
        {
            throw std::invalid_argument("exx mixing error. exx_separate_loop==false, mixing_mode!=plain or pulay");
        }
        GlobalC::rpa_exx_lcao.Hexx_para.mixing_beta = GlobalC::CHR.mixing_beta;
    }
#endif

    GlobalC::rpa_exx_lcao.lcaos = Exx_Abfs::Construct_Orbs::change_orbs(GlobalC::ORB, GlobalC::rpa_exx_lcao.kmesh_times);
#ifdef __MPI
    Exx_Abfs::Util::bcast(GlobalC::rpa_exx_lcao.info.files_abfs, 0, MPI_COMM_WORLD);
#endif

    const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> abfs_same_atom
        = Exx_Abfs::Construct_Orbs::abfs_same_atom(GlobalC::rpa_exx_lcao.lcaos,
                                                   GlobalC::rpa_exx_lcao.kmesh_times,
                                                   GlobalC::rpa_exx_lcao.info.pca_threshold); // Peize Lin test
    if (GlobalC::rpa_exx_lcao.info.files_abfs.empty())
    {
        GlobalC::rpa_exx_lcao.abfs = abfs_same_atom;
    }
    else
    {
        GlobalC::rpa_exx_lcao.abfs = Exx_Abfs::IO::construct_abfs(abfs_same_atom,
                                                          GlobalC::ORB,
                                                          GlobalC::rpa_exx_lcao.info.files_abfs,
                                                          GlobalC::rpa_exx_lcao.kmesh_times);
    }

    switch (GlobalC::rpa_exx_lcao.info.hybrid_type)
    {
    case Exx_Global::Hybrid_Type::HF:
        GlobalC::rpa_exx_lcao.abfs_ccp
            = Conv_Coulomb_Pot_K::cal_orbs_ccp(GlobalC::rpa_exx_lcao.abfs, Conv_Coulomb_Pot_K::Ccp_Type::Hf, {}, GlobalC::rpa_exx_lcao.info.ccp_rmesh_times);
        break;
    case Exx_Global::Hybrid_Type::No:
        GlobalC::rpa_exx_lcao.abfs_ccp
            = Conv_Coulomb_Pot_K::cal_orbs_ccp(GlobalC::rpa_exx_lcao.abfs, Conv_Coulomb_Pot_K::Ccp_Type::Hf, {}, GlobalC::rpa_exx_lcao.info.ccp_rmesh_times);
        break;
    case Exx_Global::Hybrid_Type::PBE0:
        GlobalC::rpa_exx_lcao.abfs_ccp
            = Conv_Coulomb_Pot_K::cal_orbs_ccp(GlobalC::rpa_exx_lcao.abfs, Conv_Coulomb_Pot_K::Ccp_Type::Ccp, {}, GlobalC::rpa_exx_lcao.info.ccp_rmesh_times);
        break;
    case Exx_Global::Hybrid_Type::HSE:
        GlobalC::rpa_exx_lcao.abfs_ccp = Conv_Coulomb_Pot_K::cal_orbs_ccp(GlobalC::rpa_exx_lcao.abfs,
                                                    Conv_Coulomb_Pot_K::Ccp_Type::Hse,
                                                    {{"hse_omega", GlobalC::rpa_exx_lcao.info.hse_omega}},
                                                    GlobalC::rpa_exx_lcao.info.ccp_rmesh_times);
        break;
    default:
        throw std::domain_error(ModuleBase::GlobalFunc::TO_STRING(__FILE__) + " line "
                                + ModuleBase::GlobalFunc::TO_STRING(__LINE__));
        break;
    }

    for (size_t T = 0; T != GlobalC::rpa_exx_lcao.abfs.size(); ++T)
    {
        Exx_Abfs::Lmax = std::max(Exx_Abfs::Lmax, static_cast<int>(GlobalC::rpa_exx_lcao.abfs[T].size()) - 1);
    }


    const ModuleBase::Element_Basis_Index::Range &&range_lcaos = Exx_Abfs::Abfs_Index::construct_range(GlobalC::rpa_exx_lcao.lcaos);
    GlobalC::rpa_exx_lcao.index_lcaos = ModuleBase::Element_Basis_Index::construct_index(range_lcaos);

    const ModuleBase::Element_Basis_Index::Range &&range_abfs = Exx_Abfs::Abfs_Index::construct_range(GlobalC::rpa_exx_lcao.abfs);
    GlobalC::rpa_exx_lcao.index_abfs = ModuleBase::Element_Basis_Index::construct_index(range_abfs);

    GlobalC::rpa_exx_lcao.m_abfs_abfs.init(2, GlobalC::rpa_exx_lcao.kmesh_times, (1 + GlobalC::rpa_exx_lcao.info.ccp_rmesh_times) / 2.0);
   
    GlobalC::rpa_exx_lcao.m_abfs_abfs.init_radial(GlobalC::rpa_exx_lcao.abfs_ccp, GlobalC::rpa_exx_lcao.abfs);
  
  
    GlobalC::rpa_exx_lcao.m_abfslcaos_lcaos.init(1, GlobalC::rpa_exx_lcao.kmesh_times, 1);
  
    GlobalC::rpa_exx_lcao.m_abfslcaos_lcaos.init_radial(GlobalC::rpa_exx_lcao.abfs_ccp, GlobalC::rpa_exx_lcao.lcaos, GlobalC::rpa_exx_lcao.lcaos);
  
    GlobalC::rpa_exx_lcao.Born_von_Karman_period = ModuleBase::Vector3<int>{GlobalC::kv.nmp[0], GlobalC::kv.nmp[1], GlobalC::kv.nmp[2]};
  
}

void DFT_RPA_interface::exx_cal_ions()
{
    cout<<"rpa_exx_cal_ions!!!"<<endl;
    auto cal_atom_centres_core = [](const std::vector<std::pair<size_t, size_t>> &atom_pairs_core) -> std::set<size_t> {
        std::set<size_t> atom_centres_core;
        for (const std::pair<size_t, size_t> &atom_pair: atom_pairs_core)
        {
            atom_centres_core.insert(atom_pair.first);
            atom_centres_core.insert(atom_pair.second);
        }
        return atom_centres_core;
    };

#ifdef __MPI
    if (GlobalC::rpa_exx_lcao.atom_pairs_core_origin.empty())
        switch (GlobalC::rpa_exx_lcao.info.distribute_type)
        {
        case Exx_Lcao::Distribute_Type::Htime:
            GlobalC::rpa_exx_lcao.atom_pairs_core_origin
                = Exx_Abfs::Parallel::Distribute::Htime::distribute(GlobalC::rpa_exx_lcao.Born_von_Karman_period, GlobalC::rpa_exx_lcao.info.ccp_rmesh_times);
            break;
        case Exx_Lcao::Distribute_Type::Kmeans2:
            GlobalC::rpa_exx_lcao.atom_pairs_core_origin = Exx_Abfs::Parallel::Distribute::Kmeans::distribute_kmeans2(MPI_COMM_WORLD);
            break;
        case Exx_Lcao::Distribute_Type::Kmeans1:
            GlobalC::rpa_exx_lcao.atom_pairs_core_origin
                = Exx_Abfs::Parallel::Distribute::Kmeans::distribute_kmeans1(MPI_COMM_WORLD, GlobalC::rpa_exx_lcao.info.ccp_rmesh_times);
            break;
        case Exx_Lcao::Distribute_Type::Order:
            GlobalC::rpa_exx_lcao.atom_pairs_core_origin = Exx_Abfs::Parallel::Distribute::Order::distribute(GlobalC::rpa_exx_lcao.info.ccp_rmesh_times);
            break;
        default:
            throw std::domain_error(ModuleBase::GlobalFunc::TO_STRING(__FILE__) + " line "
                                    + ModuleBase::GlobalFunc::TO_STRING(__LINE__));
            break;
         
        }
#endif
   
    GlobalC::rpa_exx_lcao.init_radial_table_ions(cal_atom_centres_core(GlobalC::rpa_exx_lcao.atom_pairs_core_origin), GlobalC::rpa_exx_lcao.atom_pairs_core_origin);
 
    GlobalC::rpa_exx_lcao.Vs = Abfs::cal_Vs(GlobalC::rpa_exx_lcao.atom_pairs_core_origin, GlobalC::rpa_exx_lcao.m_abfs_abfs, GlobalC::rpa_exx_lcao.index_abfs, GlobalC::rpa_exx_lcao.info.ccp_rmesh_times, GlobalC::rpa_exx_lcao.info.v_threshold, GlobalC::rpa_exx_lcao.Vws);

    Abfs::delete_empty_ptrs(GlobalC::rpa_exx_lcao.Vws);

    GlobalC::rpa_exx_lcao.Vps = Abfs::cal_mps(GlobalC::rpa_exx_lcao.Born_von_Karman_period, GlobalC::rpa_exx_lcao.Vs);
   
    GlobalC::rpa_exx_lcao.atom_pairs_core = Abfs::get_atom_pair(GlobalC::rpa_exx_lcao.Vps);

    const std::set<size_t> atom_centres_core = cal_atom_centres_core(GlobalC::rpa_exx_lcao.atom_pairs_core);

    GlobalC::rpa_exx_lcao.H_atom_pairs_core = Abfs::get_H_pairs_core(GlobalC::rpa_exx_lcao.atom_pairs_core);
   
    // GlobalC::rpa_exx_lcao.Cs = Abfs::cal_Cs(atom_centres_core,
    //                   GlobalC::rpa_exx_lcao.m_abfs_abfs,
    //                   GlobalC::rpa_exx_lcao.m_abfslcaos_lcaos,
    //                   GlobalC::rpa_exx_lcao.index_abfs,
    //                   GlobalC::rpa_exx_lcao.index_lcaos,
    //                   GlobalC::rpa_exx_lcao.info.c_threshold,
    //                   GlobalC::rpa_exx_lcao.Cws,
    //                   GlobalC::rpa_exx_lcao.Vws);

    Abfs::delete_empty_ptrs(GlobalC::rpa_exx_lcao.Cws); 
 //   GlobalC::rpa_exx_lcao.Cps = Abfs::cal_mps(GlobalC::rpa_exx_lcao.Born_von_Karman_period, GlobalC::rpa_exx_lcao.Cs);
 
}


void DFT_RPA_interface::print_matrix(char *desc, const ModuleBase::matrix &mat)
{
    int nr = mat.nr;
    int nc = mat.nc;
    printf("\n %s\n", desc);
    for (int i = 0; i < nr; i++)
    {
        for (int j = 0; j < nc; j++)
            printf(" %9.5f", mat.c[i * nc + j]);
        printf("\n");
    }
}
void DFT_RPA_interface::print_complex_matrix(char *desc, const ModuleBase::ComplexMatrix &mat)
{
    int nr = mat.nr;
    int nc = mat.nc;
    printf("\n %s\n", desc);
    for (int i = 0; i < nr; i++)
    {
        for (int j = 0; j < nc; j++)
            printf("%10.6f,%8.6f ", mat.c[i * nc + j].real(), mat.c[i * nc + j].imag());
        printf("\n");
    }
}

} // namespace ModuleRPA