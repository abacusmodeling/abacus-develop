#include "../input.h"
#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"

void Input::Print(const std::string &fn) const
{
    if (GlobalV::MY_RANK != 0)
        return;

    ModuleBase::TITLE("Input", "Print");

    std::ofstream ofs(fn.c_str());

    //----------------------------------
    // output the information in INPUT.
    //----------------------------------
    ofs << "INPUT_PARAMETERS" << std::endl;
    ofs << std::setiosflags(ios::left);

    ofs << "#Parameters (1.General)" << std::endl;
    ModuleBase::GlobalFunc::OUTP(ofs, "suffix", suffix, "the name of main output directory");
    ModuleBase::GlobalFunc::OUTP(ofs, "latname", latname, "the name of lattice name");
    ModuleBase::GlobalFunc::OUTP(ofs,
                                 "stru_file",
                                 GlobalV::stru_file,
                                 "the filename of file containing atom positions");
    ModuleBase::GlobalFunc::OUTP(ofs,
                                 "kpoint_file",
                                 GlobalV::global_kpoint_card,
                                 "the name of file containing k points");
    ModuleBase::GlobalFunc::OUTP(ofs,
                                 "pseudo_dir",
                                 GlobalV::global_pseudo_dir,
                                 "the directory containing pseudo files");
    ModuleBase::GlobalFunc::OUTP(ofs,
                                 "orbital_dir",
                                 GlobalV::global_orbital_dir,
                                 "the directory containing orbital files");
    // ModuleBase::GlobalFunc::OUTP(ofs, "pseudo_type", GlobalV::global_pseudo_type, "the type pseudo files");
    ModuleBase::GlobalFunc::OUTP(ofs, "pseudo_rcut", pseudo_rcut, "cut-off radius for radial integration");
    ModuleBase::GlobalFunc::OUTP(ofs,
                                 "pseudo_mesh",
                                 pseudo_mesh,
                                 "0: use our own mesh to do radial renormalization; 1: use mesh as in QE");
    ModuleBase::GlobalFunc::OUTP(ofs, "lmaxmax", lmaxmax, "maximum of l channels used");
    ModuleBase::GlobalFunc::OUTP(ofs, "dft_functional", dft_functional, "exchange correlation functional");
    ModuleBase::GlobalFunc::OUTP(ofs, "xc_temperature", xc_temperature, "temperature for finite temperature functionals");
    ModuleBase::GlobalFunc::OUTP(ofs, "calculation", calculation, "test; scf; relax; nscf; ienvelope; istate; sto-scf; sto-md; ofdft; of-md");
    ModuleBase::GlobalFunc::OUTP(ofs, "ntype", ntype, "atom species number");
    ModuleBase::GlobalFunc::OUTP(ofs, "nspin", nspin, "1: single spin; 2: up and down spin; 4: noncollinear spin");
    ModuleBase::GlobalFunc::OUTP(ofs, "kspacing", kspacing, "unit in 1/bohr, should be > 0, default is 0 which means read KPT file");
    ModuleBase::GlobalFunc::OUTP(ofs, "min_dist_coef", min_dist_coef, "factor related to the allowed minimum distance between two atoms");
    ModuleBase::GlobalFunc::OUTP(ofs, "nbands", nbands, "number of bands");
    ModuleBase::GlobalFunc::OUTP(ofs, "nbands_sto", nbands_sto, "number of stochastic bands");
    ModuleBase::GlobalFunc::OUTP(ofs,
                                 "nbands_istate",
                                 nbands_istate,
                                 "number of bands around Fermi level for istate calulation");
    ModuleBase::GlobalFunc::OUTP(ofs, "symmetry", symmetry, "the control of symmetry");
    ModuleBase::GlobalFunc::OUTP(ofs, "init_vel", init_vel, "read velocity from STRU or not");
    ModuleBase::GlobalFunc::OUTP(ofs,
                                 "symmetry_prec",
                                 symmetry_prec,
                                 "accuracy for symmetry"); // LiuXh add 2021-08-12, accuracy for symmetry
    ModuleBase::GlobalFunc::OUTP(ofs, "nelec", nelec, "input number of electrons");
    ModuleBase::GlobalFunc::OUTP(ofs, "tot_magnetization", tot_magnetization, "total magnetization of the system");
    ModuleBase::GlobalFunc::OUTP(ofs, "out_mul", GlobalV::out_mul, " mulliken  charge or not"); // qifeng add 2019/9/10
    ModuleBase::GlobalFunc::OUTP(ofs, "noncolin", noncolin, "using non-collinear-spin");
    ModuleBase::GlobalFunc::OUTP(ofs, "lspinorb", lspinorb, "consider the spin-orbit interaction");
    ModuleBase::GlobalFunc::OUTP(ofs, "kpar", kpar, "devide all processors into kpar groups and k points will be distributed among each group");
    ModuleBase::GlobalFunc::OUTP(ofs, "bndpar", bndpar, "devide all processors into bndpar groups and bands will be distributed among each group");
    ModuleBase::GlobalFunc::OUTP(ofs, "out_freq_elec", out_freq_elec, "the frequency ( >= 0) of electronic iter to output charge density and wavefunction. 0: output only when converged");
    ModuleBase::GlobalFunc::OUTP(ofs, "dft_plus_dmft", dft_plus_dmft, "true:DFT+DMFT; false: standard DFT calcullation(default)");
    ModuleBase::GlobalFunc::OUTP(ofs, "printe", printe, "Print out energy for each band for every printe steps");
    ModuleBase::GlobalFunc::OUTP(ofs, "mem_saver", mem_saver, "Only for nscf calculations. if set to 1, then a memory saving technique will be used for many k point calculations.");
    ModuleBase::GlobalFunc::OUTP(ofs, "diago_proc", diago_proc, "the number of procs used to do diagonalization");
    ModuleBase::GlobalFunc::OUTP(ofs, "nbspline", nbspline, "the order of B-spline basis");
    ModuleBase::GlobalFunc::OUTP(ofs, "wannier_card", wannier_card, "input card for wannier functions");
    ModuleBase::GlobalFunc::OUTP(ofs, "soc_lambda", soc_lambda, "The fraction of averaged SOC pseudopotential is given by (1-soc_lambda)");
    ModuleBase::GlobalFunc::OUTP(ofs, "cal_force", cal_force, "if calculate the force at the end of the electronic iteration");
    ModuleBase::GlobalFunc::OUTP(ofs, "out_freq_ion", out_freq_ion, "the frequency ( >= 0 ) of ionic step to output charge density and wavefunction. 0: output only when ion steps are finished");

    ofs << "\n#Parameters (2.PW)" << std::endl;
    ModuleBase::GlobalFunc::OUTP(ofs, "ecutwfc", ecutwfc, "#energy cutoff for wave functions");
    if (ks_solver == "cg")
    {
        ModuleBase::GlobalFunc::OUTP(ofs, "pw_diag_nmax", pw_diag_nmax, "max iteration number for cg");
        ModuleBase::GlobalFunc::OUTP(ofs, "diago_cg_prec", diago_cg_prec, "diago_cg_prec");
    }
    else if (ks_solver == "dav")
    {
        ModuleBase::GlobalFunc::OUTP(ofs, "pw_diag_ndim", pw_diag_ndim, "max dimension for davidson");
    }
    ModuleBase::GlobalFunc::OUTP(ofs,
                                 "pw_diag_thr",
                                 pw_diag_thr,
                                 "threshold for eigenvalues is cg electron iterations");
    ModuleBase::GlobalFunc::OUTP(ofs, "scf_thr", scf_thr, "charge density error");
    ModuleBase::GlobalFunc::OUTP(ofs, "init_wfc", init_wfc, "start wave functions are from 'atomic', 'atomic+random', 'random' or 'file'");
    ModuleBase::GlobalFunc::OUTP(ofs, "init_chg", init_chg, "start charge is from 'atomic' or file");
    ModuleBase::GlobalFunc::OUTP(ofs,
                                 "chg_extrap",
                                 chg_extrap,
                                 "atomic; first-order; second-order; dm:coefficients of SIA");
    ModuleBase::GlobalFunc::OUTP(ofs, "out_chg", out_chg, ">0 output charge density for selected electron steps");
    ModuleBase::GlobalFunc::OUTP(ofs, "out_pot", out_pot, "output realspace potential");
    ModuleBase::GlobalFunc::OUTP(ofs, "out_wfc_pw", out_wfc_pw, "output wave functions");
    ModuleBase::GlobalFunc::OUTP(ofs, "out_wfc_r", out_wfc_r, "output wave functions in realspace");
    ModuleBase::GlobalFunc::OUTP(ofs, "out_dos", out_dos, "output energy and dos");
    ModuleBase::GlobalFunc::OUTP(ofs, "out_band", out_band, "output energy and band structure");
    ModuleBase::GlobalFunc::OUTP(ofs, "out_proj_band", out_proj_band, "output projected band structure");
    ModuleBase::GlobalFunc::OUTP(ofs, "restart_save", restart_save, "print to disk every step for restart");
    ModuleBase::GlobalFunc::OUTP(ofs, "restart_load", restart_load, "restart from disk");
    ModuleBase::GlobalFunc::OUTP(ofs, "read_file_dir", read_file_dir, "directory of files for reading");
    ModuleBase::GlobalFunc::OUTP(ofs, "nx", nx, "number of points along x axis for FFT grid");
    ModuleBase::GlobalFunc::OUTP(ofs, "ny", ny, "number of points along y axis for FFT grid");
    ModuleBase::GlobalFunc::OUTP(ofs, "nz", nz, "number of points along z axis for FFT grid");
    ModuleBase::GlobalFunc::OUTP(ofs,
                                 "cell_factor",
                                 cell_factor,
                                 "used in the construction of the pseudopotential tables");
    ModuleBase::GlobalFunc::OUTP(ofs, "pw_seed", pw_seed, "random seed for initializing wave functions");
    
    ofs << "\n#Parameters (3.Stochastic DFT)" << std::endl;
    ModuleBase::GlobalFunc::OUTP(ofs, "method_sto", method_sto, "1: slow and save memory, 2: fast and waste memory");
    ModuleBase::GlobalFunc::OUTP(ofs, "npart_sto", npart_sto, "Reduce memory when calculating Stochastic DOS");
    ModuleBase::GlobalFunc::OUTP(ofs, "nbands_sto", nbands_sto, "number of stochstic orbitals");
    ModuleBase::GlobalFunc::OUTP(ofs, "nche_sto", nche_sto, "Chebyshev expansion orders");
    ModuleBase::GlobalFunc::OUTP(ofs, "emin_sto", emin_sto, "trial energy to guess the lower bound of eigen energies of the Hamitonian operator");
    ModuleBase::GlobalFunc::OUTP(ofs, "emax_sto", emax_sto, "trial energy to guess the upper bound of eigen energies of the Hamitonian operator");
    ModuleBase::GlobalFunc::OUTP(ofs, "seed_sto", seed_sto, "the random seed to generate stochastic orbitals");
    ModuleBase::GlobalFunc::OUTP(ofs, "initsto_freq", initsto_freq, "frequency to generate new stochastic orbitals when running md");
    ModuleBase::GlobalFunc::OUTP(ofs, "cal_cond", cal_cond, "calculate electronic conductivities");
    ModuleBase::GlobalFunc::OUTP(ofs, "cond_nche", cond_nche, "orders of Chebyshev expansions for conductivities");
    ModuleBase::GlobalFunc::OUTP(ofs, "cond_dw", cond_dw, "frequency interval for conductivities");
    ModuleBase::GlobalFunc::OUTP(ofs, "cond_wcut", cond_wcut, "cutoff frequency (omega) for conductivities");
    ModuleBase::GlobalFunc::OUTP(ofs, "cond_wenlarge", cond_wenlarge, "control the t interval: dt = PI/wcut/cond_wenlarge");
    ModuleBase::GlobalFunc::OUTP(ofs, "cond_fwhm", cond_fwhm, "FWHM for conductivities");
    ModuleBase::GlobalFunc::OUTP(ofs, "cond_nonlocal", cond_nonlocal, "Nonlocal effects for conductivities");
    
    ofs << "\n#Parameters (4.Relaxation)" << std::endl;
    ModuleBase::GlobalFunc::OUTP(ofs,
                                 "ks_solver",
                                 GlobalV::KS_SOLVER,
                                 "cg; dav; lapack; genelpa; hpseps; scalapack_gvx; cusolver");
    ModuleBase::GlobalFunc::OUTP(ofs, "scf_nmax", scf_nmax, "#number of electron iterations");
    ModuleBase::GlobalFunc::OUTP(ofs, "out_force", out_force, "output the out_force or not");
    ModuleBase::GlobalFunc::OUTP(ofs, "relax_nmax", relax_nmax, "number of ion iteration steps");
    ModuleBase::GlobalFunc::OUTP(ofs, "out_stru", out_stru, "output the structure files after each ion step");
    ModuleBase::GlobalFunc::OUTP(ofs, "force_thr", force_thr, "force threshold, unit: Ry/Bohr");
    ModuleBase::GlobalFunc::OUTP(ofs,
                                 "force_thr_ev",
                                 force_thr * 13.6058 / 0.529177,
                                 "force threshold, unit: eV/Angstrom");
    ModuleBase::GlobalFunc::OUTP(ofs, "force_thr_ev2", force_thr_ev2, "force invalid threshold, unit: eV/Angstrom");
    ModuleBase::GlobalFunc::OUTP(ofs,
                                 "relax_cg_thr",
                                 relax_cg_thr,
                                 "threshold for switching from cg to bfgs, unit: eV/Angstrom");
    ModuleBase::GlobalFunc::OUTP(ofs, "stress_thr", stress_thr, "stress threshold");
    ModuleBase::GlobalFunc::OUTP(ofs, "press1", press1, "target pressure, unit: KBar");
    ModuleBase::GlobalFunc::OUTP(ofs, "press2", press2, "target pressure, unit: KBar");
    ModuleBase::GlobalFunc::OUTP(ofs, "press3", press3, "target pressure, unit: KBar");
    ModuleBase::GlobalFunc::OUTP(ofs, "relax_bfgs_w1", relax_bfgs_w1, "wolfe condition 1 for bfgs");
    ModuleBase::GlobalFunc::OUTP(ofs, "relax_bfgs_w2", relax_bfgs_w2, "wolfe condition 2 for bfgs");
    ModuleBase::GlobalFunc::OUTP(ofs, "relax_bfgs_rmax", relax_bfgs_rmax, "maximal trust radius, unit: Bohr");
    ModuleBase::GlobalFunc::OUTP(ofs, "relax_bfgs_rmin", relax_bfgs_rmin, "minimal trust radius, unit: Bohr");
    ModuleBase::GlobalFunc::OUTP(ofs, "relax_bfgs_init", relax_bfgs_init, "initial trust radius, unit: Bohr");
    ModuleBase::GlobalFunc::OUTP(ofs, "cal_stress", cal_stress, "calculate the stress or not");
    ModuleBase::GlobalFunc::OUTP(ofs, "fixed_axes", fixed_axes, "which axes are fixed");
    ModuleBase::GlobalFunc::OUTP(ofs, "relax_method", relax_method, "bfgs; sd; cg; cg_bfgs;"); // pengfei add 2013-08-15
    ModuleBase::GlobalFunc::OUTP(ofs, "out_level", out_level, "ie(for electrons); i(for ions);");
    ModuleBase::GlobalFunc::OUTP(ofs, "out_dm", out_dm, ">0 output density matrix");

    // for deepks
    ModuleBase::GlobalFunc::OUTP(ofs, "deepks_out_labels", deepks_out_labels, ">0 compute descriptor for deepks");
    ModuleBase::GlobalFunc::OUTP(ofs, "deepks_scf", deepks_scf, ">0 add V_delta to Hamiltonian");
    ModuleBase::GlobalFunc::OUTP(ofs, "deepks_bandgap", deepks_bandgap, ">0 for bandgap label");
    ModuleBase::GlobalFunc::OUTP(ofs,
                                 "deepks_out_unittest",
                                 deepks_out_unittest,
                                 "if set 1, prints intermediate quantities that shall be used for making unit test");
    ModuleBase::GlobalFunc::OUTP(ofs, "deepks_model", deepks_model, "file dir of traced pytorch model: 'model.ptg");

    ModuleBase::GlobalFunc::OUTP(ofs,
                                 "bessel_lmax",
                                 bessel_lmax,
                                 "lmax used in generating bessel functions");
    ModuleBase::GlobalFunc::OUTP(ofs,
                                 "bessel_rcut",
                                 bessel_rcut,
                                 "rcut used in generating bessel functions");
    ModuleBase::GlobalFunc::OUTP(ofs,
                                 "bessel_rcut",
                                 bessel_rcut,
                                 "tolerence level when generating bessel functions");

    ofs << "\n#Parameters (5.LCAO)" << std::endl;
    ModuleBase::GlobalFunc::OUTP(ofs, "basis_type", basis_type, "PW; LCAO in pw; LCAO");
    if (ks_solver == "HPSEPS" || ks_solver == "genelpa" || ks_solver == "scalapack_gvx" || ks_solver == "cusolver")
    {
        ModuleBase::GlobalFunc::OUTP(ofs, "nb2d", nb2d, "2d distribution of atoms");
    }
    ModuleBase::GlobalFunc::OUTP(ofs, "gamma_only", gamma_only, "Only for localized orbitals set and gamma point. If set to 1, a fast algorithm is used");
    ModuleBase::GlobalFunc::OUTP(ofs, "search_radius", search_radius, "input search radius (Bohr)");
    ModuleBase::GlobalFunc::OUTP(ofs, "search_pbc", search_pbc, "input periodic boundary condition");
    ModuleBase::GlobalFunc::OUTP(ofs, "lcao_ecut", lcao_ecut, "energy cutoff for LCAO");
    ModuleBase::GlobalFunc::OUTP(ofs, "lcao_dk", lcao_dk, "delta k for 1D integration in LCAO");
    ModuleBase::GlobalFunc::OUTP(ofs, "lcao_dr", lcao_dr, "delta r for 1D integration in LCAO");
    ModuleBase::GlobalFunc::OUTP(ofs, "lcao_rmax", lcao_rmax, "max R for 1D two-center integration table");
    ModuleBase::GlobalFunc::OUTP(ofs, "out_mat_hs", out_mat_hs, "output H and S matrix");
    ModuleBase::GlobalFunc::OUTP(ofs, "out_mat_hs2", out_mat_hs2, "output H(R) and S(R) matrix");
    ModuleBase::GlobalFunc::OUTP(ofs, "out_element_info", out_element_info, "output (projected) wavefunction of each element");
    ModuleBase::GlobalFunc::OUTP(ofs, "out_mat_r", out_mat_r, "output r(R) matrix");
    ModuleBase::GlobalFunc::OUTP(ofs, "out_wfc_lcao", out_wfc_lcao, "ouput LCAO wave functions");
    ModuleBase::GlobalFunc::OUTP(ofs, "bx", bx, "division of an element grid in FFT grid along x");
    ModuleBase::GlobalFunc::OUTP(ofs, "by", by, "division of an element grid in FFT grid along y");
    ModuleBase::GlobalFunc::OUTP(ofs, "bz", bz, "division of an element grid in FFT grid along z");

    ofs << "\n#Parameters (6.Smearing)" << std::endl;
    ModuleBase::GlobalFunc::OUTP(ofs,
                                 "smearing_method",
                                 smearing_method,
                                 "type of smearing_method: gauss; fd; fixed; mp; mp2; mv");
    ModuleBase::GlobalFunc::OUTP(ofs, "smearing_sigma", smearing_sigma, "energy range for smearing");

    ofs << "\n#Parameters (7.Charge Mixing)" << std::endl;
    ModuleBase::GlobalFunc::OUTP(ofs, "mixing_type", mixing_mode, "plain; kerker; pulay; pulay-kerker; broyden");
    ModuleBase::GlobalFunc::OUTP(ofs, "mixing_beta", mixing_beta, "mixing parameter: 0 means no new charge");
    ModuleBase::GlobalFunc::OUTP(ofs, "mixing_ndim", mixing_ndim, "mixing dimension in pulay");
    ModuleBase::GlobalFunc::OUTP(ofs, "mixing_gg0", mixing_gg0, "mixing parameter in kerker");

    ofs << "\n#Parameters (8.DOS)" << std::endl;
    ModuleBase::GlobalFunc::OUTP(ofs, "dos_emin_ev", dos_emin_ev, "minimal range for dos");
    ModuleBase::GlobalFunc::OUTP(ofs, "dos_emax_ev", dos_emax_ev, "maximal range for dos");
    ModuleBase::GlobalFunc::OUTP(ofs, "dos_edelta_ev", dos_edelta_ev, "delta energy for dos");
    ModuleBase::GlobalFunc::OUTP(ofs, "dos_scale", dos_scale, "scale dos range by");
    ModuleBase::GlobalFunc::OUTP(ofs, "dos_sigma", b_coef, "gauss b coefficeinet(default=0.07)");
    ModuleBase::GlobalFunc::OUTP(ofs, "dos_nche", dos_nche, "orders of Chebyshev expansions for dos");

	ofs << "\n#Parameters (9.Molecular dynamics)" << std::endl;
	ModuleBase::GlobalFunc::OUTP(ofs,"md_type",mdp.md_type,"choose ensemble");
    ModuleBase::GlobalFunc::OUTP(ofs,"md_thermostat",mdp.md_thermostat,"choose thermostat");
    ModuleBase::GlobalFunc::OUTP(ofs,"md_nstep",mdp.md_nstep,"md steps");
	ModuleBase::GlobalFunc::OUTP(ofs,"md_ensolver",mdp.md_ensolver,"choose potential");
	ModuleBase::GlobalFunc::OUTP(ofs,"md_dt",mdp.md_dt,"time step");
	ModuleBase::GlobalFunc::OUTP(ofs,"md_mnhc",mdp.md_mnhc,"number of Nose-Hoover chains");
	ModuleBase::GlobalFunc::OUTP(ofs,"md_tfirst",mdp.md_tfirst,"temperature first");
	ModuleBase::GlobalFunc::OUTP(ofs,"md_tlast",mdp.md_tlast,"temperature last");
	ModuleBase::GlobalFunc::OUTP(ofs,"md_dumpfreq",mdp.md_dumpfreq,"The period to dump MD information");
	ModuleBase::GlobalFunc::OUTP(ofs,"md_restartfreq",mdp.md_restartfreq,"The period to output MD restart information");
    ModuleBase::GlobalFunc::OUTP(ofs,"md_seed",mdp.md_seed,"random seed for MD");
	ModuleBase::GlobalFunc::OUTP(ofs,"md_restart",mdp.md_restart,"whether restart");
	ModuleBase::GlobalFunc::OUTP(ofs,"lj_rcut",mdp.lj_rcut,"cutoff radius of LJ potential");
	ModuleBase::GlobalFunc::OUTP(ofs,"lj_epsilon",mdp.lj_epsilon,"the value of epsilon for LJ potential");
	ModuleBase::GlobalFunc::OUTP(ofs,"lj_sigma",mdp.lj_sigma,"the value of sigma for LJ potential");
    ModuleBase::GlobalFunc::OUTP(ofs,"pot_file",mdp.pot_file,"the filename of potential files for CMD such as DP");
	ModuleBase::GlobalFunc::OUTP(ofs,"msst_direction",mdp.msst_direction,"the direction of shock wave");
	ModuleBase::GlobalFunc::OUTP(ofs,"msst_vel",mdp.msst_vel,"the velocity of shock wave");
	ModuleBase::GlobalFunc::OUTP(ofs,"msst_vis",mdp.msst_vis,"artificial viscosity");
	ModuleBase::GlobalFunc::OUTP(ofs,"msst_tscale",mdp.msst_tscale,"reduction in initial temperature");
	ModuleBase::GlobalFunc::OUTP(ofs,"msst_qmass",mdp.msst_qmass,"mass of thermostat");
	ModuleBase::GlobalFunc::OUTP(ofs,"md_tfreq",mdp.md_tfreq,"oscillation frequency, used to determine qmass of NHC");
	ModuleBase::GlobalFunc::OUTP(ofs,"md_damp",mdp.md_damp,"damping parameter (time units) used to add force in Langevin method");
    ModuleBase::GlobalFunc::OUTP(ofs,"md_nraise",mdp.md_nraise,"parameters used when md_type=0");
    ModuleBase::GlobalFunc::OUTP(ofs,"md_tolerance",mdp.md_tolerance,"tolerance for velocity rescaling (K)");

    ofs << "\n#Parameters (10.Electric field and dipole correction)" << std::endl;
    ModuleBase::GlobalFunc::OUTP(ofs,"efield_flag",efield_flag,"add electric field");
    ModuleBase::GlobalFunc::OUTP(ofs,"dip_cor_flag",dip_cor_flag,"dipole correction");
    ModuleBase::GlobalFunc::OUTP(ofs,"efield_dir",efield_dir,"the direction of the electric field or dipole correction");
    ModuleBase::GlobalFunc::OUTP(ofs,"efield_pos_max",efield_pos_max,"position of the maximum of the saw-like potential along crystal axis efield_dir");
    ModuleBase::GlobalFunc::OUTP(ofs,"efield_pos_dec",efield_pos_dec,"zone in the unit cell where the saw-like potential decreases");
    ModuleBase::GlobalFunc::OUTP(ofs,"efield_amp ",efield_amp ,"amplitude of the electric field");

    ofs << "\n#Parameters (11.Gate field)" << std::endl;
    ModuleBase::GlobalFunc::OUTP(ofs,"gate_flag",gate_flag,"compensating charge or not");
    ModuleBase::GlobalFunc::OUTP(ofs,"zgate",zgate,"position of charged plate");
    ModuleBase::GlobalFunc::OUTP(ofs,"relax",relax,"allow relaxation along the specific direction");
    ModuleBase::GlobalFunc::OUTP(ofs,"block",block,"add a block potential or not");
    ModuleBase::GlobalFunc::OUTP(ofs,"block_down",block_down,"low bound of the block");
    ModuleBase::GlobalFunc::OUTP(ofs,"block_up",block_up,"high bound of the block");
    ModuleBase::GlobalFunc::OUTP(ofs,"block_height",block_height,"height of the block");

    ofs << "\n#Parameters (12.Test)" << std::endl;
    ModuleBase::GlobalFunc::OUTP(ofs, "out_alllog", out_alllog, "output information for each processor, when parallel");
    ModuleBase::GlobalFunc::OUTP(ofs, "nurse", nurse, "for coders");
    ModuleBase::GlobalFunc::OUTP(ofs, "colour", colour, "for coders, make their live colourful");
    ModuleBase::GlobalFunc::OUTP(ofs, "t_in_h", t_in_h, "calculate the kinetic energy or not");
    ModuleBase::GlobalFunc::OUTP(ofs, "vl_in_h", vl_in_h, "calculate the local potential or not");
    ModuleBase::GlobalFunc::OUTP(ofs, "vnl_in_h", vnl_in_h, "calculate the nonlocal potential or not");
    ModuleBase::GlobalFunc::OUTP(ofs, "vh_in_h", vh_in_h, "calculate the hartree potential or not");
    ModuleBase::GlobalFunc::OUTP(ofs, "vion_in_h", vion_in_h, "calculate the local ionic potential or not");
    ModuleBase::GlobalFunc::OUTP(ofs, "test_force", test_force, "test the force");
    ModuleBase::GlobalFunc::OUTP(ofs, "test_stress", test_stress, "test the force");
    ModuleBase::GlobalFunc::OUTP(ofs, "test_skip_ewald", test_skip_ewald, "skip ewald energy");

    ofs << "\n#Parameters (13.vdW Correction)" << std::endl;
    ModuleBase::GlobalFunc::OUTP(ofs,
                                 "vdw_method",
                                 vdw_method,
                                 "the method of calculating vdw (none ; d2 ; d3_0 ; d3_bj");
    ModuleBase::GlobalFunc::OUTP(ofs, "vdw_s6", vdw_s6, "scale parameter of d2/d3_0/d3_bj");
    ModuleBase::GlobalFunc::OUTP(ofs, "vdw_s8", vdw_s8, "scale parameter of d3_0/d3_bj");
    ModuleBase::GlobalFunc::OUTP(ofs, "vdw_a1", vdw_a1, "damping parameter of d3_0/d3_bj");
    ModuleBase::GlobalFunc::OUTP(ofs, "vdw_a2", vdw_a2, "damping parameter of d3_bj");
    ModuleBase::GlobalFunc::OUTP(ofs, "vdw_d", vdw_d, "damping parameter of d2");
    ModuleBase::GlobalFunc::OUTP(ofs, "vdw_abc", vdw_abc, "third-order term?");
    ModuleBase::GlobalFunc::OUTP(ofs, "vdw_C6_file", vdw_C6_file, "filename of C6");
    ModuleBase::GlobalFunc::OUTP(ofs, "vdw_C6_unit", vdw_C6_unit, "unit of C6, Jnm6/mol or eVA6");
    ModuleBase::GlobalFunc::OUTP(ofs, "vdw_R0_file", vdw_R0_file, "filename of R0");
    ModuleBase::GlobalFunc::OUTP(ofs, "vdw_R0_unit", vdw_R0_unit, "unit of R0, A or Bohr");
    ModuleBase::GlobalFunc::OUTP(ofs,
                                 "vdw_cutoff_type",
                                 vdw_cutoff_type,
                                 "expression model of periodic structure, radius or period");
    ModuleBase::GlobalFunc::OUTP(ofs, "vdw_cutoff_radius", vdw_cutoff_radius, "radius cutoff for periodic structure");
    ModuleBase::GlobalFunc::OUTP(ofs,
                                 "vdw_radius_unit",
                                 vdw_radius_unit,
                                 "unit of radius cutoff for periodic structure");
    ModuleBase::GlobalFunc::OUTP(ofs, "vdw_cn_thr", vdw_cn_thr, "radius cutoff for cn");
    ModuleBase::GlobalFunc::OUTP(ofs, "vdw_cn_thr_unit", vdw_cn_thr_unit, "unit of cn_thr, Bohr or Angstrom");
    ofs << std::setw(20) << "vdw_cutoff_period" << vdw_cutoff_period.x << " " << vdw_cutoff_period.y << " " << vdw_cutoff_period.z
        << " #periods of periodic structure" << std::endl;

    ofs << "\n#Parameters (14.exx)" << std::endl;
    ModuleBase::GlobalFunc::OUTP(ofs, "exx_hybrid_alpha", exx_hybrid_alpha, "");
    ModuleBase::GlobalFunc::OUTP(ofs, "exx_hse_omega", exx_hse_omega, "");
    ModuleBase::GlobalFunc::OUTP(ofs, "exx_separate_loop", exx_separate_loop, "0 or 1");
    ModuleBase::GlobalFunc::OUTP(ofs, "exx_hybrid_step", exx_hybrid_step, "");
    ModuleBase::GlobalFunc::OUTP(ofs, "exx_lambda", exx_lambda, "");
    ModuleBase::GlobalFunc::OUTP(ofs, "exx_pca_threshold", exx_pca_threshold, "");
    ModuleBase::GlobalFunc::OUTP(ofs, "exx_c_threshold", exx_c_threshold, "");
    ModuleBase::GlobalFunc::OUTP(ofs, "exx_v_threshold", exx_v_threshold, "");
    ModuleBase::GlobalFunc::OUTP(ofs, "exx_dm_threshold", exx_dm_threshold, "");
    ModuleBase::GlobalFunc::OUTP(ofs, "exx_schwarz_threshold", exx_schwarz_threshold, "");
    ModuleBase::GlobalFunc::OUTP(ofs, "exx_cauchy_threshold", exx_cauchy_threshold, "");
    ModuleBase::GlobalFunc::OUTP(ofs, "exx_ccp_threshold", exx_ccp_threshold, "");
    ModuleBase::GlobalFunc::OUTP(ofs, "exx_ccp_rmesh_times", exx_ccp_rmesh_times, "");
    ModuleBase::GlobalFunc::OUTP(ofs, "exx_distribute_type", exx_distribute_type, "htime or kmeans1 or kmeans2");
    ModuleBase::GlobalFunc::OUTP(ofs, "exx_opt_orb_lmax", exx_opt_orb_lmax, "");
    ModuleBase::GlobalFunc::OUTP(ofs, "exx_opt_orb_ecut", exx_opt_orb_ecut, "");
    ModuleBase::GlobalFunc::OUTP(ofs, "exx_opt_orb_tolerence", exx_opt_orb_tolerence, "");

    ofs << "\n#Parameters (16.tddft)" << std::endl;
    ModuleBase::GlobalFunc::OUTP(ofs, "tddft", tddft, "calculate tddft or not");
    ModuleBase::GlobalFunc::OUTP(ofs, "td_scf_thr", td_scf_thr, "threshold for electronic iteration of tddft");
    ModuleBase::GlobalFunc::OUTP(ofs, "td_dt", td_dt, "time of ion step");
    ModuleBase::GlobalFunc::OUTP(ofs, "td_force_dt", td_force_dt, "time of force change");
    ModuleBase::GlobalFunc::OUTP(ofs, "td_val_elec_01", td_val_elec_01, "td_val_elec_01");
    ModuleBase::GlobalFunc::OUTP(ofs, "td_val_elec_02", td_val_elec_02, "td_val_elec_02");
    ModuleBase::GlobalFunc::OUTP(ofs, "td_val_elec_03", td_val_elec_03, "td_val_elec_03");
    ModuleBase::GlobalFunc::OUTP(ofs, "td_vext", td_vext, "add extern potential or not");
    ModuleBase::GlobalFunc::OUTP(ofs, "td_vext_dire", td_vext_dire, "extern potential direction");
    ModuleBase::GlobalFunc::OUTP(ofs, "td_timescale", td_timescale, "extern potential td_timescale");
    ModuleBase::GlobalFunc::OUTP(ofs, "td_vexttype", td_vexttype, "extern potential type");
    ModuleBase::GlobalFunc::OUTP(ofs, "td_vextout", td_vextout, "output extern potential or not");
    ModuleBase::GlobalFunc::OUTP(ofs, "td_dipoleout", td_dipoleout, "output dipole or not");
    ModuleBase::GlobalFunc::OUTP(ofs, "ocp", GlobalV::ocp, "change occupation or not");
    ModuleBase::GlobalFunc::OUTP(ofs, "ocp_set", GlobalV::ocp_set, "set occupation");

    ofs << "\n#Parameters (17.berry_wannier)" << std::endl;
    ModuleBase::GlobalFunc::OUTP(ofs, "berry_phase", berry_phase, "calculate berry phase or not");
    ModuleBase::GlobalFunc::OUTP(ofs,
                                 "gdir",
                                 gdir,
                                 "calculate the polarization in the direction of the lattice vector");
    ModuleBase::GlobalFunc::OUTP(ofs, "towannier90", towannier90, "use wannier90 code interface or not");
    ModuleBase::GlobalFunc::OUTP(ofs, "nnkpfile", NNKP, "the wannier90 code nnkp file name");
    ModuleBase::GlobalFunc::OUTP(ofs, "wannier_spin", wannier_spin, "calculate spin in wannier90 code interface");

    ofs << "\n#Parameters (18.implicit_solvation)" << std::endl;
    ModuleBase::GlobalFunc::OUTP(ofs, "imp_sol", imp_sol, "calculate implicit solvation correction or not");
    ModuleBase::GlobalFunc::OUTP(ofs, "eb_k", eb_k, "the relative permittivity of the bulk solvent");
    ModuleBase::GlobalFunc::OUTP(ofs, "tau", tau, "the effective surface tension parameter");
    ModuleBase::GlobalFunc::OUTP(ofs, "sigma_k", sigma_k, " the width of the diffuse cavity");
    ModuleBase::GlobalFunc::OUTP(ofs, "nc_k", nc_k, " the cut-off charge density");

    ofs << "\n#Parameters (19.orbital free density functional theory)" << std::endl;
    ModuleBase::GlobalFunc::OUTP(ofs, "of_kinetic", of_kinetic, "kinetic energy functional, such as tf, vw, wt");
    ModuleBase::GlobalFunc::OUTP(ofs, "of_method", of_method, "optimization method used in OFDFT, including cg1, cg2, tn (default)");
    ModuleBase::GlobalFunc::OUTP(ofs, "of_conv", of_conv, "the convergence criterion, potential, energy (default), or both");
    ModuleBase::GlobalFunc::OUTP(ofs, "of_tole", of_tole, "tolerance of the energy change (in Ry) for determining the convergence, default=2e-6 Ry");
    ModuleBase::GlobalFunc::OUTP(ofs, "of_tolp", of_tolp, "tolerance of potential for determining the convergence, default=1e-5 in a.u.");
    ModuleBase::GlobalFunc::OUTP(ofs, "of_tf_weight", of_tf_weight, "weight of TF KEDF");
    ModuleBase::GlobalFunc::OUTP(ofs, "of_vw_weight", of_vw_weight, "weight of vW KEDF");
    ModuleBase::GlobalFunc::OUTP(ofs, "of_wt_alpha", of_wt_alpha, "parameter alpha of WT KEDF");
    ModuleBase::GlobalFunc::OUTP(ofs, "of_wt_beta", of_wt_beta, "parameter beta of WT KEDF");
    ModuleBase::GlobalFunc::OUTP(ofs, "of_wt_rho0", of_wt_rho0, "the average density of system, used in WT KEDF, in Bohr^-3");
    ModuleBase::GlobalFunc::OUTP(ofs, "of_hold_rho0", of_hold_rho0, "If set to 1, the rho0 will be fixed even if the volume of system has changed, it will be set to 1 automaticly if of_wt_rho0 is not zero");
    ModuleBase::GlobalFunc::OUTP(ofs, "of_full_pw", of_full_pw, "If set to 1, ecut will be ignored when collect planewaves, so that all planewaves will be used");
    ModuleBase::GlobalFunc::OUTP(ofs, "of_full_pw_dim", of_full_pw_dim, "If of_full_pw = true, dimention of FFT is testricted to be (0) either odd or even; (1) odd only; (2) even only");
    ModuleBase::GlobalFunc::OUTP(ofs, "of_read_kernel", of_read_kernel, "If set to 1, the kernel of WT KEDF will be filled from file of_kernel_file, not from formula. Only usable for WT KEDF");
    ModuleBase::GlobalFunc::OUTP(ofs, "of_kernel_file", of_kernel_file, "The name of WT kernel file.");

    ofs << "\n#Parameters (19.dft+u)" << std::endl;
    ModuleBase::GlobalFunc::OUTP(ofs, "dft_plus_u", dft_plus_u, "true:DFT+U correction; false: standard DFT calcullation(default)");
    ModuleBase::GlobalFunc::OUTP(ofs, "dftu_type", dftu_type, "1:rotationally invarient formalism; 2:simplified form(default)");
    ModuleBase::GlobalFunc::OUTP(ofs, "double_counting", double_counting, "1:FLL(fully localized limit)(default); 2:AMF(around mean field)");
    ModuleBase::GlobalFunc::OUTP(ofs, "yukawa_lambda", yukawa_lambda, "default:0.0");
    ModuleBase::GlobalFunc::OUTP(ofs, "yukawa_potential", yukawa_potential, "default: false");
    ModuleBase::GlobalFunc::OUTP(ofs, "omc", omc, "whether turn on occupation matrix control method or not");
    ofs << std::setw(20) << "hund_j ";
    for (int i = 0; i < ntype; i++)
    {
        ofs << hund_j[i]*ModuleBase::Ry_to_eV << " ";
    }
    ofs << "#Hund exchange parameter J(ev)" << std::endl;
    ofs << std::setw(20) << "hubbard_u ";
    for (int i = 0; i < ntype; i++)
    {
        ofs << hubbard_u[i]*ModuleBase::Ry_to_eV << " ";
    }
    ofs << "#Hubbard Coulomb interaction parameter U(ev)" << std::endl;
    ofs << std::setw(20) << "orbital_corr ";
    for (int i = 0; i < ntype; i++)
    {
        ofs << orbital_corr[i] << " ";
    }
    ofs << "#which correlated orbitals need corrected ; d:2 ,f:3, do not need correction:-1" << std::endl;


    ofs.close();
    return;
}
