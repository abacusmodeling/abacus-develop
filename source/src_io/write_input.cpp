#include "../input.h"
#include "../src_pw/tools.h"

void Input::Print(const std::string &fn)const
{
    if (GlobalV::MY_RANK!=0) return;

    ModuleBase::TITLE("Input","Print");

    std::ofstream ofs(fn.c_str());

	//----------------------------------
	// output the information in INPUT.
	//----------------------------------
    ofs << "INPUT_PARAMETERS" << std::endl;
	ofs << std::setiosflags(ios::left);
	
	ofs << "#Parameters (1.General)" << std::endl;
	ModuleBase::GlobalFunc::OUTP(ofs,"suffix",suffix,"the name of main output directory");
	ModuleBase::GlobalFunc::OUTP(ofs,"latname",latname,"the name of lattice name");
	ModuleBase::GlobalFunc::OUTP(ofs,"atom_file",GlobalV::global_atom_card,"the filename of file containing atom positions");
	ModuleBase::GlobalFunc::OUTP(ofs,"kpoint_file",GlobalV::global_kpoint_card,"the name of file containing k points");
	ModuleBase::GlobalFunc::OUTP(ofs,"pseudo_dir",GlobalV::global_pseudo_dir,"the directory containing pseudo files");
	ModuleBase::GlobalFunc::OUTP(ofs,"orbital_dir",GlobalV::global_orbital_dir,"the directory containing orbital files");
	ModuleBase::GlobalFunc::OUTP(ofs,"pseudo_type",GlobalV::global_pseudo_type,"the type pseudo files");
	ModuleBase::GlobalFunc::OUTP(ofs,"pseudo_rcut",pseudo_rcut,"cut-off radius for radial integration");
	ModuleBase::GlobalFunc::OUTP(ofs,"renormwithmesh",renormwithmesh,"0: use our own mesh to do radial renormalization; 1: use mesh as in QE");
	ModuleBase::GlobalFunc::OUTP(ofs,"lmaxmax",lmaxmax,"maximum of l channels used");
	ModuleBase::GlobalFunc::OUTP(ofs,"dft_functional",dft_functional,"exchange correlation functional");
	ModuleBase::GlobalFunc::OUTP(ofs,"calculation",calculation,"test; scf; relax; nscf; ienvelope; istate;");
	ModuleBase::GlobalFunc::OUTP(ofs,"ntype",ntype,"atom species number");
	ModuleBase::GlobalFunc::OUTP(ofs,"nspin",nspin,"1: single spin; 2: up and down spin; 4: noncollinear spin");
	ModuleBase::GlobalFunc::OUTP(ofs,"nbands",nbands,"number of bands");
	ModuleBase::GlobalFunc::OUTP(ofs,"nbands_sto",nbands_sto,"number of stochastic bands");
	ModuleBase::GlobalFunc::OUTP(ofs,"nbands_istate",nbands_istate,"number of bands around Fermi level for istate calulation");
	ModuleBase::GlobalFunc::OUTP(ofs,"nche_sto",nche_sto,"number of orders for Chebyshev expansion in stochastic DFT");
	ModuleBase::GlobalFunc::OUTP(ofs,"symmetry",symmetry,"turn symmetry on or off");	
	ModuleBase::GlobalFunc::OUTP(ofs,"set_vel",set_vel,"read velocity from STRU or not");
	ModuleBase::GlobalFunc::OUTP(ofs,"symmetry_prec",symmetry_prec,"accuracy for symmetry"); // LiuXh add 2021-08-12, accuracy for symmetry
	ModuleBase::GlobalFunc::OUTP(ofs,"nelec",nelec,"input number of electrons");
	ModuleBase::GlobalFunc::OUTP(ofs,"tot_magnetization",tot_magnetization,"total magnetization of the system");

	ofs << "\n#Parameters (2.PW)" << std::endl;
	ModuleBase::GlobalFunc::OUTP(ofs,"ecutwfc",ecutwfc,"#energy cutoff for wave functions");
	if(ks_solver=="cg")
	{
		ModuleBase::GlobalFunc::OUTP(ofs,"diago_cg_maxiter",diago_cg_maxiter,"max iteration number for cg");
		ModuleBase::GlobalFunc::OUTP(ofs,"diago_cg_prec",diago_cg_prec,"diago_cg_prec");
	}
	else if(ks_solver=="dav")
	{
		ModuleBase::GlobalFunc::OUTP(ofs,"diago_david_ndim",diago_david_ndim,"max dimension for davidson");
	}
	ModuleBase::GlobalFunc::OUTP(ofs,"ethr",ethr,"threshold for eigenvalues is cg electron iterations");
	ModuleBase::GlobalFunc::OUTP(ofs,"dr2",dr2,"charge density error");
	ModuleBase::GlobalFunc::OUTP(ofs,"start_wfc",start_wfc,"start wave functions are from 'atomic' or 'file'");
	ModuleBase::GlobalFunc::OUTP(ofs,"start_charge",start_pot,"start charge is from 'atomic' or file");
	ModuleBase::GlobalFunc::OUTP(ofs,"charge_extrap",charge_extrap,"atomic; first-order; second-order; dm:coefficients of SIA");
	ModuleBase::GlobalFunc::OUTP(ofs,"out_charge",out_charge,">0 output charge density for selected electron steps");
	ModuleBase::GlobalFunc::OUTP(ofs,"out_potential",out_potential,"output realspace potential");
	ModuleBase::GlobalFunc::OUTP(ofs,"out_wf",out_wf,"output wave functions");
	ModuleBase::GlobalFunc::OUTP(ofs,"out_wf_r",out_wf_r,"output wave functions in realspace");
	ModuleBase::GlobalFunc::OUTP(ofs,"out_dos",out_dos,"output energy and dos");
	ModuleBase::GlobalFunc::OUTP(ofs,"out_band",out_band,"output energy and band structure");
	ModuleBase::GlobalFunc::OUTP(ofs,"restart_save",restart_save,"print to disk every step for restart");
	ModuleBase::GlobalFunc::OUTP(ofs,"restart_load",restart_load,"restart from disk");
	ModuleBase::GlobalFunc::OUTP(ofs,"read_file_dir",read_file_dir,"directory of files for reading");
	ModuleBase::GlobalFunc::OUTP(ofs,"nx",nx,"number of points along x axis for FFT grid");
	ModuleBase::GlobalFunc::OUTP(ofs,"ny",ny,"number of points along y axis for FFT grid");
	ModuleBase::GlobalFunc::OUTP(ofs,"nz",nz,"number of points along z axis for FFT grid");	
	ModuleBase::GlobalFunc::OUTP(ofs,"cell_factor",cell_factor,"used in the construction of the pseudopotential tables");	
	
	ofs << "\n#Parameters (3.Relaxation)" << std::endl;
	ModuleBase::GlobalFunc::OUTP(ofs,"ks_solver",GlobalV::KS_SOLVER,"cg; dav; lapack; genelpa; hpseps; scalapack_gvx");
	ModuleBase::GlobalFunc::OUTP(ofs,"niter",niter,"#number of electron iterations");
	ModuleBase::GlobalFunc::OUTP(ofs,"force_set",force_set,"output the force_set or not"); 
	ModuleBase::GlobalFunc::OUTP(ofs,"nstep",nstep,"number of ion iteration steps");
	ModuleBase::GlobalFunc::OUTP(ofs,"out_stru",out_stru,"output the structure files after each ion step");
	ModuleBase::GlobalFunc::OUTP(ofs,"force_thr",force_thr,"force threshold, unit: Ry/Bohr");
	ModuleBase::GlobalFunc::OUTP(ofs,"force_thr_ev",force_thr*13.6058/0.529177,"force threshold, unit: eV/Angstrom");
	ModuleBase::GlobalFunc::OUTP(ofs,"force_thr_ev2",force_thr_ev2,"force invalid threshold, unit: eV/Angstrom");
	ModuleBase::GlobalFunc::OUTP(ofs,"cg_threshold",cg_threshold,"threshold for switching from cg to bfgs, unit: eV/Angstrom");
	ModuleBase::GlobalFunc::OUTP(ofs,"stress_thr",stress_thr,"stress threshold");
	ModuleBase::GlobalFunc::OUTP(ofs,"press1",press1,"target pressure, unit: KBar");
	ModuleBase::GlobalFunc::OUTP(ofs,"press2",press2,"target pressure, unit: KBar");
	ModuleBase::GlobalFunc::OUTP(ofs,"press3",press3,"target pressure, unit: KBar");
	ModuleBase::GlobalFunc::OUTP(ofs,"bfgs_w1",bfgs_w1,"wolfe condition 1 for bfgs");
	ModuleBase::GlobalFunc::OUTP(ofs,"bfgs_w2",bfgs_w2,"wolfe condition 2 for bfgs");
	ModuleBase::GlobalFunc::OUTP(ofs,"trust_radius_max", trust_radius_max,"maximal trust radius, unit: Bohr");
	ModuleBase::GlobalFunc::OUTP(ofs,"trust_radius_min", trust_radius_min,"minimal trust radius, unit: Bohr");
	ModuleBase::GlobalFunc::OUTP(ofs,"trust_radius_ini", trust_radius_ini,"initial trust radius, unit: Bohr");
	ModuleBase::GlobalFunc::OUTP(ofs,"stress",stress,"calculate the stress or not");
	ModuleBase::GlobalFunc::OUTP(ofs,"fixed_axes",fixed_axes,"which axes are fixed");
	ModuleBase::GlobalFunc::OUTP(ofs,"move_method",ion_dynamics,"bfgs; sd; cg; cg_bfgs;"); //pengfei add 2013-08-15
	ModuleBase::GlobalFunc::OUTP(ofs,"out_level",out_level,"ie(for electrons); i(for ions);");
	ModuleBase::GlobalFunc::OUTP(ofs,"out_dm",out_dm,">0 output density matrix");

	// for deepks
	ModuleBase::GlobalFunc::OUTP(ofs,"out_descriptor",out_descriptor,">0 compute descriptor for deepks");
	ModuleBase::GlobalFunc::OUTP(ofs,"lmax_descriptor",lmax_descriptor,">0 lmax used in descriptor for deepks");

	ofs << "\n#Parameters (4.LCAO)" << std::endl;
	ModuleBase::GlobalFunc::OUTP(ofs,"basis_type",basis_type,"PW; LCAO in pw; LCAO");
	ModuleBase::GlobalFunc::OUTP(ofs,"new_dm",new_dm,"Type of density matrix; 0: old 1: new");
	if(ks_solver=="HPSEPS" || ks_solver=="genelpa" || ks_solver=="scalapack_gvx")
	{
		ModuleBase::GlobalFunc::OUTP(ofs,"nb2d",nb2d,"2d distribution of atoms");
	}
	ModuleBase::GlobalFunc::OUTP(ofs,"search_radius",search_radius,"input search radius (Bohr)");
	ModuleBase::GlobalFunc::OUTP(ofs,"search_pbc",search_pbc,"input periodic boundary condition");
	ModuleBase::GlobalFunc::OUTP(ofs,"lcao_ecut",lcao_ecut,"energy cutoff for LCAO");
	ModuleBase::GlobalFunc::OUTP(ofs,"lcao_dk",lcao_dk,"delta k for 1D integration in LCAO");
	ModuleBase::GlobalFunc::OUTP(ofs,"lcao_dr",lcao_dr,"delta r for 1D integration in LCAO");
	ModuleBase::GlobalFunc::OUTP(ofs,"lcao_rmax",lcao_rmax,"max R for 1D two-center integration table");
	ModuleBase::GlobalFunc::OUTP(ofs,"out_hs",out_hs,"output H and S matrix");
	ModuleBase::GlobalFunc::OUTP(ofs,"out_hs2",out_hs2,"output H(R) and S(R) matrix");
	ModuleBase::GlobalFunc::OUTP(ofs,"out_r",out_r_matrix,"output r(R) matrix");
	ModuleBase::GlobalFunc::OUTP(ofs,"out_lowf",out_lowf,"ouput LCAO wave functions");
	ModuleBase::GlobalFunc::OUTP(ofs,"bx",bx,"division of an element grid in FFT grid along x");
	ModuleBase::GlobalFunc::OUTP(ofs,"by",by,"division of an element grid in FFT grid along y");
	ModuleBase::GlobalFunc::OUTP(ofs,"bz",bz,"division of an element grid in FFT grid along z");

	ofs << "\n#Parameters (5.Smearing)" << std::endl;
	ModuleBase::GlobalFunc::OUTP(ofs,"smearing",smearing,"type of smearing: gauss; fd; fixed; mp; mp2; mv");
	ModuleBase::GlobalFunc::OUTP(ofs,"sigma",degauss,"energy range for smearing");
	
	ofs << "\n#Parameters (6.Charge Mixing)" << std::endl;
	ModuleBase::GlobalFunc::OUTP(ofs,"mixing_type",mixing_mode,"plain; kerker; pulay; pulay-kerker; broyden");
	ModuleBase::GlobalFunc::OUTP(ofs,"mixing_beta",mixing_beta,"mixing parameter: 0 means no new charge");
	ModuleBase::GlobalFunc::OUTP(ofs,"mixing_ndim",mixing_ndim,"mixing dimension in pulay");
	ModuleBase::GlobalFunc::OUTP(ofs,"mixing_gg0",mixing_gg0,"mixing parameter in kerker");

	ofs << "\n#Parameters (7.DOS)" << std::endl;
	ModuleBase::GlobalFunc::OUTP(ofs,"dos_emin_ev",dos_emin_ev,"minimal range for dos");
	ModuleBase::GlobalFunc::OUTP(ofs,"dos_emax_ev",dos_emax_ev,"maximal range for dos");
	ModuleBase::GlobalFunc::OUTP(ofs,"dos_edelta_ev",dos_edelta_ev,"delta energy for dos");
	ModuleBase::GlobalFunc::OUTP(ofs,"dos_scale",dos_scale,"scale dos range by");
	ModuleBase::GlobalFunc::OUTP(ofs,"dos_sigma",b_coef,"gauss b coefficeinet(default=0.07)");

	ofs << "\n#Parameters (8.Technique)" << std::endl;
	ModuleBase::GlobalFunc::OUTP(ofs,"gamma_only",gamma_only,"gamma only, only used in LCAO basis");
	ModuleBase::GlobalFunc::OUTP(ofs,"diago_proc",GlobalV::DIAGO_PROC,"number of proc used to diago");
	ModuleBase::GlobalFunc::OUTP(ofs,"npool",npool,"number of pools for k points, pw only");
	ModuleBase::GlobalFunc::OUTP(ofs,"mem_saver",mem_saver,"memory saver for many k points used");
	ModuleBase::GlobalFunc::OUTP(ofs,"printe",printe,"print band energy for selectively ionic steps");

	ofs << "\n#Parameters (9.SIAO)" << std::endl;
	ModuleBase::GlobalFunc::OUTP(ofs,"selinv_npole",selinv_npole,"number of selected poles");
	ModuleBase::GlobalFunc::OUTP(ofs,"selinv_temp",selinv_temp,"temperature for Fermi-Dirac distribution");
	ModuleBase::GlobalFunc::OUTP(ofs,"selinv_gap",selinv_gap,"supposed gap in the calculation");
	ModuleBase::GlobalFunc::OUTP(ofs,"selinv_deltae",selinv_deltae,"expected energy range");
	ModuleBase::GlobalFunc::OUTP(ofs,"selinv_mu",selinv_mu,"chosen mu as Fermi energy");
	ModuleBase::GlobalFunc::OUTP(ofs,"selinv_threshold",selinv_threshold,"threshold for calculated electron number");
	ModuleBase::GlobalFunc::OUTP(ofs,"selinv_niter",selinv_niter,"max number of steps to update mu");

	ofs << "\n#Parameters (10.Molecular dynamics)" << std::endl;
	ModuleBase::GlobalFunc::OUTP(ofs,"md_mdtype",mdp.mdtype,"choose ensemble");
	ModuleBase::GlobalFunc::OUTP(ofs,"md_potential",mdp.md_potential,"choose potential");
	ModuleBase::GlobalFunc::OUTP(ofs,"md_dt",mdp.dt,"time step");
	ModuleBase::GlobalFunc::OUTP(ofs,"mnhc",mdp.MNHC,"number of Nose-Hoover chains");
	ModuleBase::GlobalFunc::OUTP(ofs,"md_qmass",mdp.Qmass,"mass of thermostat");
	ModuleBase::GlobalFunc::OUTP(ofs,"md_tfirst",mdp.tfirst,"temperature first");
	ModuleBase::GlobalFunc::OUTP(ofs,"md_tlast",mdp.tlast,"temperature last");
	ModuleBase::GlobalFunc::OUTP(ofs,"md_dumpmdfred",mdp.recordFreq,"The period to dump MD information for monitoring and restarting MD");
	ModuleBase::GlobalFunc::OUTP(ofs,"md_mdoutpath",mdp.mdoutputpath,"output path of md");
	ModuleBase::GlobalFunc::OUTP(ofs,"md_rstmd",mdp.rstMD,"whether restart");
	ModuleBase::GlobalFunc::OUTP(ofs,"md_fixtemperature",mdp.fixTemperature,"period to change temperature");
	ModuleBase::GlobalFunc::OUTP(ofs,"md_ediff",mdp.ediff,"parameter for constraining total energy change");
	ModuleBase::GlobalFunc::OUTP(ofs,"md_ediffg",mdp.ediffg,"parameter for constraining max force change");
	ModuleBase::GlobalFunc::OUTP(ofs,"NVT_tau",mdp.NVT_tau,"parameter for adjust effect of thermostat");
	ModuleBase::GlobalFunc::OUTP(ofs,"NVT_control",mdp.NVT_control,"choose which thermostat used in NVT ensemble");
	ModuleBase::GlobalFunc::OUTP(ofs,"rcut_lj",mdp.rcut_lj,"cutoff radius of LJ potential");
	ModuleBase::GlobalFunc::OUTP(ofs,"epsilon_lj",mdp.epsilon_lj,"the value of epsilon for LJ potential");
	ModuleBase::GlobalFunc::OUTP(ofs,"sigma_lj",mdp.sigma_lj,"the value of sigma for LJ potential");

	ofs << "\n#Parameters (11.Efield)" << std::endl;
	ModuleBase::GlobalFunc::OUTP(ofs,"efield",efield,"add electric field");
	ModuleBase::GlobalFunc::OUTP(ofs,"edir",edir,"add electric field");
	ModuleBase::GlobalFunc::OUTP(ofs,"emaxpos",emaxpos,"maximal position of efield [0,1)");
	ModuleBase::GlobalFunc::OUTP(ofs,"eopreg",eopreg,"where sawlike potential decrease");
	ModuleBase::GlobalFunc::OUTP(ofs,"eamp",eamp,"amplitute of the efield, unit is a.u.");
	ModuleBase::GlobalFunc::OUTP(ofs,"eamp_v",eamp*51.44,"amplitute of the efield, unit is V/A");

	ofs << "\n#Parameters (12.Test)" << std::endl;
	ModuleBase::GlobalFunc::OUTP(ofs,"out_alllog",out_alllog,"output information for each processor, when parallel");
	ModuleBase::GlobalFunc::OUTP(ofs,"nurse", nurse,"for coders");
	ModuleBase::GlobalFunc::OUTP(ofs,"colour", colour,"for coders, make their live colourful");
	ModuleBase::GlobalFunc::OUTP(ofs,"t_in_h", t_in_h,"calculate the kinetic energy or not");
	ModuleBase::GlobalFunc::OUTP(ofs,"vl_in_h", vl_in_h,"calculate the local potential or not");
	ModuleBase::GlobalFunc::OUTP(ofs,"vnl_in_h", vnl_in_h,"calculate the nonlocal potential or not");
	ModuleBase::GlobalFunc::OUTP(ofs,"vh_in_h", vh_in_h,"calculate the hartree potential or not");
	ModuleBase::GlobalFunc::OUTP(ofs,"vxc_in_h", vxc_in_h,"calculate the xc potential or not");
	ModuleBase::GlobalFunc::OUTP(ofs,"vion_in_h", vion_in_h,"calculate the local ionic potential or not");
	ModuleBase::GlobalFunc::OUTP(ofs,"test_force", test_force, "test the force");
	ModuleBase::GlobalFunc::OUTP(ofs,"test_stress", test_stress, "test the force");
	
	ofs << "\n#Parameters (13.Other Methods)" << std::endl;
	ModuleBase::GlobalFunc::OUTP(ofs,"mlwf_flag",mlwf_flag,"turn MLWF on or off");
	ModuleBase::GlobalFunc::OUTP(ofs,"opt_epsilon2",opt_epsilon2,"calculate the dielectic function");
	ModuleBase::GlobalFunc::OUTP(ofs,"opt_nbands",opt_nbands,"number of bands for optical calculation");
	
	ofs << "\n#Parameters (14.VdW Correction)" << std::endl;								
	ModuleBase::GlobalFunc::OUTP(ofs,"vdw_method",vdw_method,"the method of calculating vdw (none ; d2 ; d3_0 ; d3_bj");
	ModuleBase::GlobalFunc::OUTP(ofs,"vdw_s6",vdw_s6,"scale parameter of d2/d3_0/d3_bj");
    ModuleBase::GlobalFunc::OUTP(ofs,"vdw_s8",vdw_s8,"scale parameter of d3_0/d3_bj");
    ModuleBase::GlobalFunc::OUTP(ofs,"vdw_a1",vdw_a1,"damping parameter of d3_0/d3_bj");
    ModuleBase::GlobalFunc::OUTP(ofs,"vdw_a2",vdw_a2,"damping parameter of d3_bj");
    ModuleBase::GlobalFunc::OUTP(ofs,"vdw_d",vdw_d,"damping parameter of d2");		
    ModuleBase::GlobalFunc::OUTP(ofs,"vdw_abc",vdw_abc,"third-order term?");
	ModuleBase::GlobalFunc::OUTP(ofs,"vdw_C6_file",vdw_C6_file,"filename of C6");
    ModuleBase::GlobalFunc::OUTP(ofs,"vdw_C6_unit",vdw_C6_unit,"unit of C6, Jnm6/mol or eVA6");
    ModuleBase::GlobalFunc::OUTP(ofs,"vdw_R0_file",vdw_R0_file,"filename of R0");
    ModuleBase::GlobalFunc::OUTP(ofs,"vdw_R0_unit",vdw_R0_unit,"unit of R0, A or Bohr");
	ModuleBase::GlobalFunc::OUTP(ofs,"vdw_model",vdw_model,"expression model of periodic structure, radius or period");
    ModuleBase::GlobalFunc::OUTP(ofs,"vdw_radius",vdw_radius,"radius cutoff for periodic structure");
    ModuleBase::GlobalFunc::OUTP(ofs,"vdw_radius_unit",vdw_radius_unit,"unit of radius cutoff for periodic structure");
    ModuleBase::GlobalFunc::OUTP(ofs,"vdw_cn_thr",vdw_cn_thr,"radius cutoff for cn");
    ModuleBase::GlobalFunc::OUTP(ofs,"vdw_cn_thr_unit",vdw_cn_thr_unit,"unit of cn_thr, Bohr or Angstrom");
	ofs << std::setw(20) << "vdw_period" << vdw_period.x 
			<< " " << vdw_period.y << " " 
			<< vdw_period.z<< " #periods of periodic structure" << std::endl;
	
	
	ofs << "\n#Parameters (15.spectrum)" << std::endl;              // pengfei Li add 2016-11-23
	ModuleBase::GlobalFunc::OUTP(ofs,"spectral_type",spectral_type,"the type of the calculated spectrum");
	ModuleBase::GlobalFunc::OUTP(ofs,"spectral_method",spectral_method,"0: tddft(linear response)");
	ModuleBase::GlobalFunc::OUTP(ofs,"kernel_type",kernel_type,"the kernel type: rpa, tdlda ...");
	ModuleBase::GlobalFunc::OUTP(ofs,"eels_method",eels_method,"0: hilbert_transform method; 1: standard method");
	ModuleBase::GlobalFunc::OUTP(ofs,"absorption_method",absorption_method,"0: vasp's method  1: pwscf's method");
	ModuleBase::GlobalFunc::OUTP(ofs,"system",system_type,"the calculate system");
	ModuleBase::GlobalFunc::OUTP(ofs,"eta",eta,"eta(Ry)");
	ModuleBase::GlobalFunc::OUTP(ofs,"domega",domega,"domega(Ry)");
	ModuleBase::GlobalFunc::OUTP(ofs,"nomega",nomega,"nomega");
	ModuleBase::GlobalFunc::OUTP(ofs,"ecut_chi",ecut_chi,"the dimension of chi matrix");
	ofs << std::setw(20) <<"q_start"<<q_start[0]<<"   "
		<<q_start[1]<<"   "<<q_start[2]
		<<"  #the position of the first q point in direct coordinate" <<std::endl;
	ofs << std::setw(20) <<"q_direction"<<q_direct[0]<<"   "<<q_direct[1]<<"   "<<q_direct[2]<<"  #the q direction" <<std::endl;
	ModuleBase::GlobalFunc::OUTP(ofs,"nq",nq,"the total number of qpoints for calculation");
	ModuleBase::GlobalFunc::OUTP(ofs,"out_epsilon",out_epsilon,"output epsilon or not");
	ModuleBase::GlobalFunc::OUTP(ofs,"out_chi",out_chi,"output chi or not");
	ModuleBase::GlobalFunc::OUTP(ofs,"out_chi0",out_chi0,"output chi0 or not");
	ModuleBase::GlobalFunc::OUTP(ofs,"fermi_level",fermi_level,"the change of the fermi_level(Ry)");
	ModuleBase::GlobalFunc::OUTP(ofs,"coulomb_cutoff",coulomb_cutoff," turn on the coulomb_cutoff or not");
	ModuleBase::GlobalFunc::OUTP(ofs,"kmesh_interpolation",kmesh_interpolation,"calculting <i,0|j,R>");
	for(int i=0; i<nq; i++)
	{
		ofs << std::setw(20) <<"qcar" << qcar[i][0] <<"   "<< qcar[i][1] <<"   "<<qcar[i][2]<<"  #(unit: 2PI/lat0)" << std::endl;
	}
	ModuleBase::GlobalFunc::OUTP(ofs,"ocp",GlobalV::ocp,"change occupation or not");
	ModuleBase::GlobalFunc::OUTP(ofs,"ocp_set",GlobalV::ocp_set,"set occupation");
	//ModuleBase::GlobalFunc::OUTP(ofs,"ocp_n",ocp_n,"number of occupation");
	// for(int i=0; i<ocp_n; i++)
	// {
		// ofs << std::setw(20) <<"ocp_kb" << GlobalV::ocp_kb[i]<< std::endl;
	// }
	ofs << std::setw(20) <<"lcao_box"<<lcao_box[0]
		<<"   "<<lcao_box[1]<<"   "
		<<lcao_box[2]<<"  #the scale for searching the existence of the overlap <i,0|j,R>" <<std::endl;
	ModuleBase::GlobalFunc::OUTP(ofs," mulliken", GlobalV::mulliken," mulliken  charge or not");//qifeng add 2019/9/10
	
	//ModuleBase::GlobalFunc::OUTP(ofs,"epsilon0",epsilon0,"calculate the macroscopic dielectric constant or not");
	ModuleBase::GlobalFunc::OUTP(ofs,"intrasmear",intrasmear,"Eta");
	ModuleBase::GlobalFunc::OUTP(ofs,"shift",shift,"shift");
	ModuleBase::GlobalFunc::OUTP(ofs,"metalcalc",metalcalc,"metal or not");
	ModuleBase::GlobalFunc::OUTP(ofs,"eps_degauss",eps_degauss,"degauss in calculating epsilon0");
	ModuleBase::GlobalFunc::OUTP(ofs,"noncolin",noncolin,"using non-collinear-spin");
	ModuleBase::GlobalFunc::OUTP(ofs,"lspinorb",lspinorb,"consider the spin-orbit interaction");
	
	//ModuleBase::GlobalFunc::OUTP(ofs,"epsilon0_choice",epsilon0_choice,"0: vasp's method  1: pwscf's method");
	
	ofs << "\n#Parameters (17.exx)" << std::endl;
	ModuleBase::GlobalFunc::OUTP(ofs,"exx_hybrid_type",exx_hybrid_type,"no, hf, pbe0, hse or opt_orb");
	ModuleBase::GlobalFunc::OUTP(ofs,"exx_hybrid_alpha",exx_hybrid_alpha,"");
	ModuleBase::GlobalFunc::OUTP(ofs,"exx_hse_omega",exx_hse_omega,"");
	ModuleBase::GlobalFunc::OUTP(ofs,"exx_separate_loop",exx_separate_loop,"0 or 1");
	ModuleBase::GlobalFunc::OUTP(ofs,"exx_hybrid_step",exx_hybrid_step,"");
	ModuleBase::GlobalFunc::OUTP(ofs,"exx_lambda",exx_lambda,"");
	ModuleBase::GlobalFunc::OUTP(ofs,"exx_pca_threshold",exx_pca_threshold,"");
	ModuleBase::GlobalFunc::OUTP(ofs,"exx_c_threshold",exx_c_threshold,"");
	ModuleBase::GlobalFunc::OUTP(ofs,"exx_v_threshold",exx_v_threshold,"");
	ModuleBase::GlobalFunc::OUTP(ofs,"exx_dm_threshold",exx_dm_threshold,"");
	ModuleBase::GlobalFunc::OUTP(ofs,"exx_schwarz_threshold",exx_schwarz_threshold,"");
	ModuleBase::GlobalFunc::OUTP(ofs,"exx_cauchy_threshold",exx_cauchy_threshold,"");
	ModuleBase::GlobalFunc::OUTP(ofs,"exx_ccp_threshold",exx_ccp_threshold,"");
	ModuleBase::GlobalFunc::OUTP(ofs,"exx_ccp_rmesh_times",exx_ccp_rmesh_times,"");
	ModuleBase::GlobalFunc::OUTP(ofs,"exx_distribute_type",exx_distribute_type,"htime or kmeans1 or kmeans2");
	ModuleBase::GlobalFunc::OUTP(ofs,"exx_opt_orb_lmax",exx_opt_orb_lmax,"");
	ModuleBase::GlobalFunc::OUTP(ofs,"exx_opt_orb_ecut",exx_opt_orb_ecut,"");
	ModuleBase::GlobalFunc::OUTP(ofs,"exx_opt_orb_tolerence",exx_opt_orb_tolerence,"");

	ofs << "\n#Parameters (17.tddft)" << std::endl;
	ModuleBase::GlobalFunc::OUTP(ofs,"tddft",tddft,"calculate tddft or not");
	ModuleBase::GlobalFunc::OUTP(ofs,"td_dr2",td_dr2,"threshold for electronic iteration of tddft");
	ModuleBase::GlobalFunc::OUTP(ofs,"td_dt",td_dt,"time of ion step");
	ModuleBase::GlobalFunc::OUTP(ofs,"td_force_dt",td_force_dt,"time of force change");
	ModuleBase::GlobalFunc::OUTP(ofs,"td_val_elec_01",td_val_elec_01,"td_val_elec_01");
	ModuleBase::GlobalFunc::OUTP(ofs,"td_val_elec_02",td_val_elec_02,"td_val_elec_02");
	ModuleBase::GlobalFunc::OUTP(ofs,"td_val_elec_03",td_val_elec_03,"td_val_elec_03");
	ModuleBase::GlobalFunc::OUTP(ofs,"td_vext",td_vext,"add extern potential or not");
	ModuleBase::GlobalFunc::OUTP(ofs,"td_vext_dire",td_vext_dire,"extern potential direction");
	ModuleBase::GlobalFunc::OUTP(ofs,"td_timescale",td_timescale,"extern potential td_timescale");
	ModuleBase::GlobalFunc::OUTP(ofs,"td_vexttype",td_vexttype,"extern potential type");
	ModuleBase::GlobalFunc::OUTP(ofs,"td_vextout",td_vextout,"output extern potential or not");
	ModuleBase::GlobalFunc::OUTP(ofs,"td_dipoleout",td_dipoleout,"output dipole or not");

	
	ofs << "\n#Parameters (18.berry_wannier)" << std::endl;
	ModuleBase::GlobalFunc::OUTP(ofs,"berry_phase",berry_phase,"calculate berry phase or not");
	ModuleBase::GlobalFunc::OUTP(ofs,"gdir",gdir,"calculate the polarization in the direction of the lattice std::vector");
	ModuleBase::GlobalFunc::OUTP(ofs,"towannier90",towannier90,"use wannier90 code interface or not");
	ModuleBase::GlobalFunc::OUTP(ofs,"nnkpfile",NNKP,"the wannier90 code nnkp file name");
	ModuleBase::GlobalFunc::OUTP(ofs,"wannier_spin",wannier_spin,"calculate spin in wannier90 code interface");
	
    ofs.close();
    return;
}
