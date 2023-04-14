#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "module_io/input_conv.h"
#include "module_base/global_variable.h"
#include "for_testing_input_conv.h"

/************************************************
 *  unit test of input_conv.cpp
 ***********************************************/

/**
 * - Tested Functions:
 *   - Convert()
 */

#define private public
#include "module_io/input.h"

class InputConvTest : public testing::Test
{
protected:
	std::string output;
};

TEST_F(InputConvTest, Conv)
{
	INPUT.Default();
	std::string input_file = "./support/INPUT";
	INPUT.Read(input_file);
	Input_Conv::Convert();
	EXPECT_EQ(GlobalV::stru_file,INPUT.stru_file);
	EXPECT_EQ(GlobalV::global_wannier_card,"none");
	EXPECT_EQ(GlobalV::global_kpoint_card,INPUT.kpoint_file);
	EXPECT_EQ(GlobalV::global_pseudo_dir,INPUT.pseudo_dir+"/");
	EXPECT_EQ(GlobalV::global_orbital_dir,INPUT.orbital_dir+"/");
	EXPECT_EQ(GlobalC::ucell.latName,"none");
	EXPECT_EQ(GlobalC::ucell.ntype,1);
	EXPECT_EQ(GlobalC::ucell.init_vel,false);
	EXPECT_EQ(GlobalC::ucell.lc[0],1);
	EXPECT_EQ(GlobalC::ucell.lc[1],1);
	EXPECT_EQ(GlobalC::ucell.lc[2],1);
	EXPECT_EQ(GlobalV::KSPACING,0.0);
	EXPECT_DOUBLE_EQ(GlobalV::MIN_DIST_COEF,0.2);
	EXPECT_EQ(GlobalV::NBANDS,8);
	EXPECT_EQ(GlobalC::wf.pw_seed,1);
	EXPECT_EQ(GlobalV::NBANDS_ISTATE,5);
	EXPECT_EQ(GlobalV::device_flag,"cpu");
	EXPECT_EQ(GlobalV::KPAR,1);
	EXPECT_EQ(GlobalV::NSTOGROUP,1);
	EXPECT_EQ(GlobalV::precision_flag,"double");
	EXPECT_EQ(GlobalV::CALCULATION,"scf");
	EXPECT_EQ(GlobalV::ESOLVER_TYPE,"ksdft");
	EXPECT_DOUBLE_EQ(GlobalV::PSEUDORCUT,15.0);
	EXPECT_EQ(GlobalV::PSEUDO_MESH,false);
	EXPECT_EQ(GlobalV::DFT_FUNCTIONAL,"hse");
	EXPECT_DOUBLE_EQ(GlobalV::XC_TEMPERATURE,0.0);
	EXPECT_EQ(GlobalV::NSPIN,1);
	EXPECT_EQ(GlobalV::CURRENT_SPIN,0);
	EXPECT_EQ(GlobalV::CAL_FORCE,0);
	EXPECT_NEAR(GlobalV::FORCE_THR,0.001,0.0000001);
	EXPECT_DOUBLE_EQ(GlobalV::STRESS_THR,0.01);
	EXPECT_EQ(GlobalV::PRESS1,0);
	EXPECT_EQ(GlobalV::PRESS2,0);
	EXPECT_EQ(GlobalV::PRESS3,0);
	EXPECT_EQ(GlobalV::out_element_info,0);
	EXPECT_EQ(Force_Stress_LCAO::force_invalid_threshold_ev,0);
	EXPECT_DOUBLE_EQ(BFGS_Basic::relax_bfgs_w1,0.01);
	EXPECT_DOUBLE_EQ(BFGS_Basic::relax_bfgs_w2,0.5);
	EXPECT_DOUBLE_EQ(Ions_Move_Basic::relax_bfgs_rmax,0.8);
	EXPECT_DOUBLE_EQ(Ions_Move_Basic::relax_bfgs_rmin,0.00001);
	EXPECT_EQ(Ions_Move_Basic::relax_bfgs_init,0.5);
	EXPECT_EQ(Ions_Move_Basic::out_stru,0);
	EXPECT_EQ(Lattice_Change_Basic::out_stru,0);
	EXPECT_EQ(GlobalV::CAL_STRESS,0);
	EXPECT_EQ(GlobalV::RELAX_METHOD,"cg");
	EXPECT_DOUBLE_EQ(GlobalV::relax_scale_force,0.5);
	EXPECT_EQ(GlobalV::relax_new,true);
	EXPECT_EQ(GlobalV::OUT_LEVEL,"ie");
	EXPECT_DOUBLE_EQ(Ions_Move_CG::RELAX_CG_THR,0.5);
	EXPECT_EQ(ModuleSymmetry::Symmetry::symm_flag,1);
	EXPECT_DOUBLE_EQ(GlobalC::symm.epsilon,0.00001);
	EXPECT_EQ(GlobalV::BASIS_TYPE,"lcao");
	EXPECT_EQ(GlobalV::KS_SOLVER,"genelpa");
	EXPECT_EQ(GlobalV::SEARCH_RADIUS,-1);
	EXPECT_EQ(GlobalV::SEARCH_PBC,1);
	EXPECT_EQ(GlobalV::GAMMA_ONLY_LOCAL,true);
	EXPECT_EQ(GlobalV::DIAGO_PROC,4);
	EXPECT_EQ(GlobalV::PW_DIAG_NMAX,50);
	EXPECT_EQ(GlobalV::DIAGO_CG_PREC,1);
	EXPECT_EQ(GlobalV::PW_DIAG_NDIM,4);
	EXPECT_DOUBLE_EQ(GlobalV::PW_DIAG_THR,0.01);
	EXPECT_EQ(GlobalV::NB2D,0);
	EXPECT_EQ(GlobalV::NURSE,0);
	EXPECT_EQ(GlobalV::COLOUR,0);
	EXPECT_EQ(GlobalV::T_IN_H,1);
	EXPECT_EQ(GlobalV::VL_IN_H,1);
	EXPECT_EQ(GlobalV::VNL_IN_H,1);
	EXPECT_EQ(GlobalV::VH_IN_H,1);
	EXPECT_EQ(GlobalV::VION_IN_H,1);
	EXPECT_EQ(GlobalV::TEST_FORCE,0);
	EXPECT_EQ(GlobalV::TEST_STRESS,0);
	EXPECT_EQ(GlobalV::test_skip_ewald,0);
	EXPECT_DOUBLE_EQ(GlobalV::SCF_THR,0.00000001);
	EXPECT_EQ(GlobalC::wf.init_wfc,"atomic");
	EXPECT_EQ(GlobalC::wf.mem_saver,0);
	EXPECT_EQ(GlobalC::en.printe,100);
	EXPECT_EQ(GlobalV::LSPINORB,false);
	EXPECT_EQ(GlobalV::NONCOLIN,false);
	EXPECT_EQ(GlobalV::DOMAG,false);
	EXPECT_EQ(GlobalV::DOMAG_Z,false);
	EXPECT_EQ(GlobalV::NPOL,1);
	EXPECT_EQ(GlobalV::EFIELD_FLAG,0);
	EXPECT_EQ(GlobalV::DIP_COR_FLAG,0);
	EXPECT_EQ(elecstate::Efield::efield_dir,2);
	EXPECT_DOUBLE_EQ(elecstate::Efield::efield_pos_max,0.5);
	EXPECT_DOUBLE_EQ(elecstate::Efield::efield_pos_dec,0.1);
	EXPECT_EQ(elecstate::Efield::efield_amp,0);
	EXPECT_EQ(GlobalV::GATE_FLAG,0);
	EXPECT_EQ(GlobalV::nelec,0);
	EXPECT_DOUBLE_EQ(elecstate::Gatefield::zgate,0.5);
	EXPECT_EQ(elecstate::Gatefield::relax,0);
	EXPECT_EQ(elecstate::Gatefield::block,0);
	EXPECT_DOUBLE_EQ(elecstate::Gatefield::block_down,0.45);
	EXPECT_DOUBLE_EQ(elecstate::Gatefield::block_up,0.55);
	EXPECT_DOUBLE_EQ(elecstate::Gatefield::block_height,0.1);

	EXPECT_EQ(ELEC_evolve::td_force_dt,0.02);
	EXPECT_EQ(ELEC_evolve::td_val_elec_01,1);
	EXPECT_EQ(ELEC_evolve::td_val_elec_02,1);
	EXPECT_EQ(ELEC_evolve::td_val_elec_03,1);
	EXPECT_EQ(ELEC_evolve::td_vext,false);
	EXPECT_EQ(ELEC_evolve::out_dipole,false);
	EXPECT_EQ(ELEC_evolve::out_efield,false);
	EXPECT_EQ(ELEC_evolve::td_print_eij,-1.0);
	EXPECT_EQ(ELEC_evolve::td_edm,0);
	EXPECT_EQ(GlobalV::ocp,false);
	EXPECT_EQ(GlobalV::ocp_set,INPUT.ocp_set);
	EXPECT_EQ(GlobalV::out_mul,0);
	EXPECT_EQ(GlobalC::ppcell.cell_factor,1.2);
	EXPECT_EQ(GlobalV::SCF_NMAX,50);
	EXPECT_EQ(GlobalV::RELAX_NMAX,1);
	EXPECT_EQ(GlobalV::MD_NSTEP,10);
	EXPECT_EQ(GlobalV::OUT_FREQ_ELEC,0);
	EXPECT_EQ(GlobalV::OUT_FREQ_ION,0);
	EXPECT_EQ(GlobalV::init_chg,"atomic");
	EXPECT_EQ(GlobalV::chg_extrap,"atomic");
	EXPECT_EQ(GlobalV::out_chg,false);
	EXPECT_EQ(GlobalV::nelec,0.0);
	EXPECT_EQ(GlobalV::out_pot,2);
	EXPECT_EQ(GlobalC::wf.out_wfc_pw,0);
	EXPECT_EQ(GlobalC::wf.out_wfc_r,false);
	EXPECT_EQ(GlobalC::en.out_dos,0);
	EXPECT_EQ(GlobalC::en.out_band,false);
	EXPECT_EQ(GlobalC::en.out_proj_band,false);
	EXPECT_EQ(GlobalV::out_app_flag,false);
	EXPECT_EQ(GlobalV::out_bandgap,false);
	EXPECT_EQ(Local_Orbital_Charge::out_dm,false);
	EXPECT_EQ(Local_Orbital_Charge::out_dm1,false);
	EXPECT_EQ(hsolver::HSolverLCAO::out_mat_hs,false);
	EXPECT_EQ( hsolver::HSolverLCAO::out_mat_hsR,false);
	EXPECT_EQ(hsolver::HSolverLCAO::out_mat_t,false);
	EXPECT_EQ(hsolver::HSolverLCAO::out_mat_dh,INPUT.out_mat_dh);
	EXPECT_EQ(hsolver::HSolverLCAO::out_hsR_interval,1);
	EXPECT_EQ(elecstate::ElecStateLCAO::out_wfc_lcao,false);
	EXPECT_EQ(GlobalC::en.dos_emin_ev,-15);
	EXPECT_EQ(GlobalC::en.dos_emax_ev,15);
	EXPECT_EQ(GlobalC::en.dos_edelta_ev,0.01);
	EXPECT_EQ(GlobalC::en.dos_scale,0.01);
	EXPECT_EQ(GlobalC::en.bcoeff,0.07);
	EXPECT_EQ(berryphase::berry_phase_flag,false);
	EXPECT_EQ(GlobalV::imp_sol,false);
	EXPECT_EQ(GlobalV::eb_k,80);
	EXPECT_EQ(GlobalV::tau,INPUT.tau);
	EXPECT_EQ(GlobalV::sigma_k,0.6);
	EXPECT_EQ(GlobalV::sigma_k,0.6);
	EXPECT_EQ( GlobalV::nc_k,0.00037);
	EXPECT_EQ( GlobalV::of_kinetic,"vw");
	EXPECT_EQ(GlobalV::of_method,"tn");
	EXPECT_EQ(GlobalV::of_conv,"energy");
	EXPECT_EQ(GlobalV::of_tole,0.000001);
	EXPECT_EQ(GlobalV::of_tolp,0.00001);
	EXPECT_EQ(GlobalV::of_tf_weight,1.0);
	EXPECT_EQ(GlobalV::of_vw_weight,1.0);
	EXPECT_EQ(GlobalV::of_wt_alpha,0.833333);
	EXPECT_EQ(GlobalV::of_wt_beta,0.833333);
	EXPECT_EQ(GlobalV::of_wt_rho0,1.0);
	EXPECT_EQ(GlobalV::of_hold_rho0,false);
    EXPECT_EQ(GlobalV::of_full_pw,false);
	EXPECT_EQ(GlobalV::of_full_pw_dim,0);
	EXPECT_EQ(GlobalV::of_read_kernel,false);
	EXPECT_EQ(GlobalV::of_kernel_file,"WTkernel.txt");
	EXPECT_EQ(GlobalV::global_readin_dir,GlobalV::global_out_dir);
}

TEST_F(InputConvTest, ConvRelax)
{
	INPUT.Default();
	std::string input_file = "./support/INPUT";
	INPUT.Read(input_file);
	INPUT.calculation="relax";
	INPUT.fixed_ibrav=1;
	INPUT.relax_new=false;
	std::string output2;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(Input_Conv::Convert(), ::testing::ExitedWithCode(0),"");
	output2 = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output2,testing::HasSubstr(" \n !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
		"                         NOTICE                           \n !"
		"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n \n fixed_ibrav only available for relax_new = 1\n"
		" CHECK IN FILE : warning.log\n \n !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
		"                         NOTICE                           \n"
		" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n  |CLASS_NAME---------|NAME---------------|TIME(Sec)-----|CALLS----|AVG------|PER%-------\n -------"
		"---------------------------------------------------------------------------------\n See output information in : \n"
		));

	INPUT.Read(input_file);
	INPUT.calculation="relax";
	INPUT.latname="none";
	INPUT.fixed_ibrav=1;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(Input_Conv::Convert(), ::testing::ExitedWithCode(0),"");
	output2 = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output2,testing::HasSubstr(" \n !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
		"                         NOTICE                           \n !"
		"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n \n to use fixed_ibrav, latname must be provided\n"
		" CHECK IN FILE : warning.log\n \n !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
		"                         NOTICE                           \n"
		" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n  |CLASS_NAME---------|NAME---------------|TIME(Sec)-----|CALLS----|AVG------|PER%-------\n -------"
		"---------------------------------------------------------------------------------\n See output information in : \n"
		));

	INPUT.Read(input_file);
	INPUT.calculation="relax";
	INPUT.fixed_atoms=1;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(Input_Conv::Convert(), ::testing::ExitedWithCode(0),"");
	output2 = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output2,testing::HasSubstr(" \n !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
		"                         NOTICE                           \n !"
		"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n \n fixed_atoms is not meant to be used for calculation = relax\n"
		" CHECK IN FILE : warning.log\n \n !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
		"                         NOTICE                           \n"
		" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n  |CLASS_NAME---------|NAME---------------|TIME(Sec)-----|CALLS----|AVG------|PER%-------\n -------"
		"---------------------------------------------------------------------------------\n See output information in : \n"
		));

	INPUT.Read(input_file);
	INPUT.calculation="relax";
	INPUT.relax_new=false;
	INPUT.fixed_axes="shape";
	testing::internal::CaptureStdout();
	EXPECT_EXIT(Input_Conv::Convert(), ::testing::ExitedWithCode(0),"");
	output2 = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output2,testing::HasSubstr(" \n !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
		"                         NOTICE                           \n !"
		"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n \n fixed shape and fixed volume only supported for relax_new = 1\n"
		" CHECK IN FILE : warning.log\n \n !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
		"                         NOTICE                           \n"
		" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n  |CLASS_NAME---------|NAME---------------|TIME(Sec)-----|CALLS----|AVG------|PER%-------\n -------"
		"---------------------------------------------------------------------------------\n See output information in : \n"
		));

	INPUT.Default();
	INPUT.Read(input_file);
	INPUT.calculation="relax";
	INPUT.relax_new=true;
	INPUT.relax_method="sd";
	Input_Conv::Convert();
	EXPECT_EQ(INPUT.relax_new,false);
}

    TEST_F(InputConvTest, ConvGpu)
    {
    INPUT.Default();
   std::string input_file = "./support/INPUT";
    INPUT.Read(input_file);
    INPUT.device="gpu";
    INPUT.ks_solver="cg";//lapack
    INPUT.basis_type="pw";//LCAO
	std::string output2;
		testing::internal::CaptureStdout();
		EXPECT_EXIT(Input_Conv::Convert(), ::testing::ExitedWithCode(0),"");
		output2 = testing::internal::GetCapturedStdout();
		EXPECT_THAT(output2,testing::HasSubstr("INPUT device setting does not match the request!"
			"\n Input device = gpu"
			"\n Input basis_type = pw"
			"\n Input ks_solver = cg"
			"\n Compile setting = host"
			"\n Environment device_num = -1"
			"\n"));
  }

TEST_F(InputConvTest, ConvRelax2)
{
	INPUT.Default();
	std::string input_file = "./support/INPUT";
	INPUT.Read(input_file);
	INPUT.calculation="relax";
	INPUT.chg_extrap="first-order";
	std::string output2;
	testing::internal::CaptureStdout();
	Input_Conv::Convert();
	output2 = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output2,testing::HasSubstr(" For relaxation, charge extrapolation is set to atomic.\n"));
	EXPECT_EQ(INPUT.chg_extrap,"atomic");
}

TEST_F(InputConvTest, dftplus)
{
	INPUT.Default();
	std::string input_file = "./support/INPUT";
	INPUT.Read(input_file);
	INPUT.dft_plus_u=true;
	Input_Conv::Convert();
	EXPECT_EQ(GlobalV::dft_plus_u,true);
	EXPECT_EQ(GlobalC::dftu.Yukawa,false);
	EXPECT_EQ(GlobalC::dftu.omc,false);//
	EXPECT_EQ(GlobalC::dftu.orbital_corr,INPUT.orbital_corr);
	EXPECT_EQ(GlobalC::dftu.mixing_dftu,false);
	EXPECT_EQ(GlobalC::dftu.U,INPUT.hubbard_u);
}

TEST_F(InputConvTest, nspin)
{
	INPUT.Default();
	std::string input_file = "./support/INPUT";
	INPUT.Read(input_file);
	INPUT.noncolin=true;
	Input_Conv::Convert();
	EXPECT_EQ(GlobalV::NSPIN,4);
	EXPECT_EQ(GlobalV::NONCOLIN,true);
	EXPECT_EQ(GlobalV::NPOL,2);
	EXPECT_EQ(GlobalV::DOMAG,false);
	EXPECT_EQ(GlobalV::DOMAG_Z,true);
	EXPECT_EQ(GlobalV::LSPINORB,false);
	EXPECT_EQ(GlobalV::soc_lambda,INPUT.soc_lambda);
}

TEST_F(InputConvTest, nspinbeta)
{
	INPUT.Default();
	std::string input_file = "./support/INPUT";
	INPUT.Read(input_file);
	INPUT.noncolin=true;
	INPUT.cal_stress=true;
	std::string output2;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(Input_Conv::Convert(), ::testing::ExitedWithCode(0),"");
	output2 = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output2,testing::HasSubstr(" \n !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
		"                         NOTICE                           \n !"
		"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n \n force & stress not ready for soc yet!\n"
		" CHECK IN FILE : warning.log\n \n !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
		"                         NOTICE                           \n"
		" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n  |CLASS_NAME---------|NAME---------------|TIME(Sec)-----|CALLS----|AVG------|PER%-------\n -------"
		"---------------------------------------------------------------------------------\n See output information in : \n"
		));
}

TEST_F(InputConvTest, nupdown)
{
	INPUT.Default();
	std::string input_file = "./support/INPUT";
	INPUT.Read(input_file);
	INPUT.nupdown=0.1;
	Input_Conv::Convert();
	EXPECT_EQ(GlobalV::TWO_EFERMI,true);
	EXPECT_EQ(GlobalV::nupdown,0.1);
}

TEST_F(InputConvTest,restart_save )
{
	INPUT.Default();
	std::string input_file = "./support/INPUT";
	INPUT.Read(input_file);
	INPUT.restart_save=true;
	GlobalV::MY_RANK=0;
	INPUT.dft_functional="hf";
	Input_Conv::Convert();
	EXPECT_EQ(GlobalC::restart.folder,GlobalV::global_readin_dir + "restart/");
	int a=access(GlobalC::restart.folder.c_str(),0);
	EXPECT_EQ(a,0);
	EXPECT_EQ(GlobalC::restart.info_save.save_charge,true);
	EXPECT_EQ(GlobalC::restart.info_save.save_H,true);
}

TEST_F(InputConvTest,restart_save2 )
{
	INPUT.Default();
	std::string input_file = "./support/INPUT";
	INPUT.Read(input_file);
	INPUT.restart_save=true;
	INPUT.dft_functional = "default";
	Input_Conv::Convert();
	EXPECT_EQ(GlobalC::restart.info_save.save_charge,true);
}

TEST_F(InputConvTest, restart_load)
{
	INPUT.Default();
	std::string input_file = "./support/INPUT";
	INPUT.Read(input_file);
	INPUT.restart_load=true;
	INPUT.dft_functional = "hf";
	Input_Conv::Convert();
	EXPECT_EQ( GlobalC::restart.folder,GlobalV::global_readin_dir + "restart/");
	EXPECT_EQ(GlobalC::restart.info_load.load_charge,true);
}

TEST_F(InputConvTest,restart_load2 )
{
	INPUT.Default();
	std::string input_file = "./support/INPUT";
	INPUT.Read(input_file);
	INPUT.restart_load=true;
	INPUT.dft_functional="b3lyp";
	Input_Conv::Convert();
	EXPECT_EQ(GlobalC::restart.info_load.load_charge,true);
	EXPECT_EQ(GlobalC::restart.info_load.load_H,true);
}

TEST_F(InputConvTest,cell_factor  )
{
	INPUT.Default();
	std::string input_file = "./support/INPUT";
	INPUT.Read(input_file);
	INPUT.calculation="cell-relax";
	INPUT.cell_factor=1.0;
	Input_Conv::Convert();
	EXPECT_EQ(INPUT.cell_factor,2.0);
}

TEST_F(InputConvTest,mixing_beta2 )
{
	INPUT.Default();
	std::string input_file = "./support/INPUT";
	INPUT.Read(input_file);
	INPUT.mixing_beta=2.0;
	std::string output2;
	testing::internal::CaptureStdout();
	Input_Conv::Convert();
	output2 = testing::internal::GetCapturedStdout();
	EXPECT_EQ(GlobalC::CHR_MIX.autoset,false);
}

TEST_F(InputConvTest,mixing_beta )
{
	INPUT.Default();
	std::string input_file = "./support/INPUT";
	INPUT.Read(input_file);
	INPUT.mixing_beta=-10.0;
	Input_Conv::Convert();
	EXPECT_EQ(GlobalC::CHR_MIX.autoset,true);
}

TEST_F(InputConvTest,neighbour  )
{
	INPUT.Default();
	std::string input_file = "./support/INPUT";
	INPUT.Read(input_file);
	INPUT.calculation="nscf";
	INPUT.towannier90=false;
	INPUT.berry_phase=false;
	Input_Conv::Convert();
	EXPECT_EQ(elecstate::ElecStateLCAO::need_psi_grid,false);
}

TEST_F(InputConvTest,neighbour2  )
{
	INPUT.Default();
	std::string input_file = "./support/INPUT";
	INPUT.Read(input_file);
	INPUT.calculation="test_neighbour";
	GlobalV::NPROC =2;
	std::string output2;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(Input_Conv::Convert(), ::testing::ExitedWithCode(0),"");
	output2 = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output2,testing::HasSubstr(" \n !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
		"                         NOTICE                           \n !"
		"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n \n test_neighbour must be done with 1 processor\n"
		" CHECK IN FILE : warning.log\n \n !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
		"                         NOTICE                           \n"
		" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n  |CLASS_NAME---------|NAME---------------|TIME(Sec)-----|CALLS----|AVG------|PER%-------\n -------"
		"---------------------------------------------------------------------------------\n See output information in : \n"
		));
}

TEST_F(InputConvTest, compile)
{
	INPUT.Default();
	std::string input_file = "./support/INPUT";
	INPUT.Read(input_file);
	INPUT.deepks_scf=true;
	std::string output2;
	testing::internal::CaptureStdout();
	EXPECT_EXIT(Input_Conv::Convert(), ::testing::ExitedWithCode(0),"");
	output2 = testing::internal::GetCapturedStdout();
	EXPECT_THAT(output2,testing::HasSubstr(" \n !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
		"                         NOTICE                           \n !"
		"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n \n please compile with DeePKS\n"
		" CHECK IN FILE : warning.log\n \n !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
		"                         NOTICE                           \n"
		" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n  |CLASS_NAME---------|NAME---------------|TIME(Sec)-----|CALLS----|AVG------|PER%-------\n -------"
		"---------------------------------------------------------------------------------\n See output information in : \n"
		));
}

TEST_F(InputConvTest,globalReadinDir )
{
	INPUT.Default();
	std::string input_file = "./support/INPUT";
	INPUT.Read(input_file);
	INPUT.read_file_dir="/root";
	Input_Conv::Convert();
	EXPECT_EQ(GlobalV::global_readin_dir,"/root/");
}

TEST_F(InputConvTest,parse )
{
	INPUT.Default();
	std::string input_file = "./support/INPUT";
	INPUT.Read(input_file);
	ELEC_evolve::td_vext=true;
	Input_Conv::Convert();
	EXPECT_EQ(ELEC_evolve::td_vext_dire_case.size(),0);
}

TEST_F(InputConvTest,parse2 )
{
	INPUT.Default();
	std::string input_file = "./support/INPUT";
	INPUT.Read(input_file);
	GlobalV::ocp=1;
	Input_Conv::Convert();
	EXPECT_EQ(GlobalV::ocp_kb.size(),0);
}

TEST_F(InputConvTest,ParseExpressionDouble)
{
	std::vector<double> vec;
	std::string input = "2*3.5 1.2 0.5";
	Input_Conv::parse_expression(input,vec);
	EXPECT_EQ(vec.size(),4);
	EXPECT_DOUBLE_EQ(vec[0],3.5);
	EXPECT_DOUBLE_EQ(vec[1],3.5);
	EXPECT_DOUBLE_EQ(vec[2],1.2);
	EXPECT_DOUBLE_EQ(vec[3],0.5);
}

#undef private
