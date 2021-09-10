#include "FORCE_gamma.h"
#include "../src_pw/global.h"
#include "dftu.h"  //Quxin add for DFT+U on 20201029

Force_LCAO_gamma::Force_LCAO_gamma ()
{}

Force_LCAO_gamma::~Force_LCAO_gamma ()
{}

// be called in force_lo.cpp
void Force_LCAO_gamma::ftable_gamma (
	const bool isforce,
	const bool isstress,
	ModuleBase::matrix& foverlap,
	ModuleBase::matrix& ftvnl_dphi,
	ModuleBase::matrix& fvnl_dbeta,	
	ModuleBase::matrix& fvl_dphi,
	ModuleBase::matrix& soverlap,
	ModuleBase::matrix& stvnl_dphi,
	ModuleBase::matrix& svnl_dbeta,
	ModuleBase::matrix& svl_dphi)
{
    ModuleBase::TITLE("Force_LCAO_gamma", "ftable");
    ModuleBase::timer::tick("Force_LCAO_gamma","ftable_gamma");
    
    // allocate DSloc_x, DSloc_y, DSloc_z
    // allocate DHloc_fixed_x, DHloc_fixed_y, DHloc_fixed_z
    this->allocate_gamma();

    // calculate the 'energy density matrix' here.
    this->cal_foverlap(isforce, isstress, foverlap, soverlap);

    if(INPUT.new_dm>0)
    {
        this->cal_ftvnl_dphi(GlobalC::LOC.wfc_dm_2d.dm_gamma, isforce, isstress, ftvnl_dphi, stvnl_dphi);
        this->cal_fvnl_dbeta(GlobalC::LOC.wfc_dm_2d.dm_gamma, isforce, isstress, fvnl_dbeta, svnl_dbeta);
        this->cal_fvl_dphi(GlobalC::LOC.wfc_dm_2d.dm_gamma, isforce, isstress, fvl_dphi, svl_dphi);

        //quxin added for DFT+U
        if(INPUT.dft_plus_u) GlobalC::dftu.force_stress();
    }
    else
    {
        ModuleBase::timer::tick("Force_LCAO_gamma","cal_dm_grid");
        // calculate the 'density matrix' here.
        ModuleBase::matrix dm2d;
		dm2d.create(GlobalV::NSPIN, GlobalC::ParaO.nloc);
        ModuleBase::Memory::record ("Force_LCAO_gamma", "dm2d", GlobalC::ParaO.nloc*GlobalV::NSPIN, "double");    

        bool with_energy = false;

		this->set_EDM_gamma(dm2d, with_energy);

        ModuleBase::timer::tick("Force_LCAO_gamma","cal_dm_grid");

        this->cal_ftvnl_dphi(dm2d, isforce, isstress, ftvnl_dphi, stvnl_dphi);
        this->cal_fvnl_dbeta(dm2d, isforce, isstress, fvnl_dbeta, svnl_dbeta);

        //quxin added for DFT+U
        if(INPUT.dft_plus_u) GlobalC::dftu.force_stress();

        // calculate < dphi | V | phi > on real space grid.
        this->cal_fvl_dphi(dm2d, isforce, isstress, fvl_dphi, svl_dphi);

    }
	if(isforce)
	{
        Parallel_Reduce::reduce_double_pool( foverlap.c, foverlap.nr * foverlap.nc);
        Parallel_Reduce::reduce_double_pool( ftvnl_dphi.c, ftvnl_dphi.nr * ftvnl_dphi.nc);
        Parallel_Reduce::reduce_double_pool( fvnl_dbeta.c, fvnl_dbeta.nr * fvnl_dbeta.nc);
        Parallel_Reduce::reduce_double_pool( fvl_dphi.c, fvl_dphi.nr * fvl_dphi.nc);
	}
    if(isstress)
    {
        Parallel_Reduce::reduce_double_pool( soverlap.c, soverlap.nr * soverlap.nc);
        Parallel_Reduce::reduce_double_pool( stvnl_dphi.c, stvnl_dphi.nr * stvnl_dphi.nc);
        Parallel_Reduce::reduce_double_pool( svnl_dbeta.c, svnl_dbeta.nr * svnl_dbeta.nc);
        Parallel_Reduce::reduce_double_pool( svl_dphi.c, svl_dphi.nr * svl_dphi.nc);
    }

    // delete DSloc_x, DSloc_y, DSloc_z
    // delete DHloc_fixed_x, DHloc_fixed_y, DHloc_fixed_z
    this->finish_ftable_gamma();

    ModuleBase::timer::tick("Force_LCAO_gamma","ftable_gamma");
    return;
}

void Force_LCAO_gamma::allocate_gamma(void)
{
    ModuleBase::TITLE("Force_LCAO_gamma","allocate_gamma");
    ModuleBase::timer::tick("Force_LCAO_gamma","allocate_gamma");

    // need to calculate the derivative in build_ST_new
    bool cal_deri = true;

    //calculate dS in LCAO
    //liaochen add on 2010/7/12
    //save the results in dense matrix by now.
    //GlobalC::ParaO.nloc: number of H elements in this proc.
    GlobalC::LM.DSloc_x = new double [GlobalC::ParaO.nloc];
    GlobalC::LM.DSloc_y = new double [GlobalC::ParaO.nloc];
    GlobalC::LM.DSloc_z = new double [GlobalC::ParaO.nloc];
    ModuleBase::GlobalFunc::ZEROS(GlobalC::LM.DSloc_x, GlobalC::ParaO.nloc);
    ModuleBase::GlobalFunc::ZEROS(GlobalC::LM.DSloc_y, GlobalC::ParaO.nloc);
    ModuleBase::GlobalFunc::ZEROS(GlobalC::LM.DSloc_z, GlobalC::ParaO.nloc);
    //allocate stress part in gamma_only-line, added by zhengdy-stress
    if(GlobalV::STRESS)
    {
        GlobalC::LM.DSloc_11 = new double [GlobalC::ParaO.nloc];
        GlobalC::LM.DSloc_12 = new double [GlobalC::ParaO.nloc];
        GlobalC::LM.DSloc_13 = new double [GlobalC::ParaO.nloc];
        GlobalC::LM.DSloc_22 = new double [GlobalC::ParaO.nloc];
        GlobalC::LM.DSloc_23 = new double [GlobalC::ParaO.nloc];
        GlobalC::LM.DSloc_33 = new double [GlobalC::ParaO.nloc];
        ModuleBase::GlobalFunc::ZEROS(GlobalC::LM.DSloc_11, GlobalC::ParaO.nloc);
        ModuleBase::GlobalFunc::ZEROS(GlobalC::LM.DSloc_12, GlobalC::ParaO.nloc);
        ModuleBase::GlobalFunc::ZEROS(GlobalC::LM.DSloc_13, GlobalC::ParaO.nloc);
        ModuleBase::GlobalFunc::ZEROS(GlobalC::LM.DSloc_22, GlobalC::ParaO.nloc);
        ModuleBase::GlobalFunc::ZEROS(GlobalC::LM.DSloc_23, GlobalC::ParaO.nloc);
        ModuleBase::GlobalFunc::ZEROS(GlobalC::LM.DSloc_33, GlobalC::ParaO.nloc);
        GlobalC::LM.DHloc_fixed_11 = new double [GlobalC::ParaO.nloc];
        GlobalC::LM.DHloc_fixed_12 = new double [GlobalC::ParaO.nloc];
        GlobalC::LM.DHloc_fixed_13 = new double [GlobalC::ParaO.nloc];
        GlobalC::LM.DHloc_fixed_22 = new double [GlobalC::ParaO.nloc];
        GlobalC::LM.DHloc_fixed_23 = new double [GlobalC::ParaO.nloc];
        GlobalC::LM.DHloc_fixed_33 = new double [GlobalC::ParaO.nloc];
        ModuleBase::GlobalFunc::ZEROS (GlobalC::LM.DHloc_fixed_11, GlobalC::ParaO.nloc);
        ModuleBase::GlobalFunc::ZEROS (GlobalC::LM.DHloc_fixed_12, GlobalC::ParaO.nloc);
        ModuleBase::GlobalFunc::ZEROS (GlobalC::LM.DHloc_fixed_13, GlobalC::ParaO.nloc);
        ModuleBase::GlobalFunc::ZEROS (GlobalC::LM.DHloc_fixed_22, GlobalC::ParaO.nloc);
        ModuleBase::GlobalFunc::ZEROS (GlobalC::LM.DHloc_fixed_23, GlobalC::ParaO.nloc);
        ModuleBase::GlobalFunc::ZEROS (GlobalC::LM.DHloc_fixed_33, GlobalC::ParaO.nloc);
    }
    //calculate dS in LCAO basis
    // tips: build_ST_new --> GlobalC::ParaO.set_force 
    //ModuleBase::timer::tick("Force_LCAO_gamma","build_S_new");
    GlobalC::UHM.genH.build_ST_new ('S', cal_deri, GlobalC::ucell);
    //ModuleBase::timer::tick("Force_LCAO_gamma","build_S_new");

    ModuleBase::Memory::record("force_lo", "dS", GlobalC::ParaO.nloc*3, "double");

    //calculate dT in LCAP
    //allocation dt
    //liaochen add on 2010/7/12
    GlobalC::LM.DHloc_fixed_x = new double [GlobalC::ParaO.nloc];
    GlobalC::LM.DHloc_fixed_y = new double [GlobalC::ParaO.nloc];
    GlobalC::LM.DHloc_fixed_z = new double [GlobalC::ParaO.nloc];
    ModuleBase::GlobalFunc::ZEROS (GlobalC::LM.DHloc_fixed_x, GlobalC::ParaO.nloc);
    ModuleBase::GlobalFunc::ZEROS (GlobalC::LM.DHloc_fixed_y, GlobalC::ParaO.nloc);
    ModuleBase::GlobalFunc::ZEROS (GlobalC::LM.DHloc_fixed_z, GlobalC::ParaO.nloc);
    
    //calculate dT
    //calculate T + VNL(P1) in LCAO basis
    //ModuleBase::timer::tick("Force_LCAO_gamma","build_T_new");
    GlobalC::UHM.genH.build_ST_new ('T', cal_deri, GlobalC::ucell);
    //ModuleBase::timer::tick("Force_LCAO_gamma","build_T_new");
    //test_gamma(GlobalC::LM.DHloc_fixed_x, "dHloc_fixed_x T part");
    
    //GlobalC::UHM.genH.build_Nonlocal_beta (cal_deri);
    //ModuleBase::timer::tick("Force_LCAO_gamma","build_Nonlocal_mu");
    GlobalC::UHM.genH.build_Nonlocal_mu (cal_deri);
    //ModuleBase::timer::tick("Force_LCAO_gamma","build_Nonlocal_mu");
    //test_gamma(GlobalC::LM.DHloc_fixed_x, "dHloc_fixed_x Vnl part");

    ModuleBase::Memory::record("force_lo", "dTVNL", GlobalC::ParaO.nloc*3, "double");

    ModuleBase::timer::tick("Force_LCAO_gamma","allocate_gamma");
    return;
}

void Force_LCAO_gamma::finish_ftable_gamma(void)
{
    delete [] GlobalC::LM.DSloc_x;
    delete [] GlobalC::LM.DSloc_y;
    delete [] GlobalC::LM.DSloc_z;
    delete [] GlobalC::LM.DHloc_fixed_x;
    delete [] GlobalC::LM.DHloc_fixed_y;
    delete [] GlobalC::LM.DHloc_fixed_z;
    if(GlobalV::STRESS)//added by zhengdy-stress
    {
        delete [] GlobalC::LM.DSloc_11;
        delete [] GlobalC::LM.DSloc_12;
        delete [] GlobalC::LM.DSloc_13;
        delete [] GlobalC::LM.DHloc_fixed_11;
        delete [] GlobalC::LM.DHloc_fixed_12;
        delete [] GlobalC::LM.DHloc_fixed_13;
        delete [] GlobalC::LM.DSloc_22;
        delete [] GlobalC::LM.DSloc_23;
        delete [] GlobalC::LM.DSloc_33;
        delete [] GlobalC::LM.DHloc_fixed_22;
        delete [] GlobalC::LM.DHloc_fixed_23;
        delete [] GlobalC::LM.DHloc_fixed_33;
    }
    return;
}


void Force_LCAO_gamma::test_gamma(double* mm, const std::string &name)
{
    std::cout << "\n PRINT " << name << std::endl;
    std::cout << std::setprecision(6) << std::endl;
    for(int i=0; i<GlobalV::NLOCAL; i++)
    {
        for(int j=0; j<GlobalV::NLOCAL; j++)
        {
            if( abs(mm[i*GlobalV::NLOCAL+j])>1.0e-5)
            {
                std::cout << std::setw(12) << mm[i*GlobalV::NLOCAL+j];
            }
            else
            {
                std::cout << std::setw(12) << "0";
            }
        }
        std::cout << std::endl;
    }
    return;
}
