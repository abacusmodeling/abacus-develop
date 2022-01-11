#include "LCAO_deepks_test.h"

namespace Test_Deepks
{
	Grid_Driver GridD(GlobalV::test_deconstructor,GlobalV::test_grid_driver,GlobalV::test_grid);
}

test_deepks::test_deepks()
{}

test_deepks::~test_deepks()
{}


void test_deepks::check_dstable(void)
{
	OGT.talpha.print_Table_DSR(ORB);
}

void test_deepks::check_psialpha(void)
{
	std::vector<int> na;
	na.resize(ucell.ntype);
	for(int it=0;it<ucell.ntype;it++)
	{
		na[it] = ucell.atoms[it].na;
	}
	ld.init(ORB,
		ucell.nat,
		ucell.ntype,
		na);

	ld.build_psialpha(GlobalV::FORCE,
		ucell,
		ORB,
		Test_Deepks::GridD,
		ParaO.trace_loc_row,
		ParaO.trace_loc_col,
		OGT);

	ld.check_psialpha(GlobalV::FORCE,
		ucell,
		ORB,
		Test_Deepks::GridD,
		ParaO.trace_loc_row,
		ParaO.trace_loc_col,
		OGT);
}

void test_deepks::read_dm(void)
{
	ifstream ifs("dm");
	dm.create(GlobalV::NLOCAL,GlobalV::NLOCAL);

	for (int mu=0;mu<GlobalV::NLOCAL;mu++)
	{
		for (int nu=0;nu<GlobalV::NLOCAL;nu++)
		{
			double c;
			ifs >> c;
			dm(mu,nu)=c;
		}
	}
}

void test_deepks::read_dm_k(const int nks)
{
	dm_k.resize(nks);
	stringstream ss;
	for(int ik=0;ik<nks;ik++)
	{
        ss.str("");
        ss<<"dm_"<<ik;
        ifstream ifs(ss.str().c_str());
		dm_k[ik].create(GlobalV::NLOCAL,GlobalV::NLOCAL);

		for (int mu=0;mu<GlobalV::NLOCAL;mu++)
		{
			for (int nu=0;nu<GlobalV::NLOCAL;nu++)
			{
				std::complex<double> c;
				ifs >> c;
				dm_k[ik](mu,nu)=c;
			}
		}
	}
}

void test_deepks::check_pdm(void)
{
	if(GlobalV::GAMMA_ONLY_LOCAL)
	{
		this->read_dm();
		this->ld.cal_projected_DM(dm,
			ucell,
			ORB,
			Test_Deepks::GridD,
			ParaO.trace_loc_row,
			ParaO.trace_loc_col);
	}
	else
	{
		this->read_dm_k(kv.nkstot);
		this->ld.cal_projected_DM_k(dm_k,
			ucell,
			ORB,
			Test_Deepks::GridD,
			ParaO.trace_loc_row,
			ParaO.trace_loc_col,
			kv.nkstot,
			kv.kvec_d);		
	}
	this->ld.check_projected_dm();
}

void test_deepks::check_gdmx(void)
{
	this->ld.init_gdmx(ucell.nat);
	if(GlobalV::GAMMA_ONLY_LOCAL)
	{
		this->ld.cal_gdmx(dm,
			ucell,
			ORB,
			Test_Deepks::GridD,
			ParaO.trace_loc_row,
			ParaO.trace_loc_col);
	}
	else
	{
		this->ld.cal_gdmx_k(dm_k,
			ucell,
			ORB,
			Test_Deepks::GridD,
			ParaO.trace_loc_row,
			ParaO.trace_loc_col,
			kv.nkstot,
			kv.kvec_d);			
	}
	this->ld.check_gdmx(ucell.nat);
}

void test_deepks::check_descriptor(void)
{
	this->ld.cal_descriptor();
	this->ld.check_descriptor(ucell.nat);
}

void test_deepks::check_gvx(void)
{
	this->ld.cal_gvx(ucell.nat);
	this->ld.check_gvx(ucell.nat);
}

void test_deepks::check_edelta(void)
{
	this->ld.load_model("model.ptg");
	if(GlobalV::GAMMA_ONLY_LOCAL)
	{
		this->ld.allocate_V_delta(ucell.nat,ParaO.nloc);
	}
	else
	{
		this->ld.allocate_V_delta(ucell.nat,ParaO.nloc,kv.nkstot);
	}
	this->ld.cal_gedm(ucell.nat);

	ofstream ofs("E_delta");
	ofs << std::setprecision(10) << this->ld.E_delta;
	
	this->ld.check_gedm();
}

void test_deepks::check_e_deltabands(void)
{
	if(GlobalV::GAMMA_ONLY_LOCAL)
	{
		this->ld.cal_e_delta_band(dm,
			ParaO.trace_loc_row,
			ParaO.trace_loc_col,
			ParaO.nrow);
	}
	else
	{
		this->folding_nnr(kv);
		this->ld.cal_e_delta_band_k(dm_k,
			ParaO.trace_loc_row,
			ParaO.trace_loc_col,
			kv.nkstot,
			ParaO.nrow,
			ParaO.ncol);
	}

	ofstream ofs("E_delta_bands.dat");
	ofs << std::setprecision(10) << this->ld.e_delta_band;
}

void test_deepks::check_v_delta()
{
	if(GlobalV::GAMMA_ONLY_LOCAL)
	{
		this->ld.add_v_delta(ucell,
            ORB,
            Test_Deepks::GridD,
            ParaO.trace_loc_row,
            ParaO.trace_loc_col,
            ParaO.nrow,
            ParaO.ncol);
		this->ld.check_v_delta(ParaO.nrow,ParaO.ncol);
	}
	else
	{
		this->cal_nnr();
		this->ld.allocate_V_deltaR(nnr);
		this->ld.add_v_delta_k(ucell,
        	ORB,
            Test_Deepks::GridD,
            ParaO.trace_loc_row,
			ParaO.trace_loc_col,
			nnr);
		this->ld.check_v_delta_k(nnr);
	}
}

void test_deepks::check_f_delta()
{
	ModuleBase::matrix svnl_dalpha;
	svnl_dalpha.create(3,3);
	if(GlobalV::GAMMA_ONLY_LOCAL)
	{
		ld.cal_f_delta_gamma(dm,
            ucell,
            ORB,
            Test_Deepks::GridD,
            ParaO.trace_loc_row,
			ParaO.trace_loc_col,
            1, svnl_dalpha);
	}
	else
	{
		ld.cal_f_delta_k(dm_k,
			ucell,
            ORB,
            Test_Deepks::GridD,
            ParaO.trace_loc_row,
			ParaO.trace_loc_col,
			kv.nkstot,
			kv.kvec_d,
			1,svnl_dalpha);
	}
	ld.check_f_delta(ucell.nat, svnl_dalpha);
}