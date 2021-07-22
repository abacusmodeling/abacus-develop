#ifdef __MPI
	#include <mpi.h>
#endif
#include <cstring>
#include <iostream>
#include "winput.h"
#include "../src_pw/global.h"


string winput::target;//add 2008-06-04
string winput::wlmr_dir;
double winput::rcut;//a.u.
bool winput::before_iter;//"1" stop before iteration
bool winput::after_iter;
bool winput::begin_stop_flag;
bool winput::end_flag;
string winput::wf_type;
bool winput::build_wf;
int winput::imp_pao;
bool winput::b_out_wf;
bool winput::b_fftwan;//add 2008-07-20
bool winput::b_plot_build;//add 2008-06-04
bool winput::b_plot_atomic;//add 2008-06-04
string winput::trial;//"atomic" or "gauss" 
double winput::bs;//parameters for gauss orbit
double winput::bp;
double winput::px;
double winput::g1;
double winput::g2;
double winput::trunc_ao;
double winput::trunc_wlmr;
double winput::trunc_wan;
double winput::fermi_t;
double winput::clm2_lowest;
int winput::bloch_begin;
int winput::bloch_end;
int winput::sph_proj;//"1" spherical project,"2": first minus atomic orbitals
bool winput::sph_type;//0:Rewrite 1:Skip
bool winput::b_recon;//"1" reconstruction of wannier function
bool winput::b_mix_wf;// add 2008-06-15
double winput::mix_wf;//add 2008-06-13
bool winput::recon_wanq;
bool winput::speed_mode;
bool winput::b_near_atom;
double winput::range0;
double winput::range1;
int winput::L_start;
int winput::L_end;
int winput::atom_start;
int winput::atom_end;
bool winput::plot_wanq;//add 2008-01-26
string winput::plot_option;//(110),[110] etc.
int winput::n_unitcell;//number of unitcell to plot
bool winput::out_all;
bool winput::out_charge;
string winput::charge_type;
bool winput::cal_bands; //for wan  	   wan basis + wan charge
bool winput::cal_bands2;//for semi-wan ;pw basis + wan charge add 2008-4-11
bool winput::cal_dos;
int winput::mesh;
double winput::dr;

bool winput::no_center;
int winput::sum_lm;
bool winput::compare_atomic;
int winput::out_spillage;
string winput::spillage_outdir = "./"; // mohan add 2010-06-07
	
winput::winput()
{
}


winput::~winput()
{
}


void winput::Init(const string &fn)
{
	Default();
	if(test_winput) TITLE("winput","Init");
	//==========================================
	// First readin and check value in root cpu
	// and then bcast value
	//==========================================
	Read(fn);
	Check();

#ifdef __MPI
	Bcast();
#endif
	return;
}

void winput::Read(const string &fn)
{
	TITLE("winput","Read");
	
	if(MY_RANK!=0) return;

	ifstream ifs(fn.c_str(), ios::in);
	if (!ifs)
	{
		//xiaohui move warning 2015-09-30
		//WARNING("winput::Read","Can't find wannier input file.");
		return;
	}
	else
	{
		ofs_running << "Open file : " << fn << endl;
	}

	ifs.clear();
	ifs.seekg(0);

	char word[80], word1[80];
	int ierr = 0;
	ifs.rdstate();

	while (ifs.good())
	{
		ifs >> word;
		ifs.ignore(75, '\n');

		if (std::strcmp(word , "WANNIER_PARAMETERS") == 0){ierr = 1;break;}

		ifs.rdstate();
	}

	if (ierr == 0)
	{
		WARNING("winput::Read","error parameteters title, should be WANNIER_PARAMETERS");
	}

	ifs.rdstate();

	while (ifs.good())
	{
		ifs >> word1;
		INPUT.strtolower(word1, word);
		//parameters for <ESP.wannier> users

		if (strcmp("target",         word) == 0)       {read_value(ifs, target);}
		else if (strcmp("trial",         word) == 0)       {read_value(ifs, trial);}

		else if (strcmp("near_atom",      word) == 0)       {read_value(ifs, b_near_atom);}
		else if (strcmp("range0",		 word) == 0)	{read_value(ifs,range0);}
		else if (strcmp("range1",		 word) == 0)	{read_value(ifs,range1);}

		else if (strcmp("wlmr_dir",      word) == 0)       {read_value(ifs,  wlmr_dir);}
		else if (strcmp("sph_proj",      word) == 0)       {read_value(ifs,  sph_proj);}
		else if (strcmp("sph_type",      word) == 0)       {read_value(ifs,  sph_type);}
		else if (strcmp("recon",		word) == 0)       {read_value(ifs,  b_recon);}
		else if (strcmp("speed_mode",	word) == 0)       {read_value(ifs,  speed_mode);}
		else if (strcmp("recon_wanq",	word) == 0)       {read_value(ifs,  recon_wanq);}
		else if (strcmp("b_mix_wf",		word) == 0)       {read_value(ifs,  b_mix_wf);}
		else if (strcmp("mix_wf",		word) == 0)       {read_value(ifs,  mix_wf);}

		else if (strcmp("wf_type",       word) == 0)       {read_value(ifs,  wf_type);}
		else if (strcmp("build_wf",       word) == 0)       {read_value(ifs,  build_wf);}
		else if (strcmp("imp_pao",       word) == 0)       {read_value(ifs,  imp_pao);}
		else if (strcmp("b_out_wf",         word) == 0)       {read_value(ifs,  b_out_wf);}
		else if (strcmp("b_plot_build",word) == 0)       {read_value(ifs,  b_plot_build);}
		else if (strcmp("b_plot_atomic",word) == 0)       {read_value(ifs,  b_plot_atomic);}

		else if (strcmp("l_start",   	word) == 0)       {read_value(ifs,  L_start);}
		else if (strcmp("l_end", 		word) == 0)       {read_value(ifs,  L_end);}
		else if (strcmp("atom_start",   	word) == 0)   {read_value(ifs,  atom_start);}
		else if (strcmp("atom_end", 		word) == 0)       {read_value(ifs,  atom_end);}

		else if (strcmp("trunc_ao",    word) == 0)       {read_value(ifs,  trunc_ao);}
		else if (strcmp("trunc_wlmr",    word) == 0)       {read_value(ifs,  trunc_wlmr);}
		else if (strcmp("trunc_wan",    word) == 0)       {read_value(ifs,  trunc_wan);}
		else if (strcmp("fermi_t",    	word) == 0)       {read_value(ifs,  fermi_t);}
		else if (strcmp("clm2_lowest",    word) == 0)       {read_value(ifs,  clm2_lowest);}
		else if (strcmp("rcut",          word) == 0)       {read_value(ifs,  rcut);}

		else if (strcmp("plotflag",      word) == 0)       {read_value(ifs,  plot_option);}
		else if (strcmp("unit",          word) == 0)       {read_value(ifs,  n_unitcell);}

		else if (strcmp("b_s",           word) == 0)       {read_value(ifs,  bs);}
		else if (strcmp("b_p",           word) == 0)       {read_value(ifs,  bp);}
		else if (strcmp("px",            word) == 0)       {read_value(ifs,  px);}
		else if (strcmp("g1",            word) == 0)       {read_value(ifs,  g1);}
		else if (strcmp("g2",            word) == 0)       {read_value(ifs,  g2);}

		else if (strcmp("before_iter",	word) == 0)	{read_value(ifs,  before_iter);}
		else if (strcmp("after_iter",      word) == 0) {read_value(ifs,  after_iter);}

		else if (strcmp("begin_stop_flag",	word) == 0)	{read_value(ifs,  begin_stop_flag);}
		else if (strcmp("end_flag",		word) == 0)	{read_value(ifs,  end_flag);}

		else if (strcmp("out_all",		word) == 0)	{read_value(ifs,  out_all);}
		else if (strcmp("out_charge",	word) == 0)	{read_value(ifs,  out_charge);}
		else if (strcmp("compare_atomic",word) == 0)	{read_value(ifs,  compare_atomic);}
		//add 2008-1-26
		else if (strcmp("plot_wanq",	word) == 0)	{read_value(ifs,  plot_wanq);}

		//add 2008-3-10
		else if (strcmp("cal_bands",	word) == 0)	{read_value(ifs,  cal_bands);}
		else if (strcmp("cal_bands2",	word) == 0)	{read_value(ifs,  cal_bands2);}
		else if (strcmp("chgtype",	word) == 0)  {read_value(ifs,  charge_type);}
		else if (strcmp("cal_dos",	word) == 0)	{read_value(ifs,  cal_dos);}

		//add 2008-3-17
		else if (strcmp("bloch_begin",	word) == 0)	{read_value(ifs,  bloch_begin);}
		else if (strcmp("bloch_end",	word) == 0)	{read_value(ifs,  bloch_end);}

		else if (strcmp("mesh",	word) == 0)	{read_value(ifs,  mesh);}
		else if (strcmp("dr",	word) == 0)	{read_value(ifs,  dr);}

		// add 2008-07-20
		else if (strcmp("b_fftwan",	word) == 0)	{read_value(ifs,  b_fftwan);}
		// add 2009-04-19
		else if (strcmp("out_spillage",	word) == 0)	{read_value(ifs,  out_spillage);}
		// add 2010-06-07
		else if (strcmp("spillage_outdir", word) == 0) {read_value(ifs,  spillage_outdir);}
		else
		{
			cout << " The parametr name '" << word << "' is not used." << endl;
			ifs.ignore(150, '\n');
		}

		ifs.rdstate();

		if (ifs.eof() != 0){break;}
		else if (ifs.bad() != 0){cout << "bad parameter. " << endl;break;}
		else if (ifs.fail() != 0){cout << " bad parameter. " << endl;break;}
		else if (ifs.good() == 0){break;}
	}

}

void winput::Default()
{
	if(test_winput) TITLE("winput","Default");
	//========================
	//	part1 : control
	//========================
	target		    = "test";
	rcut		    = 10;
	before_iter	    = false;
	after_iter	    = false;
	begin_stop_flag	= false;
	end_flag		= false;

	//=======================
	//	part2 : Build wf
	//=======================
	wf_type		      = "V";
	build_wf	      = 0;
	imp_pao			  = 0;
	b_out_wf		  = false;
	b_fftwan		  = false;
	b_plot_build	  = false;
	b_plot_atomic	  = false;
		//=========================
		// part2.1 select trial wf
		//=========================
		trial 		= "atomic";//atomic || gauss
		bs			= 2.5; // 
		bp			= 2.0; //  gausss para,eters 
		px			= 2.0; // 
		g1			= 3.0; // D orbital center
		g2 			= 3.0; // D orbital center
		//=======================
		// part2.2 select bands
		//=======================
		bloch_begin = 0;
		bloch_end = 0;


	//========================
	//	part3: spheri & recon
	//========================
	no_center   = false;
	sph_proj	= 0;
	sph_type	= 0;//Rewrite mode
	b_recon		= 0;
	speed_mode	= 1;
	recon_wanq  = 0;
	b_mix_wf	= 0;
	mix_wf		= 0;
		//===============================
		//	part2.1: multi-center shperi
		//===============================
		b_near_atom	= 0;
		range0		= 0;
		range1		= 0;
		//==========================
		//	part2.2: select L,atom
		//==========================
		L_start 	= 0;
		L_end		= 2;
		atom_start	= 0;
		atom_end	= 1;
		//=======================
		//	part2.3: truncation
		//=======================
		trunc_ao	= 6;
		trunc_wlmr	= 14;
		trunc_wan	= 6;
		fermi_t		= 1;
		clm2_lowest = 10e-8;


	//==============================
	//	part3 : Plot wavefunction
	//==============================
	plot_wanq = 0;//add 2008-1-26
	plot_option	= "(110)";
	n_unitcell	= 2;


	//===============================
	// part4 : out_all || out_charge
	//===============================
	out_all		= 0;
	out_charge  = 0;
	compare_atomic = 0;


	//=======================================
	//	part5 : other functions: bands & dos
	//=======================================
	//	add 2008-3-10
	cal_bands = 0;
	cal_bands2= 0;
	charge_type = "planewave";
	cal_dos = 0;

	out_spillage = 0;
	spillage_outdir = "./";

	//add 2008-3-17
	//=========================================
	//	part6 : Uniform mesh ,but not used now
	//=========================================
	mesh = 999;
	dr = 0.01;

	return;
}

void winput::Check(void)
{
	if(test_winput) TITLE("winput","Check");

	if(MY_RANK!=0) return;

	if(CALCULATION=="nscf")
	{
		if(out_charge)
		{
			out_charge = false;
			AUTO_SET("winput::out_charge",out_charge);
		}
	}
	//=====================
	// calculate sum_lm
	//=====================
	for(int il=L_start;il<=L_end;il++)
	{
		sum_lm += il * 2 + 1;
	}

	//xiaohui modify 2013-09-02
	//if(LOCAL_BASIS==1)
	//{
	//	WARNING("winput::Check","Using local basis.");
	//	// turn on
	//	if(!b_recon && !recon_wanq )
	//	{
	//		WARNING("winput::Check","Auto start reconstruction operation.");
	//		b_recon = true;
	//		AUTO_SET("b_recon",b_recon);
	//	}
	//	// turn off
	//	if(before_iter)
	//	{
	//		WARNING("winput::Check","Auto turn down 'before_iter'.");
	//		before_iter = false;
	//		AUTO_SET("before_iter",before_iter);
	//	}
	//	if(after_iter)
	//	{
	//		WARNING("winput::Check","Auto turn down after_iter.");
	//		after_iter = false;
	//		AUTO_SET("after_iter",after_iter);
	//	}
	//	if(build_wf)
	//	{
	//		WARNING("winput::Check","Not available to build wannier functions in local basis");
	//		build_wf = 0;
	//		AUTO_SET("build_wf",build_wf);
	//	}
	//	// some meaning
	//	if(sph_proj>0)
	//	{
	//		if(sph_proj==2)
	//		{
	//			WARNING("winput::Check","Add self site wave functions during reconstruction.");
	//			no_center = true;
	//			AUTO_SET("no_center",no_center);
	//		}
	//		else if(sph_proj==1)
	//		{
	//			WARNING("winput::Check","Not meaning in reconstruction if sph_proj==1");
	//		}
	//		sph_proj=0;
	//		AUTO_SET("sph_proj",sph_proj);
	//	}
	//	if(imp_pao>0)
	//	{
	//		if(imp_pao==2)
	//		{
	//			WARNING("winput::Check","Use improve_pao 2 method, add the center wave functions");
	//		}
	//		if(trunc_wan > 0)
	//		{
	//			WARNING("winput::Check","Use real space truncation. So we must get fft_init started.");
	//			b_fftwan = true;
	//			AUTO_SET("b_fftwan",b_fftwan);
	//		}
	//		clm2_lowest = 0.0;
	//		AUTO_SET("clm2_lowest",clm2_lowest);
	//	}
	//}
	//else if(LOCAL_BASIS==0)
	//{
	//	WARNING("winput::Check","Use plane wave basis.");
	//	// turn off
	//	if(b_recon)
	//	{
	//		WARNING("winput::Check","Auto turn off the reconstruction.");
	//		AUTO_SET("b_recon",b_recon);
	//		b_recon = 0;
	//	}
	//	if(recon_wanq)
	//	{
	//		WARNING("winput::Check","Auto turn off the recon_wanq");
	//		AUTO_SET("recon_wanq",recon_wanq);
	//		recon_wanq = 0;
	//	}
	//	// if turn on 
	//	if(after_iter == true)
	//	{
	//		if(imp_pao == 2)
	//		{
	//			WARNING("winput::Check","Use improve_pao method.");
	//			OUT(ofs_warning,"imp_pao",imp_pao);

	//			WARNING("winput::Check","If use imp_pao>0 ,sph_proj must be 0.");
	//			sph_proj=0;
	//			AUTO_SET("sph_proj",sph_proj);

	//			WARNING("winput::Check","If use imp_pao>0 ,build is no need, can be 0.");
	//			build_wf=0;
	//			AUTO_SET("build_wf",build_wf);

	//			if(trunc_wan > 0 && b_recon)
	//			{
	//				WARNING("winput::Check","Use real space truncation. So we must get fft_init started.");
	//				b_fftwan = true;
	//				AUTO_SET("b_fftwan",b_fftwan);
	//			}
	//			clm2_lowest = 0.0;
	//			AUTO_SET("clm2_lowest",clm2_lowest);
	//		}
	//		else if(sph_proj>0)
	//		{
	//			WARNING("winput::Check","Auto build localized wave functions.");
	//			build_wf = true;
	//			AUTO_SET("build_wf",build_wf);
	//			if(sph_proj == 1) 
	//			{
	//				no_center = false;
	//				AUTO_SET("no_center",no_center);
	//			}
	//			else if(sph_proj == 2) 
	//			{
	//				WARNING("winput::Check","Searching Adjacent without self site.");
	//				no_center = true;
	//				AUTO_SET("no_center",no_center);
	//			}
	//			else
	//			{
	//				WARNING_QUIT("winput::Check","sph_proj must be 0 1 or 2");
	//			}
	//		}

	//		if(build_wf == true)
	//		{
	//			if(bloch_end == 0 || bloch_begin > bloch_end)
	//			{
	//				WARNING_QUIT("winput::Check","Please check your bloch_end");
	//			}
	//			if(bloch_end > NBANDS)
	//			{
	//				WARNING_QUIT("winput::Check","Bloch_end > NBANDS, reset either of them");
	//			}
	//		}
	//	}// end after_iter
	//}

	if(out_charge == true || cal_bands == true || 
			out_all == true || cal_bands2 == true || 
			cal_dos == true)
	{
		end_flag=true;
	}

	if(b_plot_build == true || b_plot_atomic == true)
	{
		if(b_plot_build == true && b_plot_atomic == true)
		{
			WARNING_QUIT("winput::Check()","Plot atomic or plot build wannier functions?");
		}
		else
		{
			b_fftwan = true;
			AUTO_SET("b_fftwan",b_fftwan);
		}
	}

	return;
}

void winput::Print(const string &fn)
{
	if(test_winput) TITLE("winput","Print");

	if(MY_RANK!=0) return;
 
	ofstream ofs(fn.c_str());
	ofs << setiosflags(ios::left);
	ofs << "WANNIER_PARAMETERS" << endl;

	ofs << "#Parameters (General)" << endl;
	OUTP(ofs,"target",target);
	OUTP(ofs,"wlmr_dir",wlmr_dir);
	OUTP(ofs,"rcut",rcut);
	OUTP(ofs,"before_iter",before_iter);
	OUTP(ofs,"after_iter",after_iter);
	OUTP(ofs,"begin_stop_flag",begin_stop_flag);
	OUTP(ofs,"end_flag",end_flag);
	
	ofs << "#Parameters (Build Wannier Functions)" << endl;
	OUTP(ofs,"wf_type",wf_type);
	OUTP(ofs,"build_wf",build_wf);
	OUTP(ofs,"imp_pao",imp_pao);
	OUTP(ofs,"b_out_wf",b_out_wf);
	OUTP(ofs,"b_fftwan",b_fftwan);
	OUTP(ofs,"b_plot_build",b_plot_build);
	OUTP(ofs,"b_plot_atomic",b_plot_atomic);
	
	ofs << "#Parameters (Select trial wave functions)" << endl;
	OUTP(ofs,"trial",trial);
	OUTP(ofs,"bs",bs);
	OUTP(ofs,"bp",bp);
	OUTP(ofs,"px",px);
	OUTP(ofs,"g1",g1);
	OUTP(ofs,"g2",g2);

	ofs << "#Parameters (Select bands)" << endl;
	OUTP(ofs,"bloch_begin",bloch_begin);
	OUTP(ofs,"bloch_end",bloch_end);

	ofs << "#Parameters (Spheri & recon)" << endl;
	OUTP(ofs,"no_center",no_center);	
	OUTP(ofs,"sph_proj",sph_proj);
	OUTP(ofs,"sph_type",sph_type);
	OUTP(ofs,"b_recon",b_recon);
	OUTP(ofs,"speed_mode",speed_mode);
	OUTP(ofs,"recon_wanq",recon_wanq);
	OUTP(ofs,"b_mix_wf",b_mix_wf);
	OUTP(ofs,"mix_wf",mix_wf);

	ofs << "#Parameters (Multi-center spheri)" << endl;
	OUTP(ofs,"b_near_atom",b_near_atom);
	OUTP(ofs,"range0",range0);
	OUTP(ofs,"range1",range1);

	ofs << "#Parameters (Select L, atom)" << endl;
	OUTP(ofs,"L_start",L_start);
	OUTP(ofs,"L_end",L_end);
	OUTP(ofs,"atom_start",atom_start);
	OUTP(ofs,"atom_end",atom_end);

	ofs << "#Parameters (Truncation)" << endl;
	OUTP(ofs,"trunc_ao",trunc_ao);
	OUTP(ofs,"trunc_wlmr",trunc_wlmr);
	OUTP(ofs,"trunc_wan",trunc_wan);
	OUTP(ofs,"fermi_t",fermi_t);
	OUTP(ofs,"clm2_lowest",clm2_lowest);

	ofs << "#Parameters (Plot wave functions)" << endl;
	OUTP(ofs,"plot_wanq",plot_wanq);
	OUTP(ofs,"plot_option",plot_option);
	OUTP(ofs,"n_unitcell",n_unitcell);

	ofs << "#Parameters (out_all || out_charge)" << endl;
	OUTP(ofs,"out_all",out_all);
	OUTP(ofs,"out_charge",out_charge);
	OUTP(ofs,"compare_atomic",compare_atomic);

	ofs << "#Parameters (Other functions: bands & dos)" << endl;
	OUTP(ofs,"cal_bands",cal_bands);
	OUTP(ofs,"cal_bands2",cal_bands2);
	OUTP(ofs,"charge_type",charge_type);
	OUTP(ofs,"cal_dos",cal_dos);
	OUTP(ofs,"out_spillage",out_spillage);
	OUTP(ofs,"spillage_outdir",spillage_outdir);
	
	ofs << "#Parameters (Uniform mesh)" << endl;
	OUTP(ofs,"mesh",mesh);
	OUTP(ofs,"dr",dr);

	OUTP(ofs,"sum_lm",sum_lm);	

	ofs.close();
	return;
}

#ifdef __MPI
void winput::Bcast(void)
{
	if(test_winput) TITLE("winput","Bcast");

	Parallel_Common::bcast_string( target );
	Parallel_Common::bcast_bool( before_iter );
	Parallel_Common::bcast_bool( after_iter );
	Parallel_Common::bcast_bool( begin_stop_flag );
	Parallel_Common::bcast_bool( end_flag );

	Parallel_Common::bcast_double( rcut );
	Parallel_Common::bcast_double( trunc_ao );
	Parallel_Common::bcast_double( trunc_wlmr );
	Parallel_Common::bcast_double( trunc_wan );
	
	Parallel_Common::bcast_string( wlmr_dir );
	Parallel_Common::bcast_string( wf_type );
	Parallel_Common::bcast_bool( build_wf );
	Parallel_Common::bcast_int( imp_pao );
	Parallel_Common::bcast_bool( b_out_wf );
	Parallel_Common::bcast_bool( b_fftwan );
	Parallel_Common::bcast_bool( b_plot_build );
	Parallel_Common::bcast_bool( b_plot_atomic );

	Parallel_Common::bcast_string( trial );
	Parallel_Common::bcast_double( bs );
	Parallel_Common::bcast_double( bp );
	Parallel_Common::bcast_double( px );
	Parallel_Common::bcast_double( g1 );
	Parallel_Common::bcast_double( g2 );

	Parallel_Common::bcast_int( bloch_begin );
	Parallel_Common::bcast_int( bloch_end );
	Parallel_Common::bcast_double( fermi_t );
	Parallel_Common::bcast_double( clm2_lowest );

	Parallel_Common::bcast_int( sph_proj );
	Parallel_Common::bcast_bool( sph_type );

	Parallel_Common::bcast_bool( b_recon );
	Parallel_Common::bcast_bool( b_mix_wf );
	Parallel_Common::bcast_double( mix_wf );
	Parallel_Common::bcast_bool( recon_wanq );

	Parallel_Common::bcast_bool( speed_mode );

	Parallel_Common::bcast_bool( b_near_atom );
	Parallel_Common::bcast_double( range0 );
	Parallel_Common::bcast_double( range1 );

	Parallel_Common::bcast_int( L_start );
	Parallel_Common::bcast_int( L_end );
	Parallel_Common::bcast_int( atom_start );
	Parallel_Common::bcast_int( atom_end );

	Parallel_Common::bcast_bool( plot_wanq );
	Parallel_Common::bcast_string( plot_option );
	Parallel_Common::bcast_int( n_unitcell );
	Parallel_Common::bcast_bool( out_all );
	Parallel_Common::bcast_bool( out_charge );

	Parallel_Common::bcast_string( charge_type );
	Parallel_Common::bcast_bool( cal_bands );
	Parallel_Common::bcast_bool( cal_bands2 );
	Parallel_Common::bcast_bool( cal_dos );
	Parallel_Common::bcast_int( mesh );
	Parallel_Common::bcast_double( dr );

	Parallel_Common::bcast_int( out_spillage );
	Parallel_Common::bcast_string( spillage_outdir );

	Parallel_Common::bcast_bool( no_center );
	Parallel_Common::bcast_int( sum_lm );
	Parallel_Common::bcast_bool( compare_atomic );

	return;
}
#endif
