#include "wavefunc.h"
#include "global.h"
#include "algorithms.h"
#include "../src_pw/wavefunc_in_pw.h"
//xiaohui add 2013 -08-01
#include "winput.h"
#ifdef __FP
//#include "../src_develop/src_wannier/manipulation.h"
#endif

wavefunc::wavefunc()
{
	allocate_ekb = false;
}

wavefunc::~wavefunc()
{ 
	if(test_deconstructor)
	{
		cout << " ~wavefunc()" << endl;
	}
	if(allocate_ekb)
	{
		// bug still remains, hard to find!
		// it might be somewhere out there,
		// may be in cdiaghg.
		// I don't know why.......
		// mohan 2010-08-08
		//for(int ik=0; ik<kv.nks-1; ik++) delete[] ekb[ik];
		//delete[] ekb;
	}
}

#ifdef __FP
void wavefunc::init_local(void)
{
    TITLE("wavefunc","init_local");
    const int nks = kv.nks;
    this->npwx = this->setupIndGk(pw, nks);

	// mohan add 2009-12-24
	this->ekb = new double*[nks];
	for(int ik=0; ik<nks; ik++)
	{
		ekb[ik] = new double[NBANDS];
		ZEROS(ekb[ik],NBANDS);
	}
	this->allocate_ekb = true;

    this->wg.create(nks, NBANDS);   // the weight of each k point and band
    Memory::record("wavefunc","ekb",nks*NBANDS,"double");
    Memory::record("wavefunc","wg",nks*NBANDS,"double");

    if(test_wf)ofs_running << " Allocate : EigenValue; Weight;" << endl;
    return;
}
#endif

void wavefunc::init(const int nks)
{	
	static bool done_once = false;
	if(done_once)
	{
		return;
	}
	else
	{
		done_once = true;
	}

    TITLE("wavefunc","init");
    this->npwx = this->setupIndGk(pw, nks);
	OUT(ofs_running,"npwx",npwx);

    assert(npwx > 0);
    assert(nks > 0);
    assert(NBANDS > 0);

    delete[] g2kin;
    this->g2kin = new double[npwx];   // [npw],kinetic energy
    ZEROS(g2kin, npwx);
    Memory::record("wavefunc","g2kin",npwx,"double");
	if(test_wf)OUT(ofs_running,"g2kin allocation","Done");

    // if use spin orbital, double nks to allocate evc and wanf2.
    int prefactor = 1;
	
	this->ekb = new double*[nks];
	for(int ik=0; ik<nks; ik++)
	{
		ekb[ik] = new double[NBANDS];
		ZEROS(ekb[ik], NBANDS);
	}
	this->allocate_ekb = true;
	
    this->wg.create(nks, NBANDS);	// the weight of each k point and band
    Memory::record("wavefunc","et",nks*NBANDS,"double");
    Memory::record("wavefunc","wg",nks*NBANDS,"double");
	if(test_wf)OUT(ofs_running, "et allocation","Done"); 
	if(test_wf)OUT(ofs_running, "wg allocation","Done");

    delete[] evc;
    delete[] wanf2;

    const int nks2 = nks * prefactor;
    
	if(CALCULATION=="nscf" && wf.mem_saver==1)
	{
		// mohan add 2010-09-07
		// only available for first principle now.
		this->evc = new ComplexMatrix[1];
		this->wanf2 = new ComplexMatrix[1];
		
		evc[0].create(NBANDS*prefactor, npwx);	
#ifdef __FP
		//if(LOCAL_BASIS && LINEAR_SCALING==0) xiaohui modify 2013-09-02
		if(BASIS_TYPE=="lcao_in_pw") //xiaohui add 2013-09-02
		{
			wanf2[0].create(NLOCAL*prefactor, npwx);
    		cout << " Memory for wanf2 (MB): " << 
			Memory::record("wavefunc","wanf2",(NLOCAL*prefactor)*npwx,"complexmatrix") << endl;
		}
#endif	
    	cout << " MEMORY FOR PSI (MB)  : " << 
		Memory::record("wavefunc","evc",(NBANDS*prefactor)*npwx,"complexmatrix") << endl;
	}
	else
	{
		this->evc = new ComplexMatrix [nks2];
    	this->wanf2 = new ComplexMatrix [nks2];

    	for (int ik = 0; ik < nks2; ik++)
    	{
        	this->evc[ik].create(NBANDS*prefactor, npwx);

			//Mohan add 2010-1-10
#ifdef __FP
	    //xiaohui add 2013 -08-01 
    	    //if (LOCAL_BASIS || winput::out_spillage==2) xiaohui modify 2013-09-02
	    if((BASIS_TYPE=="lcao" || BASIS_TYPE=="lcao_in_pw") || winput::out_spillage==2) //xiaohui add 2013-09-02
	    {
		this->wanf2[ik].create(NLOCAL, npwx);
	    }
#endif	
    	};

    	cout << " MEMORY FOR PSI (MB)  : " << 
		Memory::record("wavefunc","evc",(nks*prefactor)*(NBANDS*prefactor)*npwx,"complexmatrix") << endl;
	}

	if(test_wf)
	{	
		OUT(ofs_running,"evc allocation","Done");
		//if(LOCAL_BASIS) xiaohui modify 2013-09-02
		if(BASIS_TYPE=="lcao" || BASIS_TYPE=="lcao_in_pw") //xiaohui add 2013-09-02
		{
			OUT(ofs_running,"wanf2 allocation","Done");
		}
	}
    //showMemStats();
    return;
}

//===================================================================
// This routine computes an estimate of the start_ wavefunctions
// from superposition of atomic wavefunctions or random wave functions.
//===================================================================
void wavefunc::wfcinit(void)
{
    TITLE("wavefunc","wfcinit");
    timer::tick("wavefunc","wfcinit",'C');


	this->init_at_1();	
    this->wfcinit_k();

    en.demet = 0.0;

    //================================
    // Occupations are computed here
    //================================
	//if(LOCAL_BASIS == 0 && LINEAR_SCALING == 0) xiaohui modify 2013-09-02
	if(BASIS_TYPE=="pw") //xiaohui add 2013-09-02
	{
		// mohan fix bug 2011-02-25,
		// in nscf, occupations is not needed,
		if(CALCULATION=="scf" || CALCULATION=="md" || CALCULATION=="relax") //pengfei 2014-10-13
		{
    		Occupy::calculate_weights();
		}
	}
    if (test_wf>2)
    {
        out.printrm(ofs_running, " wg  ",  wg);
        this->check_psi(evc);
    }

    timer::tick("wavefunc","wfcinit",'C');
    return;
}

int wavefunc::get_starting_nw(void)const
{
    if (start_wfc == "file")
    {
        WARNING_QUIT("wfcinit_k","\n start_ wfc from file: not implemented yet!");
        //**********************************************************************
        // ... read the wavefunction into memory (if it is not done in c_bands)
        //**********************************************************************
    }
    else if (start_wfc == "atomic")
    {
        if (ucell.natomwfc >= NBANDS)
        {
            if(test_wf)ofs_running << " Start wave functions are all pseudo atomic wave functions." << endl;
        }
        else
        {
            if(test_wf)ofs_running << " Start wave functions are atomic + "
            << NBANDS - ucell.natomwfc
            << " random wave functions." << endl;
        }
        return max(ucell.natomwfc,  NBANDS);
    }
    else if (start_wfc == "random")
    {
        if(test_wf)ofs_running << " Start wave functions are all random." << endl;
        return NBANDS;
    }
    else
    {
        WARNING_QUIT("get_starting_nw","Don't know what to do! Please Check source code!");
    }
}



// not full test.
void wavefunc::PAO_in_pw_k(const int &ik, ComplexMatrix &wvf)
{
	TITLE("wavefunc","PAO_in_pw_k");
	//assert(LINEAR_SCALING==0); xiaohui modify 2013-09-02
	//assert(LOCAL_BASIS==2); xiaohui modify 2013-09-02. Attention! Here why assert LOCAL_BASIS=2 ???

	this->atomic_wfc(ik, this->npw, ucell.lmax_ppwf, wvf, ppcell.tab_at, NQX, DQ);
}

void wavefunc::LCAO_in_pw_k(const int &ik, ComplexMatrix &wvf)
{
	TITLE("wavefunc","LCAO_in_pw_k");
	timer::tick("wavefunc","LCAO_in_pw_k",'G');

#ifdef __FP
	//assert(LOCAL_BASIS==4); xiaohui modify 2013-09-01
	assert(BASIS_TYPE=="lcao_in_pw"); //xiaohui add 2013-09-01. Attention! How about "BASIS_TYPE=="lcao""???
	
	static bool ltable = false;
	if(!ltable)
	{
		this->table_local.create(ucell.ntype, ucell.nmax_total, NQX);
		Wavefunc_in_pw::make_table_q(ORB.orbital_file, this->table_local);
		ltable = true;
	}
	
	Wavefunc_in_pw::produce_local_basis_in_pw(ik, wvf, this->table_local);

	//-------------------------------------------------------------
	// (2) diago to get wf.ekb, then the weights can be calculated.
	//-------------------------------------------------------------
	hm.init();
	hm.init_k(ik);
	
	//hm.cinitcgg(ik ,NLOCAL, NBANDS, wvf, wvf, ekb[ik]);
//	for(int ib=0; ib<NBANDS; ib++)
//	{
//		cout << " ib=" << ib << " e=" << ekb[ik][ib] << endl;
//	}

//	DONE(ofs_running,"CONSTRUCT_LOCAL_BASIS_IN_PW");
#endif
	timer::tick("wavefunc","LCAO_in_pw_k",'G');
	return;
}

void wavefunc::diago_PAO_in_pw_k(const int &ik, ComplexMatrix &wvf)
{
	TITLE("wavefunc","diago_PAO_in_pw_k");

	hm.init_k(ik);
    this->diago_PAO_in_pw_k2(ik, wvf);

	return;
}

void wavefunc::diago_PAO_in_pw_k2(const int &ik, ComplexMatrix &wvf)
{
	TITLE("wavefunc","diago_PAO_in_pw_k2");
	// (6) Prepare for atmoic orbitals or random orbitals
	const int starting_nw = this->get_starting_nw();
    assert(starting_nw > 0);

	ComplexMatrix wfcatom(starting_nw, npwx);
	if(test_wf)OUT(ofs_running, "starting_nw", starting_nw);
	if(start_wfc=="atomic")
	{
		this->atomic_wfc(ik, this->npw, ucell.lmax_ppwf, wfcatom, ppcell.tab_at, NQX, DQ);
		//====================================================
		// If not enough atomic wfc are available, complete
		// with random wfcs
		//====================================================
		this->random(wfcatom, ucell.natomwfc, NBANDS, ik);
	}
	else if(start_wfc=="random")
	{
		this->random(wfcatom,0,NBANDS,ik);
	}

	// (7) Diago with cg method.
	double *etatom  = new double[starting_nw];
	ZEROS(etatom, starting_nw);
	//if(DIAGO_TYPE == "cg") xiaohui modify 2013-09-02
	if(KS_SOLVER=="cg") //xiaohui add 2013-09-02
	{
		hm.cinitcgg(ik ,starting_nw, NBANDS, wfcatom, wfcatom, etatom);
	}

	/*
	ofs_running << " " << "ik = " << ik << " Bands(eV)" << endl;
	for (int ib=0;ib<starting_nw;ib++)
	{
		ofs_running << " " << setw(15) << etatom[ib]*Ry_to_eV << endl;
	}
	*/

	assert(wvf.nr <= wfcatom.nr);
	for (int ib=0; ib<NBANDS; ib++)
	{
		for (int ig=0; ig<this->npwx; ig++)
		{
			wvf(ib, ig) = wfcatom(ib, ig);
		}
	}

	delete[] etatom;
}

void wavefunc::wfcinit_k(void)
{
	TITLE("wavefunc","wfcinit_k");

	if(mem_saver) return;

	for(int ik=0; ik<kv.nks; ik++)
	{
		//if(LOCAL_BASIS==4) xiaohui modify 2013-09-02
		if(BASIS_TYPE=="lcao_in_pw") //xiaohui add 2013-09-02
		{
			// just get the numerical local basis wave functions
			// in plane wave basis
			this->LCAO_in_pw_k(ik, wf.wanf2[ik]);
		}
		//else if(LOCAL_BASIS==3) xiaohui modify 2013-09-02.
		//{ 
		//	// prepare for Jlq, gauss, slater-type
		//	// which is anylatic basis.
		//	WARNING_QUIT("wavefunc::wfcinit_k","not ready for local_basis=3");
		//} xiaohui modify 2013-09-02. Attention...

		//else if(LOCAL_BASIS==2) xiaohui modify 2013-09-02
		//{
		//	// prepare for using PAO directly.
		//	this->PAO_in_pw_k(ik, wf.wanf2[ik]);
		//	WARNING_QUIT("wavefunc::wfcinit_k","not ready for local_basis=2");
		//} xiaohui modify 2013-09-02. Attention...

		//else if(LOCAL_BASIS==1) xiaohui modify 2013-09-02
		//{
		//	// prepare for wannier functions.
		//	WARNING_QUIT("wavefunc::wfcinit_k","not ready for local_basis=1");
		//} xiaohui modify 2013-09-02. Attention...

		//else if(LOCAL_BASIS==0) xiaohui modify 2013-09-02
		else if (BASIS_TYPE=="pw") //xiaohui add 2013-09-02
		{
			// get the wave functions 
			// by first diagolize PAO
			// wave functions.
			this->diago_PAO_in_pw_k(ik, wf.evc[ik]);
		}
	}

    return;
}
