//==========================================================
// AUTHOR : Lixin He, mohan
// DATE : 2008-11-08
//==========================================================
#include "global.h"			// only out

#include "pseudopot_cell_vnl.h"
#include "tools.h"
#include "wavefunc.h"

pseudopot_cell_vnl::pseudopot_cell_vnl()
{
}

pseudopot_cell_vnl::~pseudopot_cell_vnl()
{
}

//==========================================================
// see allocate_nlpot.f90
//==========================================================
void pseudopot_cell_vnl::init(const int ntype, const bool allocate_vkb)
{
	TITLE("pseudopot_cell_vnl", "init");
	timer::tick("ppcell_vnl","init",'C');

	ofs_running << "\n SETUP NONLOCAL PSEUDOPOTENTIALS IN PLANE WAVE BASIS" << endl;

	int it = 0;
//----------------------------------------------------------
// MEMBER VARIABLE :
// NAME : lmaxkb(max angular momentum,(see pseudo_h))
//----------------------------------------------------------
	this->lmaxkb = - 1;
	for (it = 0;it < ntype; it++)
	{
		ofs_running << " " << ucell.atoms[it].label << " non-local projectors:" << endl;
		for (int ibeta = 0; ibeta < ucell.atoms[it].nbeta; ibeta++) 
		{
			ofs_running << " projector " << ibeta+1 << " L=" << ucell.atoms[it].lll[ibeta] <<  endl;
			this->lmaxkb = max( this->lmaxkb, ucell.atoms[it].lll[ibeta]);
		}
	}

//----------------------------------------------------------
// MEMBER VARIABLE :
// NAME : nhm(max number of different beta functions per atom)
//----------------------------------------------------------
	this->nhm = 0;
	for (it=0;it<ntype;it++)
	{	
		this->nhm = std::max(nhm, ucell.atoms[it].nh);
	}

//----------------------------------------------------------
// MEMBER VARIABLE :
// NAME : nkb(total number of beta functions, with struct.fact.)
//----------------------------------------------------------
	this->nkb = 0;
	for (it=0; it<ntype; it++)
	{
		this->nkb += ucell.atoms[it].nh * ucell.atoms[it].na;
	}

	OUT(ofs_running,"TOTAL NUMBER OF NONLOCAL PROJECTORS",nkb);

	if( this->nhm > 0 )
	{
		indv.create(ntype, this->nhm);
		nhtol.create(ntype, this->nhm);
		nhtolm.create(ntype, this->nhm);
	//	nhtoj.create(ntype, this->nhm);
		deeq.create(NSPIN, ucell.nat, this->nhm, this->nhm);
	//	qq.create(ntype, this->nhm, this->nhm);
		dvan.create(ntype, this->nhm, this->nhm);
		becsum.create(NSPIN, ucell.nat, this->nhm * (this->nhm + 1) / 2);
	}
	else
	{
		ofs_running << "\n nhm = 0, not allocate some matrix.";
	}

	// nqxq = ((sqrt(gcutm)+sqrt(xqq[1]*xqq[1]+xqq[2]*xqq[2]+xqq[3]*xqq[3])/
	// dq+4)*cell_factor;
	this->lmaxq = 2 * this->lmaxkb + 1;

	if (nkb > 0 && allocate_vkb )
	{
		vkb.create(nkb, wf.npwx);
	}

	//this->nqx = 10000;		// calculted in allocate_nlpot.f90
	NQX = this->calculate_nqx(pw.ecutwfc,DQ);
	// nqx = (sqrt(ecutwfc)/dq+4)*cell_factor;
	//
	//	nbrx is defined in constant.h
	//  max number of beta functions
	this->tab.create(ntype, nbrx, NQX);

	// nchix is defined in constant.h
	// nchix : max number of atomic wavefunctions per atom
	this->tab_at.create(ntype, nchix, NQX);

	if(test_pp > 1)
	{
		cout 
		<< "\n     ntype   = " << ntype
		<< "\n     lmaxkb  = " << this->lmaxkb
		<< "\n     nhm     = " << this->nhm  
		<< "\n     nkb     = " << this->nkb
		<< "\n     nbrx    = " << nbrx
		<< "\n     nchix   = " << nchix 
		<< endl;
	}

	timer::tick("ppcell_vnl","init",'C');
	return;
}

//----------------------------------------------------------
// MEMBER FUNCTION NAME : 
// NAME : pseudopot_cell_vnl
// Calculates beta functions (Kleinman-Bylander projectors), with
// structure factor, for all atoms, in reciprocal space
// from init_us_2.f90
//----------------------------------------------------------
void pseudopot_cell_vnl::getvnl(const int &ik)
{
	if(test_pp) TITLE("pseudopot_cell_vnl","getvnl");
	timer::tick("pp_cell_vnl","getvnl");

	if(lmaxkb < 0) 
	{
		return;
	}

	const int npw = kv.ngk[ik];
	int ig, ia, nb, ih;
	matrix vkb1(nhm, npw);
	double *vq = new double[npw];
	const int x1= (lmaxkb + 1)*(lmaxkb + 1);

	matrix ylm(x1, npw);
	Vector3<double> *gk = new Vector3<double>[npw];
	for (ig = 0;ig < npw;ig++) 
	{
		gk[ig] = wf.get_1qvec_cartesian(ik, ig);
	}


	Mathzone::Ylm_Real(x1, npw, gk, ylm);

	int jkb = 0;
	for(int it = 0;it < ucell.ntype;it++)
	{
		if(test_pp>1) OUT("it",it);
		// calculate beta in G-space using an interpolation table
		const int nbeta = ucell.atoms[it].nbeta;
		const int nh = ucell.atoms[it].nh;

		if(test_pp>1) OUT("nbeta",nbeta);

		for (nb = 0;nb < nbeta;nb++)
		{
			if(test_pp>1) OUT("ib",nb);
			for (ig = 0;ig < npw;ig++)
			{
				const double gnorm = gk[ig].norm() * ucell.tpiba;

//				cout << "\n gk[ig] = " << gk[ig].x << " " << gk[ig].y << " " << gk[ig].z;
//				cout << "\n gk.norm = " << gnorm;
				
				vq [ig] = Mathzone::Polynomial_Interpolation(
						this->tab, it, nb, NQX, DQ, gnorm );
			} // enddo

			// add spherical harmonic part
			for (ih = 0;ih < nh;ih++)
			{
				if (nb == this->indv(it, ih))
				{
					const int lm = static_cast<int>( nhtolm(it, ih) );
					for (ig = 0;ig < npw;ig++)
					{
						vkb1(ih, ig) = ylm(lm, ig) * vq [ig];
					}
				}
			} // end ih
		} // end nbeta

		// vkb1 contains all betas including angular part for type nt
		// now add the structure factor and factor (-i)^l
		for (ia=0; ia<ucell.atoms[it].na; ia++) 
		{
			complex<double> *sk = wf.get_sk(ik, it, ia);
			for (ih = 0;ih < nh;ih++)
			{
				complex<double> pref = pow( NEG_IMAG_UNIT, nhtol(it, ih));	//?
				for (ig = 0;ig < npw;ig++)
				{
					this->vkb(jkb, ig) = vkb1(ih, ig) * sk [ig] * pref;
				}
				++jkb;
			} // end ih
			delete [] sk;
		} // end ia
	} // enddo

	delete [] gk;
	delete [] vq;
	timer::tick("pp_cell_vnl","getvnl");
	return;
} // end subroutine getvnl

void pseudopot_cell_vnl::init_vnl(void)
{
	if(test_pp) TITLE("pseudopot_cell_vnl","init_vnl");
	timer::tick("ppcell_vnl","init_vnl");

	//from init_us_1
	//   a) For each non vanderbilt pseudopotential it computes the D and
	//      the betar in the same form of the Vanderbilt pseudopotential.
	//   b) It computes the indices indv which establish the correspondence
	//      nh <-> beta in the atom
	//   c) It computes the indices nhtol which establish the correspondence
	//      nh <-> angular momentum of the beta function
	//   d) It computes the indices nhtolm which establish the correspondence
	//      nh <-> combined (l,m) index for the beta function.

	//   For each pseudopotential we initialize the indices nhtol, nhtolm,
	//   nhtoj, indv, and if the pseudopotential is of KB type we initialize
	// the atomic D terms
	this->dvan.zero_out();

	for(int it=0;it<ucell.ntype;it++)
	{
		int BetaIndex=0;
		const int Nprojectors = ucell.atoms[it].nh;
		for (int ib=0; ib<ucell.atoms[it].nbeta; ib++)
		{
			const int l = ucell.atoms[it].lll [ib];
//			const int j = static_cast<int>(ucell.atoms[it].jjj [nb]);
			for(int m=0; m<2*l+1; m++)
			{
				this->nhtol(it,BetaIndex) = l;
				this->nhtolm(it,BetaIndex) = l*l + m;
//				nhtoj(it,BetaIndex) = j;
				this->indv(it,BetaIndex) = ib;
				++BetaIndex;
			}
		}

		//    From now on the only difference between KB and US pseudopotentials
		//    is in the presence of the q and Q functions.
		//    Here we initialize the D of the solid
		for (int ip=0; ip<Nprojectors; ip++)
		{
			for (int ip2=0; ip2<Nprojectors; ip2++)
			{
				if ( this->nhtol (it, ip) == nhtol (it, ip2) &&
				     this->nhtolm(it, ip) == nhtolm(it, ip2) )
				{
					const int ir = static_cast<int>( indv(it, ip ) );
					const int is = static_cast<int>( indv(it, ip2) );
					this->dvan(it, ip, ip2) = ucell.atoms[it].dion(ir, is);
				}
			} 
		} 
	} 

	//   h) It fills the interpolation table for the beta functions
	/**********************************************************
	// He Lixin: this block is used for non-local potential
	// fill the interpolation table tab
	************************************************************/

	const double pref = FOUR_PI / sqrt(ucell.omega);
	this->tab.zero_out();
	ofs_running<<"\n Init Non-Local PseudoPotential table : ";
	for (int it = 0;it < ucell.ntype;it++)  
	{
		const int nbeta = ucell.atoms[it].nbeta;
		int kkbeta = ucell.atoms[it].kkbeta;

		//mohan modify 2008-3-31
		//mohan add kkbeta>0 2009-2-27
		if ( (kkbeta%2 == 0) && kkbeta>0 )
		{
			kkbeta--;
		}

		double *jl = new double[kkbeta];
		double *aux  = new double[kkbeta];

		for (int ib = 0;ib < nbeta;ib++)
		{
			const int l = ucell.atoms[it].lll[ib];
			for (int iq=0; iq<NQX; iq++)  
			{
				const double q = iq * DQ;
				Mathzone::Spherical_Bessel(kkbeta, ucell.atoms[it].r, q, l, jl);

				for (int ir = 0;ir < kkbeta;ir++)
				{
					aux[ir] = ucell.atoms[it].betar(ib, ir) *
					          jl[ir] * ucell.atoms[it].r[ir];
				} 
				double vqint;
				Mathzone::Simpson_Integral(kkbeta, aux, ucell.atoms[it].rab, vqint);
				this->tab(it, ib, iq) = vqint * pref;
			} 
		} 
		delete[] aux;
		delete[] jl;
	}
	timer::tick("ppcell_vnl","init_vnl");
	ofs_running << "\n Init Non-Local-Pseudopotential done." << endl;
	return;
}


void pseudopot_cell_vnl::print_vnl(ofstream &ofs)
{
	out.printr3_d(ofs, " tab : ", tab);
}

int pseudopot_cell_vnl::calculate_nqx(const double &ecutwfc,const double &dq)
{
	int points_of_table = static_cast<int>( sqrt(ecutwfc)/dq + 4 ) ;//* cell_factor;
//	cout<<"\n nqx = "<< points_of_table;
//----------------------------------------------------------
// EXPLAIN : Plus 1 because the formula is transfered from 
// fortran code.
//----------------------------------------------------------
	return points_of_table + 1; 
}
