#include "atom_pseudo.h"

Atom_pseudo::Atom_pseudo()
{
	mbl = new ModuleBase::Vector3<int>[1];
	pseudo_fn = "not_init";
	mass = 0.0;

	for(int is=0;is<4;is++) this->index1_soc[is] = new int[1];
	for(int is=0;is<4;is++) this->index2_soc[is] = new int[1];

}

Atom_pseudo::~Atom_pseudo()
{
	delete[] mbl;

	for(int is=0;is<4;is++) delete[] this->index1_soc[is];
	for(int is=0;is<4;is++) delete[] this->index2_soc[is];
}

// mohan add 2021-05-07
void Atom_pseudo::set_d_so(
	ModuleBase::ComplexMatrix &d_so_in,
	const int &nproj_in,
	const int &nproj_in_so,
	const bool has_so)
{
	if (this->lmax < -1 || this->lmax > 20)
	{
		 ModuleBase::WARNING_QUIT("Numerical_Nonlocal", "bad input of lmax : should be between -1 and 20");
	}

	this->nproj = nproj_in;
	this->nproj_soc=0;

	if(has_so)
	{
		nproj_soc = nproj_in_so;
	}

	assert(nproj <= nproj_in+1); //LiuXh 2016-01-13, 2016-05-16
	assert(nproj >= 0);

//	std::cout << " has_so=" << has_so << std::endl;

	//2016-07-19 begin, LiuXh
	if(!has_so)
	{
		ModuleBase::GlobalFunc::ZEROS(this->non_zero_count_soc, 4);
	}
	else //zhengdy-soc
	{
		this->d_so.create(GlobalV::NSPIN,  nproj_soc+1,  nproj_soc+1);

		// optimize
		for(int is=0;is<4;is++)
		{
			this->non_zero_count_soc[is] = 0;
			delete[] this->index1_soc[is];
			this->index1_soc[is] = new int[nproj_soc * nproj_soc];
			delete[] this->index2_soc[is];
			this->index2_soc[is] = new int[nproj_soc * nproj_soc];
		}

//		std::cout << "lmax=" << lmax << std::endl;

		if(this->lmax > -1)
		{
			if(GlobalV::LSPINORB)
			{
				int is = 0;
				for (int is1 = 0; is1 < 2; is1++)
				{
					for (int is2 = 0; is2 < 2; is2++)
					{
						for (int L1 = 0; L1 < nproj_soc; L1++)
						{
							for (int L2 = 0; L2 < nproj_soc; L2++)
							{
								this->d_so(is, L1, L2) =
									d_so_in(L1 + nproj_soc*is1, L2 + nproj_soc*is2); 

								if(fabs(this->d_so(is, L1, L2).real())>1.0e-8 ||
										fabs(this->d_so(is, L1, L2).imag())>1.0e-8 )
								{
//									std::cout << "tt in atom is=" << is << " L1=" << L1 << " L2=" 
//									<< L2 << " " << d_so(is, L1, L2) << std::endl;

									this->index1_soc[is][non_zero_count_soc[is]] = L1;
									this->index2_soc[is][non_zero_count_soc[is]] = L2;
									this->non_zero_count_soc[is]++;
								}
							}
						}
						is++;
					}
				}
			}
			else
			{
				int is = 0;
				for (int is1 = 0; is1 < 2; is1++)
				{
					for (int is2 = 0; is2 < 2; is2++)
					{
						if(is>=GlobalV::NSPIN) break;
						for (int L1 = 0; L1 < nproj_soc; L1++)
						{
							for (int L2 = 0; L2 < nproj_soc; L2++)
							{
								if(is==1||is==2)
								{
									this->d_so(is, L1, L2) = std::complex<double>(0.0,0.0);
								}
								else
								{
									this->d_so(is, L1, L2)
										= d_so_in(L1 + nproj_soc*is1, L2 + nproj_soc*is2);
								}
								if(abs(this->d_so(is, L1, L2).real())>1.0e-8
										|| abs(this->d_so(is, L1, L2).imag())>1.0e-8)
								{
									this->index1_soc[is][non_zero_count_soc[is]] = L1;
									this->index2_soc[is][non_zero_count_soc[is]] = L2;
									this->non_zero_count_soc[is]++;
								}
							}
						}
						is++;
					}
				}

			}
		}
	}
	//2016-07-19 end, LiuXh
	
	return;
}

void Atom_pseudo::print_atom(std::ofstream &ofs)
{
	if(GlobalV::test_atom) ModuleBase::TITLE("atom_pseudo","print_atom");

	ModuleBase::GlobalFunc::OUT(ofs,"mass",mass);
	ModuleBase::GlobalFunc::OUT(ofs,"pseudo_fn",pseudo_fn);

	return;
}

#include "../src_parallel/parallel_common.h"
#ifdef __MPI
void Atom_pseudo::bcast_atom_pseudo(const int &na)
{
	ModuleBase::TITLE("Atom_pseudo","bcast_atom_pseudo");
	Parallel_Common::bcast_double( mass );
	Parallel_Common::bcast_string( pseudo_fn );

	if(GlobalV::MY_RANK!=0)
	{
		delete[] mbl;
		mbl = new ModuleBase::Vector3<int>[na];
	}

	for(int i=0;i<na;i++)
	{
		Parallel_Common::bcast_int( mbl[i].x );
		Parallel_Common::bcast_int( mbl[i].y );
		Parallel_Common::bcast_int( mbl[i].z );
	}
	return;
}

void Atom_pseudo::bcast_atom_pseudo2(void)
{
	ModuleBase::TITLE("Atom_pseudo","bcast_atom_pseudo2");
// == pseudo_h ==
//int
	Parallel_Common::bcast_int( lmax );
	Parallel_Common::bcast_int( mesh );
	Parallel_Common::bcast_int( nchi );
	Parallel_Common::bcast_int( nbeta );
	Parallel_Common::bcast_int( nv );
	Parallel_Common::bcast_int( zv );

//double
	Parallel_Common::bcast_double( etotps );
	Parallel_Common::bcast_double( ecutwfc );
	Parallel_Common::bcast_double( ecutrho );

// bool
	Parallel_Common::bcast_bool( tvanp );
	Parallel_Common::bcast_bool( nlcc );
	Parallel_Common::bcast_bool( has_so );

//std::string
	Parallel_Common::bcast_string( psd );
	Parallel_Common::bcast_string( pp_type );
	Parallel_Common::bcast_string( xc_func );

	if(GlobalV::MY_RANK!=0)
	{
		delete[] jjj;
		delete[] els;
		delete[] lchi;
		delete[] oc;
		delete[] jchi;
		delete[] nn;
		jjj = new double [nbeta];
		els = new std::string[nchi];
		lchi = new int [nchi];
		oc = new double[nchi];
		jchi = new double[nchi];
		nn = new int[nchi];
	}
	
	Parallel_Common::bcast_double( jjj, nbeta );
	Parallel_Common::bcast_string( els, nchi);
	Parallel_Common::bcast_int( lchi, nchi);
	Parallel_Common::bcast_double( oc, nchi);
	Parallel_Common::bcast_double( jchi, nchi);
	Parallel_Common::bcast_int( nn, nchi);
// == end of pseudo_h

// == pseudo_atom ==
	Parallel_Common::bcast_int( msh );
	Parallel_Common::bcast_double( rcut );
	if(GlobalV::MY_RANK != 0)
	{
		assert(mesh!=0);
		delete[] r;
		delete[] rab;
		delete[] rho_atc;
		delete[] rho_at;
		r = new double[mesh];
		rab = new double[mesh];
		rho_atc = new double[mesh];
		rho_at = new double[mesh];
		chi.create( nchi,mesh );
	}
	
	Parallel_Common::bcast_double( r, mesh );
	Parallel_Common::bcast_double( rab, mesh );
	Parallel_Common::bcast_double( rho_atc, mesh );
	Parallel_Common::bcast_double( rho_at, mesh );
	Parallel_Common::bcast_double( chi.c, nchi * mesh );
// == end of pseudo_atom ==
	
// == pseudo_vl ==
	if(GlobalV::MY_RANK != 0)
	{
		delete[] vloc_at;
		vloc_at = new double[mesh];
	}
	Parallel_Common::bcast_double( vloc_at, mesh);
// == end of pseudo_vl ==

// == pseudo_nc ==
	if(GlobalV::MY_RANK != 0)
	{
		delete[] lll;
		lll = new int[nbeta];
	}
	Parallel_Common::bcast_int( lll, nbeta );
	Parallel_Common::bcast_int( kkbeta );
	Parallel_Common::bcast_int( nh );

	int nr,nc;
	if(GlobalV::MY_RANK == 0)
	{
		nr = betar.nr;
		nc = betar.nc;
	}
	Parallel_Common::bcast_int( nr );
	Parallel_Common::bcast_int( nc );

	if(GlobalV::MY_RANK != 0)
	{
		betar.create(nr,nc);
		dion.create(nbeta, nbeta);
	}

	// below two 'bcast_double' lines of codes seem to have bugs,
	// on some computers, the code will stuck here for ever.
	// mohan note 2021-04-28
	Parallel_Common::bcast_double( dion.c , nbeta * nbeta);
	Parallel_Common::bcast_double( betar.c, nr * nc );
// == end of psesudo_nc ==

	return;
}

#endif
