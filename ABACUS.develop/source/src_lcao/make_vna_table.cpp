#include "make_vna_table.h"
#include "../src_pw/tools.h"
#include "ylm.h"
#include "lcao_orbitals.h"

double Make_Vna_Table::dr = -1.0;

Make_Vna_Table::Make_Vna_Table()
{
    ntype = 0;
    lmax = 0;
    kmesh = 0;
    Rmax = 0.0;
    dr = 0.0;
    dk = 0.0;

    nlm = 0;
    Rmesh = 0;

    kpoint = new double[1];
    r=new double[1];
    rab=new double[1];
    kab=new double[1];
}

Make_Vna_Table::~Make_Vna_Table()
{
	delete[] kpoint;
    delete[] r;
    delete[] rab;
    delete[] kab;

}

void Make_Vna_Table::allocate (
	const int &ntype_in,
	const int &lmax_in,
	const int &kmesh_in,
	const double &Rmax_in,
	const double &dr_in,
	const double &dk_in)
{
	TITLE("Make_Vna_Table", "allocate");

    this->ntype = ntype_in;// type of elements.
    this->lmax = lmax_in;
    this->kmesh = kmesh_in;
    this->Rmax = Rmax_in;
    this->dr = dr_in;
    this->dk = dk_in;

    assert(ntype > 0);
    assert(lmax >= 0);
    assert(kmesh > 0.0);
    assert(Rmax >= 0.0);
    assert(dr>0.0);
    assert(dk>0.0);

    // calculated from input parameters
    this->nlm = (2*lmax+1) * (2*lmax+1);
    this->Rmesh = static_cast<int>( Rmax/dr ) + 4;
    if(Rmesh%2==0)
    {
        ++Rmesh;
    }

	// rnew the arrays
    delete[] kpoint;
    delete[] r;
    kpoint = new double[kmesh];
    r = new double[Rmesh];

    delete[] rab;
    delete[] kab;
    kab = new double[kmesh];
    rab = new double[Rmesh];

	// set the array values
    for (int ik = 0; ik < kmesh; ik++)
    {
        this->kpoint[ik] = ik * dk_in;
        this->kab[ik] = dk_in;
    }

    for (int ir = 0; ir < Rmesh; ir++)
    {
        this->r[ir] = ir * dr;
        this->rab[ir] = dr;
    }
	return;
}

	
int Make_Vna_Table::get_rmesh(const double &R1, const double &R2)
{
    int rmesh = static_cast<int>((R1+R2)/ Make_Vna_Table::dr) + 5;
    //mohan update 2009-09-08 +1 ==> +5
    //considering interpolation or so on...
    if (rmesh % 2 == 0) rmesh ++;

    if(rmesh <= 0)
    {
        ofs_warning << "\n R1 = " << R1 << " R2 = " << R2;
        ofs_warning << "\n rmesh = " << rmesh;
        WARNING_QUIT("Make_Vna_Table::get_rmesh", "rmesh <= 0");
    }
    return rmesh;
}


void Make_Vna_Table::init_Vna_Tpair(void)
{	
	TITLE("Make_Vna_Table","init_Vna_Tpair");

	this->Vna_nTpairs = this->ntype * this->ntype;
	this->Vna_Tpair.create( this->ntype, this->ntype);
	this->Vna_L2plus1.create( this->ntype, this->ntype);

    int index = 0;
    for (int T1 = 0;  T1 < ntype ; T1++)
    {
        for (int T0 = 0 ; T0 < ntype ; T0++)
        {
             this->Vna_Tpair(T1,T0) = index;
             ++index;

             // the pair < psi | beta >
             // be careful! This is not a symmetry matrix.
             this->Vna_L2plus1(T1,T0) = std::max(ORB.Phi[T1].getLmax(), ORB.Beta[T0].getLmax() )*2+1;

             // there are special situations:
             // for example, two H atom without projector.
             // if we use s orbital,
             // Phi.getLmax = 0,
             // Beta.getLmax < 0,
             // so the value is 1.
             // however, there are no projectors.
             if(Vna_L2plus1(T1,T0) <= 0)
             {
                WARNING_QUIT("Make_Vna_Table::init_paris_nonlocal_type","NL_L2plus1<=0");
             }
        }
    }
	return;
}


void Make_Vna_Table::init_Table_Chi(void)
{
    TITLE("Make_Vna_Table", "init_Table_Chi");
    timer::tick("Make_Vna_Table", "init_Table_Chi",'D');

    // (1) allocate 1st dimension ( overlap, derivative)
    this->Table_VnaR = new double****[2];
    // (2) allocate 2nd dimension ( overlap, derivative)
    this->Table_VnaR[0] = new double*** [this->Vna_nTpairs];
    this->Table_VnaR[1] = new double*** [this->Vna_nTpairs];

    // <1Phi|2Chi>
	for (int T1 = 0;  T1 < ntype ; ++T1) // type 1 is orbital
	{
		const int Lmax1 = ORB.Phi[T1].getLmax();
		for (int T2 = 0 ; T2 < ntype ; ++T2)// type 2 is non-local projector
		{
			const int Lmax2 = ORB.Vna[T2].getLmax();

			// Tpair: type pair.
			const int Tpair=this->Vna_Tpair(T1,T2);
            const int nchi1 = ORB.Phi[T1].getTotal_nchi();

            const int nchi2 = ORB.Vna[T2].getTotal_nchi();
            const int pairs_chi = nchi1 * nchi2;

			if(pairs_chi == 0)continue;

			// init 2nd dimension
			this->Table_VnaR[0][Tpair] = new double** [ pairs_chi ];
			this->Table_VnaR[1][Tpair] = new double** [ pairs_chi ];

			const int T12_2Lplus1 = this->Vna_L2plus1(T1,T2);
			const double Rcut1 = ORB.Phi[T1].getRcut();

			for (int L1 = 0; L1 <= Lmax1; ++L1)
			{
				for (int N1 = 0; N1 < ORB.Phi[T1].getNchi(L1); ++N1)
				{
					for(int L2 = 0; L2 <= Lmax2; ++L2)
					{
						for(int N2 = 0; N2 < ORB.Vna[T2].getNchi(L2); ++N2)
						{
							const double Rcut2 = ORB.Vna[T2].rcut;
							const int Opair = this->Vna_Opair(Tpair,L1,N1,L2,N2);
							assert( Opair < pairs_chi );
							
							// init 3rd dimension
							this->Table_VnaR[0][ Tpair ][ Opair ] = new double *[T12_2Lplus1];
							this->Table_VnaR[1][ Tpair ][ Opair ] = new double *[T12_2Lplus1];
							
							const int rmesh = this->get_rmesh( Rcut1, Rcut2);
							assert( rmesh < this->Rmesh );

							//not all L in T12_2Lplus1 would function
							//const int SL = abs(L1-L2);
							//const int AL = L1+L2;

							for (int L=0; L < T12_2Lplus1 ; L++)
							{
								//Allocation
								this->Table_VnaR[0][Tpair][Opair][L] = new double[rmesh];
								this->Table_VnaR[1][Tpair][Opair][L] = new double[rmesh];

							}	
						}
					}
				}
			}
		}//end T2
	}//end T1	


	return;
}





// phi is the atomic orbital, chi is the projector of Vna
void Make_Vna_Table::init_Vna_Opair(void)
{
    const int lmax = ORB.get_lmax();
    const int nchimax = ORB.get_nchimax();
    assert(lmax+1 > 0);
    assert(nchimax > 0);
    assert(Vna_nTpairs > 0);

    this->Vna_Opair.create( this->Vna_nTpairs, lmax+1, lmax+1, nchimax, nchimax);

    for(int T1=0; T1<ntype; ++T1)
    {
        for(int T2=0; T2<ntype; ++T2)
        {
            const int dim1 = this->Vna_Tpair(T1,T2);
            int index=0;
            for(int L1=0; L1<ORB.Phi[T1].getLmax()+1; ++L1)
            {
                for(int N1=0; N1<ORB.Phi[T1].getNchi(L1); ++N1)
                {
                    for(int L2=0; L2<ORB.Vna[T2].getLmax()+1; ++L2)
                    {
                        for(int N2=0; N2<ORB.Vna[T2].getNchi(L2); ++N2)
                        {
                            this->Vna_Opair(dim1, L1, L2, N1, N2) = index;
                            ++index;
                        }
                    }
                }
            }
        }
    }
    return;
}


void Make_Vna_Table::cal_Vna_PhiChi_R(
	const int &l,
	const Numerical_Orbital_Lm &n1,
	const Neutral_Pot &n2,
	const int &rmesh,
	double* rs,
	double* drs)
{
	timer::tick("Make_Vna_Table","cal_Vna_PhiChi_R");

	/*
	assert(kmesh > 0);
	double *k1_dot_k2 = new double[kmesh];
	for (int ik = 0; ik < kmesh; ik++)
	{
	//	k1_dot_k2[ik] = n1.getPsi_k(ik) * n2.getChi_k(ik);
	}

	double* integrated_func = new double[kmesh];
	double* jl;

	for (int ir = 0; ir < rmesh; ir++)
	{
		ZEROS(integrated_func,kmesh);
		double temp = 0.0;
		jl = this->jlx[l][ir];

		for (int ik = 0; ik < kmesh; ik++)
		{
			integrated_func[ik] = jl[ik] * k1_dot_k2[ik];
		}
		// Call simpson integration
		Mathzone::Simpson_Integral(kmesh,integrated_func,kab,temp);
		rs[ir] = temp * FOUR_PI;

		//drs
		double temp1, temp2;

		if (l > 0)
		{
			jl = this->jlx[l-1][ir];
			for (int ik = 0; ik < kmesh; ik++)
			{
				integrated_func[ik] = jl[ik] * k1_dot_k2[ik] * kpoint[ik];
			}
			Mathzone::Simpson_Integral(kmesh,integrated_func,kab,temp1);
		}

        jl = this->jlx[l+1][ir];

        for (int ik = 0; ik < kmesh; ik++)
		{
            integrated_func[ik] = jl[ik] * k1_dot_k2[ik] * kpoint[ik];
		}

        Mathzone::Simpson_Integral(kmesh,integrated_func,kab,temp2);

        if (l == 0)
        {
            drs[ir] = -FOUR_PI*temp2;
        }
        else
        {
            drs[ir] = FOUR_PI*(temp1*l-(l+1)*temp2)/(2.0*l+1);
        }
    }


    //special case for R=0
    //we store Slm(R) / R**l at the fisrt point, rather than Slm(R)
    if (l > 0)
    {
        ZEROS(integrated_func,kmesh);
        double temp = 0.0;

        for (int ik = 0; ik < kmesh; ik++)
        {
            integrated_func[ik] = k1_dot_k2[ik] * pow (kpoint[ik], l);
        }

        // Call simpson integration
        Mathzone::Simpson_Integral(kmesh,integrated_func,kab,temp);
        rs[0] = FOUR_PI / Mathzone_Add1::dualfac (2*l+1) * temp;
    }

    delete[] integrated_func;
	delete[] k1_dot_k2;
	*/
	timer::tick("Make_Vna_Table","cal_Vna_PhiChi_R");
	return;
}
