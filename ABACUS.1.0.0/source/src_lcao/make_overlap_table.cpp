#include "make_overlap_table.h"
#include "lcao_orbitals.h"
#include "../src_global/sph_bessel.h"

#include <stdexcept>
double Make_Overlap_Table::dr = -1.0;

Make_Overlap_Table::Make_Overlap_Table()
{
//	cout << " \n Make_Overlap_Table::Make_Overlap_Table" << endl;
	destroy_sr = false;
	destroy_tr = false;
	destroy_nr = false;
	destroy_jlx = false;

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

Make_Overlap_Table::~Make_Overlap_Table()
{
	delete[] kpoint;
	delete[] r;
	delete[] rab;
	delete[] kab;
}

void Make_Overlap_Table::allocate
(
 	const int &ntype_in,
    const int &lmax_in,
    const int &kmesh_in,
	const double &Rmax_in,
    const double &dr_in,
    const double &dk_in
)
{
	TITLE("Make_Overlap_Table", "allocate");

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

//	OUT(ofs_running,"lmax",lmax);
//	OUT(ofs_running,"Rmax (Bohr)",Rmax);
//	OUT(ofs_running,"dr (Bohr)",dr);
//	OUT(ofs_running,"dk",dk);
//	OUT(ofs_running,"nlm",nlm);
//	OUT(ofs_running,"kmesh",kmesh);
	
	delete[] kpoint;
	delete[] r;
	kpoint = new double[kmesh];
	r = new double[Rmesh];

	delete[] rab;
	delete[] kab;
	kab = new double[kmesh];
	rab = new double[Rmesh];

	for (int ik = 0; ik < kmesh; ik++)
	{
		kpoint[ik] = ik * dk_in;
		kab[ik] = dk_in;
	}

	for (int ir = 0; ir < Rmesh; ir++)
	{
		r[ir] = ir * dr;
		rab[ir] = dr;
	}

//	OUT(ofs_running,"allocate kpoint, r, rab, kab","Done");
	return;
}

int Make_Overlap_Table::get_rmesh(const double &R1, const double &R2)
{
	int rmesh = static_cast<int>((R1+R2)/ Make_Overlap_Table::dr) + 5;
	//mohan update 2009-09-08 +1 ==> +5
	//considering interpolation or so on...
	if (rmesh % 2 == 0) rmesh ++;
	
	if(rmesh <= 0)
	{
		ofs_warning << "\n R1 = " << R1 << " R2 = " << R2;
		ofs_warning << "\n rmesh = " << rmesh;
		WARNING_QUIT("Make_Overlap_Table::get_rmesh", "rmesh <= 0");
	}
	return rmesh;
}

void Make_Overlap_Table::cal_ST_Phi12_R
(
 	const int &job,
    const int &l,
    const Numerical_Orbital_Lm &n1,
    const Numerical_Orbital_Lm &n2,
    const int &rmesh,
    double* rs,
	double* drs
) const
{
//	TITLE("Make_Overlap_Table","cal_ST_Phi12_R");
	timer::tick("Make_Overlap_Table", "cal_ST_Phi12_R");

	double* k1_dot_k2 = new double[kmesh];
	switch(job)
	{
		case 1: // calculate overlap
		for (int ik = 0; ik < kmesh; ik++)
		{
			k1_dot_k2[ik] = n1.getPsi_k(ik) * n2.getPsi_k(ik);
		}
		break;

		case 2: // calculate kinetic energy
		for (int ik = 0; ik < kmesh; ik++)
		{
			k1_dot_k2[ik] = n1.getPsi_k(ik) * n2.getPsi_k(ik) * this->kpoint[ik] * this->kpoint[ik];
		}
		break;
	}

//	Mathzone_Add1::Sbt_new (3, l, r, dr, rmesh, kpoint, dk, kmesh, k1_dot_k2, 2, rs);
//	for (int ir = 0; ir < rmesh; ir++) rs[ir] *= FOUR_PI;
	
	//test
//	double tmporg = Mathzone_Add1::uni_simpson (k1_dot_k2, kmesh, dk);
//	tmporg *= FOUR_PI;

	//liaochen modify on 2010/4/21
	//for origin 

	//Drs
	//djl = (l*j(l-1) - (l+1)j(l+1))/(2l+1)
	/*
	if (l == 0) 
	{
		double* tmp1 = new double[rmesh];
		ZEROS (tmp1, rmesh);
		Mathzone_Add1::Sbt_new (3, 1, r, dr, rmesh, kpoint, dk, kmesh, k1_dot_k2, 1, tmp1);
		for (int ir = 0; ir < rmesh; ir++) drs[ir] = -FOUR_PI*tmp1[ir];
		delete[] tmp1;
	}
	else
	{
		double* tmp1 = new double[rmesh];
		double* tmp2 = new double[rmesh];
		ZEROS (tmp1, rmesh);
		ZEROS (tmp2, rmesh);
		Mathzone_Add1::Sbt_new (3, l-1, r, dr, rmesh, kpoint, dk, kmesh, k1_dot_k2, 1, tmp1);
		Mathzone_Add1::Sbt_new (3, l+1, r, dr, rmesh, kpoint, dk, kmesh, k1_dot_k2, 1, tmp2);
		for (int ir = 0; ir < rmesh; ir++) drs[ir] = FOUR_PI*(l*tmp1[ir]-(l+1)*tmp2[ir])/(2.0*l+1);
		delete[] tmp1;
		delete[] tmp2;
	}
	*/
	
	
	//previous version
	
	
	double* integrated_func = new double[kmesh];
	
	double* jl;
//	double* jl = new double[kmesh];
	
	for (int ir = 0; ir < rmesh; ir++)
	{
		ZEROS(integrated_func,kmesh);
//		ZEROS(jl,kmesh);
		double temp = 0.0;
		// Generate Spherical Bessel Function
//		Mathzone::Spherical_Bessel(this->kmesh,this->kpoint,this->r[ir], l, jl);
		jl = this->jlx[l][ir];
		
		for (int ik = 0; ik < kmesh; ik++)
		{
			integrated_func[ik] = jl[ik] * k1_dot_k2[ik];
		}
		// Call simpson integration
		Mathzone::Simpson_Integral(kmesh,integrated_func,kab,temp);
		rs[ir] = temp * FOUR_PI ;

		//drs
		double temp1, temp2;
		
		if (l > 0)
		{
	//		ZEROS(jl,kmesh);
	//		Mathzone::Spherical_Bessel(this->kmesh,this->kpoint,this->r[ir], l-1, jl);
			jl = this->jlx[l-1][ir];
		
			for (int ik = 0; ik < kmesh; ik++)
			{
				integrated_func[ik] = jl[ik] * k1_dot_k2[ik] * kpoint[ik];
			}

			Mathzone::Simpson_Integral(kmesh,integrated_func,kab,temp1);
		}
		
//		ZEROS(jl,kmesh);
//		Mathzone::Spherical_Bessel(this->kmesh,this->kpoint,this->r[ir], l+1, jl);
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

	//liaochen modify on 2010/4/22
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
//	delete [] jl;
	delete [] integrated_func;
	

	/*
	double* integrated_func1 = new double[kmesh];
	double* integrated_func2 = new double[kmesh];
	
	for (int ir = 0; ir < rmesh; ir++)
	{
		ZEROS(integrated_func1, kmesh);
		ZEROS(integrated_func2, kmesh);
//		ZEROS(jl, kmesh);
//		ZEROS(djldr, kmesh);
		
		// Generate Spherical Bessel Function
//		Mathzone::Spherical_Bessel(this->kmesh,this->kpoint,this->r[ir], l, jl, djldr);
		for (int ik = 0; ik < kmesh; ik++)
		{
			integrated_func1[ik] = this->jlx[ik][ir][l] * k1_dot_k2[ik];
			integrated_func2[ik] = this->djlx[ik][ir][l] * k1_dot_k2[ik] * this->kpoint[ik];
		}
		// Call simpson integration
		Mathzone::Simpson_Integral(kmesh,integrated_func1,kab,rs[ir]);
		Mathzone::Simpson_Integral(kmesh,integrated_func2,kab,drs[ir]);
	
		rs[ir] *= FOUR_PI;
		drs[ir] *= FOUR_PI;
	}
	

	delete[] integrated_func1;
	delete[] integrated_func2;
//	delete[] jl;
	*/

	delete[] k1_dot_k2;
	timer::tick("Make_Overlap_Table", "cal_ST_Phi12_R");
	
	return;
}


void Make_Overlap_Table::cal_VNL_PhiBeta_R(
        const int &l,
        const Numerical_Orbital_Lm &n1,
        const Numerical_Nonlocal_Lm &n2,
        const int &rmesh,
        double *rs,
		double *drs)
{
	timer::tick ("Make_Overlap_Table", "VNL_PhiBeta_R");

	assert(kmesh > 0);

	//start calc	
    double *k1_dot_k2 = new double[kmesh];

	for (int ik = 0; ik < kmesh; ik++)
	{
		k1_dot_k2[ik] = n1.getPsi_k(ik) * n2.getBeta_k(ik);
	}

//	Mathzone_Add1::Sbt_new (3, l, r, dr, rmesh, kpoint, dk, kmesh, k1_dot_k2, 2, rs);
//	for (int ir = 0; ir < rmesh; ir++) rs[ir] *= FOUR_PI;
	
	
	//Drs
	//djl = (l*j(l-1) - (l+1)j(l+1))/(2l+1)
	/*
	if (l == 0) 
	{
		double* tmp1 = new double[rmesh];
		ZEROS (tmp1, rmesh);
		Mathzone_Add1::Sbt_new (3, 1, r, dr, rmesh, kpoint, dk, kmesh, k1_dot_k2, 1, tmp1);
		for (int ir = 0; ir < rmesh; ir++) drs[ir] = -FOUR_PI*tmp1[ir];
		delete[] tmp1;
	}
	else
	{
		double* tmp1 = new double[rmesh];
		double* tmp2 = new double[rmesh];
		ZEROS (tmp1, rmesh);
		ZEROS (tmp2, rmesh);
		Mathzone_Add1::Sbt_new (3, l-1, r, dr, rmesh, kpoint, dk, kmesh, k1_dot_k2, 1, tmp1);
		Mathzone_Add1::Sbt_new (3, l+1, r, dr, rmesh, kpoint, dk, kmesh, k1_dot_k2, 1, tmp2);
		for (int ir = 0; ir < rmesh; ir++) drs[ir] = FOUR_PI*(l*tmp1[ir]-(l+1)*tmp2[ir])/(2.0*l+1);
		delete[] tmp1;
		delete[] tmp2;
	}
	*/
	
	//previous version
	double* integrated_func = new double[kmesh];
	double* jl;
//	double* jl = new double[kmesh];
	
	for (int ir = 0; ir < rmesh; ir++)
	{
		ZEROS(integrated_func,kmesh);
//		ZEROS(jl,kmesh);
		double temp = 0.0;
		// Generate Spherical Bessel Function
//		Mathzone::Spherical_Bessel(this->kmesh,this->kpoint,this->r[ir], l, jl);
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
//			ZEROS(jl,kmesh);
//			Mathzone::Spherical_Bessel(this->kmesh,this->kpoint,this->r[ir], l-1, jl);
			jl = this->jlx[l-1][ir];
					
			for (int ik = 0; ik < kmesh; ik++)
				integrated_func[ik] = jl[ik] * k1_dot_k2[ik] * kpoint[ik];

			Mathzone::Simpson_Integral(kmesh,integrated_func,kab,temp1);
		}
		
//		ZEROS(jl,kmesh);
//		Mathzone::Spherical_Bessel(this->kmesh,this->kpoint,this->r[ir], l+1, jl);
		jl = this->jlx[l+1][ir];
				
		for (int ik = 0; ik < kmesh; ik++)
			integrated_func[ik] = jl[ik] * k1_dot_k2[ik] * kpoint[ik];
		
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
	
	//liaochen modify on 2010/4/22
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
	
	delete [] integrated_func;
//	delete [] jl;
	
	/*
	double* integrated_func1 = new double[kmesh];
	double* integrated_func2 = new double[kmesh];
	
	for (int ir = 0; ir < rmesh; ir++)
	{
		ZEROS(integrated_func1, kmesh);
		ZEROS(integrated_func2, kmesh);
//		ZEROS(jl, kmesh);
//		ZEROS(djldr, kmesh);
		
		// Generate Spherical Bessel Function
//		Mathzone::Spherical_Bessel(this->kmesh,this->kpoint,this->r[ir], l, jl, djldr);
		for (int ik = 0; ik < kmesh; ik++)
		{
			integrated_func1[ik] = this->jlx[ik][ir][l] * k1_dot_k2[ik];
			integrated_func2[ik] = this->djlx[ik][ir][l] * k1_dot_k2[ik] * this->kpoint[ik];
		}
		// Call simpson integration
		Mathzone::Simpson_Integral(kmesh,integrated_func1,kab,rs[ir]);
		Mathzone::Simpson_Integral(kmesh,integrated_func2,kab,drs[ir]);
	
		rs[ir] *= FOUR_PI;
		drs[ir] *= FOUR_PI;
	}

	delete [] integrated_func1;
	delete [] integrated_func2;
	*/

	delete[] k1_dot_k2;
	timer::tick ("Make_Overlap_Table", "VNL_PhiBeta_R");
	return;
}

void Make_Overlap_Table::init_Table( const int &job0 )
{
	TITLE("Make_Overlap_Table", "init_Table");
	timer::tick("Make_Overlap_Table", "init_Table",'D');
	const int ntype = ORB.get_ntype();
	assert( Make_Overlap_Table::dr > 0.0);
	assert( OV_nTpairs>0);

	// init 1st dimension
	switch( job0 )
	{
		case 1:
		// the second dimension stands for S(R) and dS(R)/dR
		this->Table_SR = new double****[2];
		for(int ir = 0; ir < 2; ir++)
		{
			this->Table_SR[ir] = new double***[ this->OV_nTpairs ];
		}
		break;

		case 2:
		this->Table_TR = new double****[2];
		for(int ir = 0; ir < 2; ir++)
		{
			this->Table_TR[ir] = new double***[ this->OV_nTpairs ];
		}
		break;

		case 3:
		this->Table_SR = new double****[2];
		this->Table_TR = new double****[2];
		for(int ir = 0; ir < 2; ir++)
		{
			this->Table_SR[ir] = new double***[ this->OV_nTpairs ];
			this->Table_TR[ir] = new double***[ this->OV_nTpairs ];
		}
		break;
	}

	for (int T1 = 0;  T1 < ntype ; T1++)
	{
		// Notice !! T2 start from T1
		// means that T2 >= T1
		for (int T2 = T1 ; T2 < ntype ; T2++)
		{
			// get the bigger lmax between two types
			const int Tpair=this->OV_Tpair(T1,T2);
			const int Lmax1 = ORB.Phi[T1].getLmax();
			const int Lmax2 = ORB.Phi[T2].getLmax();

			//L2plus1 could be reduced by considering Gaunt Coefficient
			//remain to be modified
			//??????
			const int lmax_now = std::max( Lmax1, Lmax2 );
			

			//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			// mohan add 2011-03-07
			// I think the lmax_now should be judged from two 
			// orbitals, not atom type!!!!!!!!!!!!!
			// there are space that can imporve the efficiency.
			//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			
			const int L2plus1 =  2*lmax_now + 1;

			const int nchi1 = ORB.Phi[T1].getTotal_nchi();
			const int nchi2 = ORB.Phi[T2].getTotal_nchi();
			const int pairs_chi = nchi1 * nchi2;

			// init 2nd dimension
			switch( job0 )
			{
				case 1:
				this->Table_SR[0][ Tpair ] = new double**[pairs_chi];
				this->Table_SR[1][ Tpair ] = new double**[pairs_chi];
				break;

				case 2:
				this->Table_TR[0][ Tpair ] = new double**[pairs_chi];
				this->Table_TR[1][ Tpair ] = new double**[pairs_chi];
				break;

				case 3:
				for(int ir = 0; ir < 2; ir++)
				{
					this->Table_SR[ir][ Tpair ] = new double**[pairs_chi];
					this->Table_TR[ir][ Tpair ] = new double**[pairs_chi];
				}
				break;
			}

			const double Rcut1 = ORB.Phi[T1].getRcut();
			const double Rcut2 = ORB.Phi[T2].getRcut();
			assert(Rcut1>0.0 && Rcut1<100);
			assert(Rcut2>0.0 && Rcut2<100);

			const int rmesh = this->get_rmesh( Rcut1, Rcut2);
			assert( rmesh < this->Rmesh );
			
			for (int L1 = 0; L1 < Lmax1 + 1; L1++)
			{
				for (int N1 = 0; N1 < ORB.Phi[T1].getNchi(L1); N1++)
				{
					for (int L2 = 0; L2 < Lmax2 + 1; L2 ++)
					{
						for (int N2 = 0; N2 < ORB.Phi[T2].getNchi(L2); N2++)
						{		
							// get the second index.
							const int Opair = this->OV_Opair(Tpair,L1,L2,N1,N2);
							
							// init 3rd dimension
							switch( job0 )
							{
								case 1:
								this->Table_SR[0][ Tpair ][ Opair ] = new double *[L2plus1];
								this->Table_SR[1][ Tpair ][ Opair ] = new double *[L2plus1];
								break;

								case 2:
								this->Table_TR[0][ Tpair ][ Opair ] = new double *[L2plus1];
								this->Table_TR[1][ Tpair ][ Opair ] = new double *[L2plus1];

								case 3:
								for(int ir = 0; ir < 2; ir++)
								{
									this->Table_SR[ir][ Tpair ][ Opair ] = new double *[L2plus1];
									this->Table_TR[ir][ Tpair ][ Opair ] = new double *[L2plus1];
								}
							}


							//L=|L1-L2|,|L1-L2|+2,...,L1+L2
							const int SL = abs(L1-L2);
							const int AL = L1+L2;
							
							for (int L=0; L < L2plus1 ; L++)
							{
								//Allocation
								switch ( job0 )
								{
									case 1:
									Table_SR[0][Tpair][Opair][L] = new double[rmesh];
									Table_SR[1][Tpair][Opair][L] = new double[rmesh];

									Memory::record("Make_Overlap_Table","Table_SR",
									2*OV_nTpairs*pairs_chi*rmesh,"double");
									break;

									case 2:
									Table_TR[0][Tpair][Opair][L] = new double[rmesh];
									Table_TR[1][Tpair][Opair][L] = new double[rmesh];

									Memory::record("Make_Overlap_Table","Table_TR",
									2*OV_nTpairs*pairs_chi*rmesh,"double");
									break;

									case 3:
									Table_SR[0][Tpair][Opair][L] = new double[rmesh];
									Table_SR[1][Tpair][Opair][L] = new double[rmesh];
									Table_TR[0][Tpair][Opair][L] = new double[rmesh];
									Table_TR[1][Tpair][Opair][L] = new double[rmesh];

									Memory::record("Make_Overlap_Table","Table_SR&TR",
									2*2*OV_nTpairs*pairs_chi*rmesh,"double");
									break;
								}
									
								//for those L whose Gaunt Coefficients = 0, we
								//assign every element in Table_SR or Table_TR as zero
								if ((L > AL) || (L < SL) || ((L-SL) % 2 == 1)) 
								{
									switch ( job0 )
									{
										case 1:
										ZEROS (Table_SR[0][Tpair][Opair][L], rmesh);
										ZEROS (Table_SR[1][Tpair][Opair][L], rmesh);
										break;

										case 2:
										ZEROS (Table_TR[0][Tpair][Opair][L], rmesh);
										ZEROS (Table_TR[1][Tpair][Opair][L], rmesh);
										break;

										case 3:
										ZEROS (Table_SR[0][Tpair][Opair][L], rmesh);
										ZEROS (Table_SR[1][Tpair][Opair][L], rmesh);
										ZEROS (Table_TR[0][Tpair][Opair][L], rmesh);
										ZEROS (Table_TR[1][Tpair][Opair][L], rmesh);
										break;
									}

									continue;
								}
								
								switch( job0 )
								{
									case 1:
									{
										this->cal_ST_Phi12_R(1,L, 
												ORB.Phi[T1].PhiLN(L1,N1),
												ORB.Phi[T2].PhiLN(L2,N2),
												rmesh,
												Table_SR[0][Tpair][Opair][L],
												Table_SR[1][Tpair][Opair][L]);
										break;
									}
									case 2:
									{

										this->cal_ST_Phi12_R(2,L, 
												ORB.Phi[T1].PhiLN(L1,N1),
												ORB.Phi[T2].PhiLN(L2,N2),
												rmesh,
												Table_TR[0][Tpair][Opair][L],
												Table_TR[1][Tpair][Opair][L]);
										break;
									}
									case 3:
									{	
										this->cal_ST_Phi12_R(1,L, 
												ORB.Phi[T1].PhiLN(L1,N1),
												ORB.Phi[T2].PhiLN(L2,N2),
												rmesh,
												Table_SR[0][Tpair][Opair][L],
												Table_SR[1][Tpair][Opair][L]);

										this->cal_ST_Phi12_R(2,L, 
												ORB.Phi[T1].PhiLN(L1,N1),
												ORB.Phi[T2].PhiLN(L2,N2),
												rmesh,
												Table_TR[0][Tpair][Opair][L],
												Table_TR[1][Tpair][Opair][L]);
										break;
									}
								}
							}//end m
						}
					}//end jl
				}
			}// end il
		}// end jt
	}// end it

	switch( job0 )
	{
		case 1:
		destroy_sr = true;
		break;

		case 2:
		destroy_tr = true;
		break;

		case 3:
		destroy_sr = true;
		destroy_tr = true;
		break;
	}
		
	timer::tick("Make_Overlap_Table", "init_Table",'D');
	return;
}

void Make_Overlap_Table::init_Table_Beta(void)
{
	TITLE("Make_Overlap_Table", "init_Table_Beta");
	timer::tick("Make_Overlap_Table", "init_Table_Beta",'D');

	// (1) allocate 1st dimension ( overlap, derivative)
	this->Table_NR = new double****[2];
	// (2) allocate 2nd dimension ( overlap, derivative)
	this->Table_NR[0] = new double*** [this->NL_nTpairs];
	this->Table_NR[1] = new double*** [this->NL_nTpairs];

	
	// <1Phi|2Beta> 
	for (int T1 = 0;  T1 < ntype ; T1++) // type 1 is orbital
	{
		for (int T2 = 0 ; T2 < ntype ; T2++)// type 2 is non-local projector
		{
			// Tpair: type pair.
			const int Tpair=this->NL_Tpair(T1,T2);
			const int Lmax1 = ORB.Phi[T1].getLmax();
			
			const int NBeta = ORB.nproj[T2];
			
			//-------------------------------------------------------------
			// how many <psi|beta_l>
			// here we count all possible psi with (L,N) index for type T1.
			//-------------------------------------------------------------
			const int pairs_chi = ORB.Phi[T1].getTotal_nchi() * NBeta;

			// CAUTION!!!
			// no matter nchi = 0 or NBeta = 0,
			// means the Tpair in this table is never used!
			if(pairs_chi == 0)continue;

			// init 2nd dimension
			this->Table_NR[0][Tpair] = new double** [ pairs_chi ];
			this->Table_NR[1][Tpair] = new double** [ pairs_chi ];

            const int T12_2Lplus1 = this->NL_L2plus1(T1,T2);

			const double Rcut1 = ORB.Phi[T1].getRcut();
			for (int L1 = 0; L1 < Lmax1 + 1; L1++)
            {
                for (int N1 = 0; N1 < ORB.Phi[T1].getNchi(L1); N1++)
				{
					
					// number of projectors.
					for (int nb = 0; nb < NBeta; nb ++)
					{
						const int L2 = ORB.Beta[T2].getL_Beta(nb);

						const double Rcut2 = ORB.Beta[T2].Proj[nb].getRcut();

						const int Opair = this->NL_Opair(Tpair,L1,N1,nb);
						assert( Opair < pairs_chi );

						// init 3rd dimension
						this->Table_NR[0][ Tpair ][ Opair ] = new double *[T12_2Lplus1];
						this->Table_NR[1][ Tpair ][ Opair ] = new double *[T12_2Lplus1];
						
						const int rmesh = this->get_rmesh( Rcut1, Rcut2);
						assert( rmesh < this->Rmesh );

						//not all L in T12_2Lplus1 would function
						const int SL = abs(L1-L2);
						const int AL = L1+L2;
							
						for (int L=0; L < T12_2Lplus1 ; L++)
						{
							//Allocation
							this->Table_NR[0][Tpair][Opair][L] = new double[rmesh];
							this->Table_NR[1][Tpair][Opair][L] = new double[rmesh];

							Memory::record("Make_Overlap_Table","Table_NR",
							2*NL_nTpairs*pairs_chi*rmesh,"double");

							//for those L whose Gaunt Coefficients = 0, we
							//assign every element in Table_NR as zero
							if ((L > AL) || (L < SL) || ((L-SL) % 2 == 1)) 
							{
								ZEROS (Table_NR[0][Tpair][Opair][L], rmesh);
								ZEROS (Table_NR[1][Tpair][Opair][L], rmesh);
								
								continue;
							}

							assert(nb < ORB.nproj[T2]);	
							this->cal_VNL_PhiBeta_R(L,
                                ORB.Phi[T1].PhiLN(L1,N1),
                                ORB.Beta[T2].Proj[nb], // mohan update 2011-03-07
                                rmesh,
								this->Table_NR[0][Tpair][Opair][L],
								this->Table_NR[1][Tpair][Opair][L]);
						}// end T12_2Lplus1
					}// end L2
				}// end N1
			}// end L1
		}// end T2
	}// end T1
	destroy_nr = true;


//	OUT(ofs_running,"allocate non-local potential matrix","Done");
	timer::tick("Make_Overlap_Table", "init_Table_Beta",'D');
	return;
}


void Make_Overlap_Table::Destroy_Table(void)
{
	if(!destroy_sr && !destroy_tr) return;
	
	const int ntype = ORB.get_ntype();
	int dim1 = 0;
	for (int ir = 0; ir < 2; ir++)
	{
    	for (int T1 = 0; T1 < ntype; T1++)
		{
			// Notice !! T2 start from T1
			// means that T2 >= T1
    	    for (int T2 = T1; T2 < ntype; T2++)
        	{
				const int Lmax1 = ORB.Phi[T1].getLmax();
				const int Lmax2 = ORB.Phi[T2].getLmax();
				const int lmax_now = std::max(Lmax1, Lmax2);
				const int pairs = ORB.Phi[T1].getTotal_nchi() * ORB.Phi[T2].getTotal_nchi();
				
				for (int dim2 = 0; dim2 < pairs; dim2++)
				{
					for (int L = 0; L < 2*lmax_now + 1; L++)
					{
						if(destroy_sr) delete [] Table_SR[ir][dim1][dim2][L];
						if(destroy_tr) delete [] Table_TR[ir][dim1][dim2][L];
                	}
                	if(destroy_sr) delete [] Table_SR[ir][dim1][dim2];
					if(destroy_tr) delete [] Table_TR[ir][dim1][dim2];
				}
            	if(destroy_sr) delete [] Table_SR[ir][dim1];
				if(destroy_tr) delete [] Table_TR[ir][dim1];
            	dim1++;

			}
        }

		dim1 = 0;
		if(destroy_sr) delete [] Table_SR[ir];
		if(destroy_tr) delete [] Table_TR[ir];
	}

	if(destroy_sr) delete[] Table_SR;
	if(destroy_tr) delete[] Table_TR;

	return;
}

void Make_Overlap_Table::Destroy_Table_Beta(void)
{
	if(!destroy_nr) return;

	const int ntype = ORB.get_ntype();
	for(int ir = 0; ir < 2; ir ++)
	{
		for(int T1=0; T1<ntype; T1++)
		{
			for(int T2=0; T2<ntype; T2++)
			{
				const int Tpair = this->NL_Tpair(T1,T2); 
				const int L2plus1 = this->NL_L2plus1(T1,T2);
				const int pairs = ORB.Phi[T1].getTotal_nchi() * ORB.nproj[T2]; 

				// mohan fix bug 2011-03-30
				if(pairs ==0) continue;
				for(int dim2=0; dim2<pairs; dim2++)
				{
					for(int L=0; L<L2plus1; L++)
					{
						delete[] Table_NR[ir][Tpair][dim2][L];
					}
					delete[] Table_NR[ir][Tpair][dim2];
				}
				delete[] Table_NR[ir][Tpair];
			}
		}
		delete[] Table_NR[ir];
	}
	delete[] Table_NR;
	return;
}


void Make_Overlap_Table::init_OV_Tpair(void)
{
	TITLE("Make_Overlap_Table","init_OV_Tpair");
    assert(ntype>0);

    this->OV_nTpairs = this->ntype * (this->ntype + 1) / 2;
    this->OV_Tpair.create(ntype, ntype);
	this->OV_L2plus1.create(ntype, ntype); // mohan fix bug 2011-03-14

//	OUT(ofs_running,"Number of type pairs",this->OV_nTpairs);

    int index = 0;
    for (int T1 = 0;  T1 < ntype ; T1++)
    {
		// Notice !! T2 start from T1
		// means that T2 >= T1
        for (int T2 = T1 ; T2 < ntype ; T2++)
        {
			// (1) pairs about atom types
			//liaochen modify 2010/8/4
			//index for T1>T2 is also needed
            this->OV_Tpair(T2, T1) = index;						
			this->OV_Tpair(T1, T2) = this->OV_Tpair(T2, T1);
            
			++index;
			// (2) pairs about lmax
			this->OV_L2plus1(T1,T2) = max(ORB.Phi[T1].getLmax(), ORB.Phi[T2].getLmax() )*2+1;
			this->OV_L2plus1(T2,T1) = this->OV_L2plus1(T1,T2);
        }
    }
    return;
}


void Make_Overlap_Table::init_NL_Tpair(void)
{
	TITLE("Make_Overlap_Table","init_NL_index");
	assert(ntype>0);
	this->NL_nTpairs = this->ntype * this->ntype;	
	this->NL_Tpair.create( this->ntype, this->ntype);
	this->NL_L2plus1.create( this->ntype, this->ntype); // mohan fix bug 2011-03-14

//	OUT(ofs_running,"Number of Nonlocal Pairs",NL_nTpairs);

	int index = 0;
	for (int T1 = 0;  T1 < ntype ; T1++)
	{
		for (int T0 = 0 ; T0 < ntype ; T0++)
		{
			 this->NL_Tpair(T1,T0) = index;
			 ++index;

			 // the pair < psi | beta >
			 // be careful! This is not a symmetry matrix.
			 this->NL_L2plus1(T1,T0) = std::max(ORB.Phi[T1].getLmax(), ORB.Beta[T0].getLmax() )*2+1;
			 
			 // there are special situations:
			 // for example, two H atom without projector.
			 // if we use s orbital, 
			 // Phi.getLmax = 0,
			 // Beta.getLmax < 0, 
			 // so the value is 1.
			 // however, there are no projectors.
			 if(NL_L2plus1(T1,T0) <= 0)
			 {
				WARNING_QUIT("Make_Overlap_Table::init_paris_nonlocal_type","NL_L2plus1<=0");
			 }
		}
	}
	return;
}


void Make_Overlap_Table::init_OV_Opair(void)
{
    const int lmax = ORB.get_lmax(); 
    const int nchimax = ORB.get_nchimax();
	assert(lmax+1 > 0);
	assert(nchimax > 0);
	assert(OV_nTpairs > 0);

    this->OV_Opair.create(OV_nTpairs, lmax+1, lmax+1, nchimax, nchimax);

    for(int T1=0; T1<ntype; T1++)
    {
		// Notice !! T2 start from T1
		// means that T2 >= T1
        for(int T2=T1; T2<ntype; T2++)
        {
			const int dim1 = this->OV_Tpair(T1,T2);
			int index=0;
            for(int L1=0; L1<ORB.Phi[T1].getLmax()+1; L1++)
            {
                for(int N1=0; N1<ORB.Phi[T1].getNchi(L1); N1++)
                {
                    for(int L2=0; L2<ORB.Phi[T2].getLmax()+1; L2++)
                    {
                        for(int N2=0; N2<ORB.Phi[T2].getNchi(L2); N2++)
                        {
                            this->OV_Opair(dim1, L1, L2, N1, N2) = index;
                            ++index;
                        }
                    }
                }
            }
        }
    }
    return;
}

void Make_Overlap_Table::init_NL_Opair(void)
{
	const int lmax = ORB.get_lmax();
	const int nchimax = ORB.get_nchimax();
	const int nprojmax = ORB.nprojmax;
	
	// may have bug if we use all H!
	if( nprojmax == 0)
	{
		WARNING("Make_Overlap_Table","nproj for nonlocal pseudopotetials are zero, it must be all H atoms");
		return;
	}
	assert( NL_nTpairs > 0);
	
	this->NL_Opair.create( this->NL_nTpairs, lmax+1, nchimax, nprojmax);
	
	// <1psi|2beta>
	// 1. orbital
	for(int T1=0; T1<ntype; T1++)
	{
		// 2. NL projector
		for(int T0=0; T0<ntype; T0++)
		{
			const int nlpair = this->NL_Tpair(T1, T0);
			int index = 0;
			for(int L1=0; L1<ORB.Phi[T1].getLmax()+1; L1++)
			{
				for(int N1=0; N1<ORB.Phi[T1].getNchi(L1); N1++)
				{
					// notice !! T0 must be Beta( Nonlocal projector)
					// mohan update 2011-03-07
					for(int ip=0; ip<ORB.nproj[T0]; ip++)
					{
						assert( nlpair < NL_nTpairs );
						assert( L1 < lmax+1 );
						assert( N1 < nchimax );
						assert( ip < nprojmax );
						this->NL_Opair(nlpair, L1, N1, ip) = index;
						++index;
					}
				}
			}
		}
	}

	return;
}

// Peize Lin update 2016-01-26
void Make_Overlap_Table::init_Lmax (const int orb_num, const int mode, int &Lmax_used, int &Lmax) const
{
	auto cal_Lmax_Phi = [&Lmax]()
	{
		//obtain maxL of all type
		const int ntype = ORB.get_ntype();
		for (int it = 0; it < ntype; it++)
		{
			const int Lmax_now = ORB.Phi[it].getLmax ();
			Lmax = max(Lmax, Lmax_now);
		}
	};

	auto cal_Lmax_Beta = [&Lmax]()
	{
		// fix bug.
		// mohan add the nonlocal part.
		// 2011-03-07
		const int ntype = ORB.get_ntype();
		for(int it=0; it< ntype; it++)
		{
			const int Lmax_now = ORB.Beta[it].getLmax();
			Lmax = max(Lmax, Lmax_now);
		}
	};

	Lmax = -1;
	
	switch( orb_num )
	{
		case 2:
			switch( mode )
			{
				case 1:			// used in <Phi|Phi> or <Beta|Phi>
					cal_Lmax_Phi();
					cal_Lmax_Beta();
					//use 2lmax+1 in dS
					Lmax_used = 2*Lmax + 1;
					break;
//				case 2:			// used in <jY|jY>
//					Lmax = max(Lmax, Exx_Abfs::Lmax);
//					Lmax_used = 2*Lmax + 1;
//					break;
				default:
					throw invalid_argument("Make_Overlap_Table::init_Lmax orb_num=2, mode error");
					break;
			}
			break;
		case 3:
			switch( mode )
			{
//				case 1:			// used in <jY|PhiPhi>
//					cal_Lmax_Phi();
//					Lmax = max(Lmax, Exx_Abfs::Lmax);
//					Lmax_used = 2*Lmax + 1;
//					Lmax_used += Exx_Abfs::Lmax;
//					break;
				default:
					throw invalid_argument("Make_Overlap_Table::init_Lmax orb_num=3, mode error");
					break;
			}
			break;
		case 4:
			switch( mode )
			{
				case 1:			// used in <PhiPhi|PhiPhi>
					cal_Lmax_Phi();
					Lmax_used = 2*( 2*Lmax + 1 );
					break;
				default:
					throw invalid_argument("Make_Overlap_Table::init_Lmax orb_num=4, mode error");
					break;
			}
			break;
		default:
			throw invalid_argument("Make_Overlap_Table::init_Lmax orb_num error");
			break;
	}

/*	if( mode==1 || mode==2)
	{
		//obtain maxL of all type

		const int ntype = ORB.get_ntype();
		for (int it = 0; it < ntype; it++)
		{
			const int Lmax_now = ORB.Phi[it].getLmax ();
			Lmax = max(Lmax, Lmax_now);
		}

		if(mode==1)
		{
			// fix bug.
			// mohan add the nonlocal part.
			// 2011-03-07
			for(int it=0; it< ntype; it++)
			{
				const int Lmax_now = ORB.Beta[it].getLmax();
				Lmax = max(Lmax, Lmax_now);
			}			
		}

		//use 2lmax+1 in dS
		Lmax_used = 2*Lmax + 1;	

		// Peize Lin add 2016-01-26
		if(mode==2)									
		{
			Lmax = max(Lmax, Use_Psi3_Center2::exx_ri_lmax);
			Lmax_used += Use_Psi3_Center2::exx_ri_lmax;
		}			
	}
	else if( mode==3 )
	{
		Lmax = Use_Psi3_Center2::exx_ri_lmax;
		Lmax_used = 2*Lmax + 1;
	}*/

	assert(Lmax_used >= 1);
}

// Peize Lin update 2016-01-26
void Make_Overlap_Table::init_Table_Spherical_Bessel (const int orb_num, const int mode, int &Lmax_used, int &Lmax)
{
	if( this->destroy_jlx )
	{
		WARNING_QUIT("Make_Overlap_Table::init_Table_Spherical_Bessel","jlx has been allocated!");
	}

	this->init_Lmax (orb_num,mode,Lmax_used,Lmax);		// Peize Lin add 2016-01-26

	// the allocation of L need to + 1,
	Sph_Bessel SB;
	this->jlx = new double**[Lmax_used+1];
	for (int l = 0; l < Lmax_used+1; l++)
	{
		this->jlx[l] = new double*[this->Rmesh];
		for (int ir = 0; ir < this->Rmesh; ir++)
		{
			this->jlx[l][ir] = new double[this->kmesh];
			SB.jlx(this->kmesh,this->kpoint,this->r[ir], l, this->jlx[l][ir]);
		}
	}

/*
// some data:
//L	x		Jl(x)old	Jl(x)web(correct)
//0	4		-0.189201	-0.18920062383
//1	11.7663	-0.0643896	-0.064389590588
//3	1.5048	0.028574	0.028573980746
//5	12.8544	-0.00829602	-0.0082960169277
//6	12.8544	-0.0776037	-0.077603690549
//7	12.8544	-0.0560009	-0.070186679825
//7	12		-0.0198184	-0.024838740722
    int lll;
    int ir;
    int ik;
    cout << " INPUT L:  " ; cin >> lll;
    cout << " INPUT ir: " ; cin >> ir;
    cout << " INPUT ik: " ; cin >> ik;
    double kr = r[ir] * kpoint[ik];
    cout <<  " L=" << lll << " kr=" << kr << " J=" << jlx[lll][ir][ik] << endl;
    goto once_again;
*/


	//memory-free flag
	this->destroy_jlx = true;

	OUT(ofs_running,"lmax used to generate Jlq",Lmax_used);
//	OUT(ofs_running,"kmesh",kmesh);
//	OUT(ofs_running,"Rmesh",Rmesh);
	Memory::record ("Make_Overlap_Table", "Jl(x)", (Lmax_used+1) * this->kmesh * this->Rmesh, "double");
}

void Make_Overlap_Table::Destroy_Table_Spherical_Bessel (const int& Lmax_used)
{
	if (!this->destroy_jlx) return;
	
	TITLE("Make_Overlap_Table","Destroy_Table_Spherical_Bessel");
//	OUT(ofs_running,"Lmax_used+1",Lmax_used+1);
//	OUT(ofs_running,"Rmesh",Rmesh);
//	OUT(ofs_running,"kmesh",kmesh);

	for (int l = 0; l < Lmax_used+1; l++)
	{
		for (int ir = 0; ir < this->Rmesh; ir++)
		{
			delete [] this->jlx[l][ir];
		}
		delete [] this->jlx[l];
	}

	delete[] this->jlx;
	destroy_jlx = false;
	
	return;
}
	
	
