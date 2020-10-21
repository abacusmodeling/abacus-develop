#include "bfield.h"
#include "../src_pw/tools.h"
#include "../src_pw/global.h"
#include "global_fp.h"

Bfield bfid;

Bfield::Bfield()
{
	tesla_x = 0.0;
	tesla_y = 0.0;
	tesla_z = 0.0;
	Rydberg_x=0.0;
	Rydberg_y=0.0;
	Rydberg_z=0.0; //sun zhiyuan updated at 2011-12-26//
	Gauge_Origin_x=0.0;
	Gauge_Origin_y=0.0;
	Gauge_Origin_z=0.0; //sun zhiyuan updated at 2011-12-26//
	allocate_tab = false;
	nlocal=0;
	fac=3.0083e-6;
	bohr_mag=1.414213562;
	Ready_A_of_Atom=0;
	fac_of_phase = 1.414213562/1.0;//'1.0' need to be updated to 'c'. sun zhiyuan added at 2011-12-31//
	c=1.0;//137.036*2.0;
}

Bfield::~Bfield()
{	

	for(int iat=0; iat<ucell.nat; iat++)
	{
		//		delete[] A_of_Atom[iat];
	}
	//	delete[] A_of_Atom;
}

void Bfield::check()
{
//=========================================================================
//	double A[3]={1.0,1.0,1.0};
//	double r[3]={1.0,1.0,1.0};
//	complex<double> chec = this->cal_phase(A,r,1);
//	cout<<"======"<<chec.real()<<"+ i"<<chec.imag()<<"==========="<<endl;
//=========================================================================
//	cout<<"Gauge_Origin_x=============="<<Gauge_Origin_x<<endl;
//	cout<<"Gauge_Origin_y=============="<<Gauge_Origin_y<<endl;
//	cout<<"Gauge_Origin_z=============="<<Gauge_Origin_z<<endl;
//=========================================================================
//	complex<double> trial = complex<double>(0,1)*1.0;//complex<double>(1,1);
//	cout<<"==========trial==="<<trial.real()<<"+i"<<trial.imag()<<endl;
}

void Bfield::convert()    //sun zhiyuan updated at 2011-12-26//
{
	this->Rydberg_x=tesla_x*this->fac;
	this->Rydberg_y=tesla_y*this->fac;
	this->Rydberg_z=tesla_z*this->fac;
}


void Bfield::cal_A(double A[3],double r[3])   //sun zhiyuan updated at 2011-12-26//
{
	double B[3],rtmp[3],DoublA[3]={0,0,0};
	B[0]=Rydberg_x;
	B[1]=Rydberg_y;
	B[2]=Rydberg_z;
	rtmp[0]=r[0]-Gauge_Origin_x;
	rtmp[1]=r[1]-Gauge_Origin_y;
	rtmp[2]=r[2]-Gauge_Origin_z;
	a_x_b(B,rtmp,DoublA);
	A[0]=DoublA[0]*0.500;
	A[1]=DoublA[1]*0.500;
	A[2]=DoublA[2]*0.500;
	return;
}
void Bfield::cal_A_of_Atom(void)  //Zhiyuan add at 2011-12-24. Former version:cal_RxB//
{
	TITLE("Bfield","cal_A_of_Atom");

	cout << " begin to prepare for bfield calculations using local orbitals." << endl;

	this->A_of_Atom = new double*[ucell.nat];
	for(int iat=0; iat<ucell.nat; iat++)
	{
		this->A_of_Atom[iat] = new double[3];
		ZEROS(A_of_Atom[iat], 3);
	}	

	double R[3];
	int iat=0;
	for(int it=0; it<ucell.ntype; it++)
	{
		Atom* atom = &ucell.atoms[it];
		for(int ia=0; ia<ucell.atoms[it].na; ia++)
		{
			R[0] = atom->tau[ia].x*ucell.lat0;
			R[1] = atom->tau[ia].y*ucell.lat0;
			R[2] = atom->tau[ia].z*ucell.lat0;
			this->cal_A(A_of_Atom[iat],R);
			++iat;
		}
	}
    this->Ready_A_of_Atom=1;
	return;
}

// exp[(+ or -) i (sqrt(2)/c)*(A dot r)],Zhiyuan add 2011-12-31//
complex<double> Bfield::cal_phase(double A[3],double r[3], int sign) 
{                                                	
	assert(sign == 1 || sign == -1);
	double A_dot_r = A[0]*r[0]+A[1]*r[1]+A[2]*r[2]; //this->a_dot_b(A,r);
	double whole = A_dot_r*this->fac_of_phase;
	if(sign == -1)
	{
		whole = - whole;
	}
	double e1 = std::cos(whole);
	double e2 = std::sin(whole);
	return complex<double>(e1, e2);
}

void Bfield::a_x_b(const double* a, const double* b, double *c)
{
	// write a program to calculate a * b, return c
	c[0] = a[1] * b[2] - a[2] * b[1];
	c[1] = a[2] * b[0] - a[0] * b[2];
	c[2] = a[0] * b[1] - a[1] * b[0];
}

double Bfield::a_dot_b(const double* a, const double* b)
{
	return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

double Bfield::a_dot_b(const Vector3<double> &a, const double* b)
{
	return a.x*b[0]+a.y*b[1]+a.z*b[2];
}

void Bfield::get_nlocal(void)
{
	this->trace=new int [NLOCAL];
	for(int i=0;i<NLOCAL;i++)
		trace[i]=-1;
	Memory::record("Bfield","trace",NLOCAL,"int");


	nlocal=0;
	for(int it=0;it<ucell.ntype;it++)
		for(int ia=0;ia<ucell.atoms[it].na;ia++)
		{
			const int start = ucell.itiaiw2iwt(it, ia, 0);
			for(int iw=0;iw<ucell.atoms[it].nw;iw++)
			{
				const int iw_all=start+iw;
				if(ParaO.trace_loc_row[iw_all]>=0|| ParaO.trace_loc_col[iw_all]>=0)
				{
					trace[iw_all]=nlocal;
					nlocal++;
				}
			}
		}
	cout<<"============basis in this processor"<<nlocal<<"==========="<<endl;
	assert(nlocal<=NLOCAL);
	//  WARNING_QUIT("QUIT IN bfield::get_nlocal","as expected");
	return;
}



void Bfield::make_table(void)
{
	TITLE("Bfield","make_table");
	this->get_nlocal();
	int nkb = ORB.nkb;

	cout << " the dimension of table is (atom, l1, n1, m1) * (atom, l, m) = " 
		<< NLOCAL << " * " << nkb << endl; 
	assert( nkb >= 0);  //in case of all Hydrogen atom, nkb=0//

	this->Tab = new complex<double>*[nlocal];
	if(nkb>0)
	{
		for(int i=0; i<nlocal; i++)
		{
			Tab[i] = new complex<double>[nkb];
			ZEROS(Tab[i], nkb);
		}
		Memory::record("Bfield","table",nlocal*nkb,"double");

		UHM.GG.cal_vnl_B(Tab);
	}

	this->allocate_tab = true;


/*
	cout << " Grid integration: table of <phi_lnm | beta_lm >" << endl;
	for(int i=0; i<nlocal; i++)
	{
		cout << " ";
		for(int j=0; j<nkb; j++)
		{
			if( norm(Tab[i][j]) > 0.001 )
			{
				cout << setw(15) << Tab[i][j].real();
				ofs_running<<setw(15)<<Tab[i][j].real();  //sun zhiyuan for test//
			}
			else
			{
				cout << setw(15) << "0";
				ofs_running<<setw(15)<<"0";  //sun zhiyuan for test//
			}
		}
		cout << endl;
	}
	*/

	return;
}

/*
   void Bfield::calculate_NL_B(void)
   {
   TITLE("Bfield","calculate_NL_B");
   assert(allocate_tab);
   for(int i=0;i<NLOCAL;i++)
   {
   if(ParaO.trace_loc_row[i]<0)  continue;
   for(int j=0;j<NLOCAL;j++)
   {
   if(ParaO.trace_loc_col[j]<0)  continue;
   int indexi=trace[i];
   int indexj=trace[j];
   complex<double> nlm=ZERO;
   this->snap_psibeta(nlm,indexi,indexj);
   LM.set_HSk(i,j,nlm,'N');
   }
   }
   }
 */	



void Bfield::calculate_NL_B(void)
{
	TITLE("Bfield","calculate_NL_B");

	assert(allocate_tab);

	for (int T1 = 0; T1 < ucell.ntype; T1++)
	{
		const Atom* atom1 = &ucell.atoms[T1];
		for (int I1 =0; I1< atom1->na; I1++)
		{
			//GridD.Find_atom( atom1->tau[I1] );
			GridD.Find_atom( atom1->tau[I1] ,T1, I1);
			const int iat1 = ucell.itia2iat(T1, I1);
			const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
			const Vector3<double> tau1 = atom1->tau[I1];

			for (int ad2=0; ad2<ucell.nat ; ad2++)
			{
				const int T2 = ucell.iat2it[ad2];
				const Atom* atom2 = &ucell.atoms[T2];

				const int I2 = ucell.iat2ia[ad2];
				const int iat2 = ucell.itia2iat(T2, I2);
				const int start2 = ucell.itiaiw2iwt(T2, I2, 0);
				const Vector3<double> tau2 = ucell.atoms[T2].tau[I2];

				// < psi1 | all projectors | psi2 >
				for (int j=0; j<atom1->nw; j++)
				{
					const int iw1_all = start1 + j;
					const int mu = ParaO.trace_loc_row[iw1_all];
					if(mu < 0)continue; 

					// fix a serious bug: atom2[T2] -> atom2
					// mohan 2010-12-20
					for (int k=0; k<atom2->nw; k++)
					{
						const int iw2_all = start2 + k;
						const int nu = ParaO.trace_loc_col[iw2_all];						
						if(nu < 0)continue;
						//(3) run over all projectors in nonlocal pseudopotential.
						for (int ad0=0; ad0 < ucell.nat; ad0++)
						{
							const int T0 = ucell.iat2it[ad0];
							const int I0 = ucell.iat2ia[ad0];

							// mohan add 2010-12-19
							//							if( ORB.nproj[T0] == 0) continue; 

							const int start0 = ucell.itiaiw2iwt(T0, I0, 0);
							const Vector3<double> tau0 = ucell.atoms[T0].tau[I0];

							complex<double> nlm=ZERO;
							this->snap_psibeta(
									nlm, 0, 
									tau1, T1, iw1_all,
									tau2, T2, iw2_all, 
									tau0, T0, I0);
							LM.set_HSk(iw1_all,iw2_all,nlm,'N');
						} // ad0
					}// k
				} // j 
				//----------------------------------------------------------------------------------
			} // ad2
		} // I1
	} // T1

	delete [] trace;

	return;
}

void Bfield::snap_psibeta(
		complex<double> &nlm,
		const int& job,
		const Vector3<double> &R1,
		const int &T1,
		const int &iw1_all,
		const Vector3<double> &R2,
		const int &T2,
		const int &iw2_all,
		const Vector3<double> &R0,// The projector.
		const int &T0,
		const int &I0) const
{
	const int nproj = ORB.nproj[T0];
	Numerical_Orbital::set_position(R1, R2);

	// (1) get distance between R1 and R2 (a.u.)
	// judge if there exist overlap
	double distance = Numerical_Orbital::get_distance()*ucell.lat0;

	const double Rcut1 = ORB.Phi[T1].getRcut();
	const double Rcut2 = ORB.Phi[T2].getRcut();

	nlm = ZERO;

	if( distance > (Rcut1 + Rcut2) ) return;

	int indexi=trace[iw1_all];
	int indexj=trace[iw2_all]; //sun zhiyuan add//

	int ib_start = ORB.itiaib2ib_all(T0,I0,0);
	for(int nb=0; nb<nproj; nb++)
	{
		const int L0 = ORB.Beta[T0].Proj[nb].getL();
		for(int m0=0; m0<2*L0+1; m0++)
		{
			nlm += conj( Tab[indexi][ib_start] ) * Tab[indexj][ib_start]* ORB.Beta[T0].getCoefficient_D(L0, L0);
			++ib_start;
		}
	}
	return;
}


void Bfield::add_zeeman()
{
    assert(NSPIN==2);
	// E = B dot Mag
	// sqrt(B)
    double Bmode=sqrt(bfid.Rydberg_x*bfid.Rydberg_x+
				bfid.Rydberg_y*bfid.Rydberg_y+
				bfid.Rydberg_z*bfid.Rydberg_z);

	// magnetism 
    double bohr_mag=bfid.bohr_mag;
    for (int is = 0;is < NSPIN;is++)
    {
        if(is == 0)
        {
            for (int i = 0;i < pw.nrxx;i++)
            {
               	pot.vrs(is, i) += bohr_mag*Bmode;
            }
        }
        if(is == 1)
        {
            for (int i = 0;i < pw.nrxx;i++)
            {
              	pot.vrs(is, i) -= bohr_mag*Bmode;
            }
        }
    }
	return;
}
