#include "record_adj.h"
#include "../src_pw/global.h"

Record_adj::Record_adj(){}
Record_adj::~Record_adj(){}

void Record_adj::delete_grid(void)
{
	for(int i=0; i<na_proc; i++)
	{
		// how many 'numerical orbital' adjacents
		// for each atom in this process.
		for(int j=0; j< na_each[i]; j++)
		{
			delete[] info[i][j];
		}
		delete[] info[i];
	}
	delete[] info;
	delete[] na_each;
}



//--------------------------------------------
// This will record the orbitals according to
// HPSEPS's 2D block division.
//--------------------------------------------
#include "LCAO_nnr.h"
void Record_adj::for_2d(void)
{
	TITLE("Record_adj","for_2d");
	timer::tick("Record_adj","for_2d");

	assert(ucell.nat>0);

	// (1) find the adjacent atoms of atom[T1,I1];
	Vector3<double> tau1, tau2, dtau;
	Vector3<double> dtau1, dtau2, tau0;

	this->na_proc = ucell.nat;

	// number of adjacents for each atom.	
	this->na_each = new int[na_proc];
	ZEROS(na_each, na_proc);
	int iat = 0;
	int irr = 0;

	
//	cout << " in for_2d" << endl;

	for (int T1 = 0; T1 < ucell.ntype; ++T1)
	{
		Atom* atom1 = &ucell.atoms[T1];
		for (int I1 = 0; I1 < atom1->na; ++I1)
		{
			tau1 = atom1->tau[I1];
			//GridD.Find_atom( tau1 );
			GridD.Find_atom(ucell,  tau1 ,T1, I1);
			const int start1 = ucell.itiaiw2iwt(T1, I1, 0);

			// (2) search among all adjacent atoms.
			for (int ad = 0; ad < GridD.getAdjacentNum()+1; ++ad)
			{
				const int T2 = GridD.getType(ad);
				const int I2 = GridD.getNatom(ad);
				const int start2 = ucell.itiaiw2iwt(T2, I2, 0);
				tau2 = GridD.getAdjacentTau(ad);
				dtau = tau2 - tau1;
				double distance = dtau.norm() * ucell.lat0;
				double rcut = ORB.Phi[T1].getRcut() + ORB.Phi[T2].getRcut();

				bool is_adj = false;
				if( distance < rcut) is_adj = true;
				else if( distance >= rcut)
				{
					for (int ad0 = 0; ad0 < GridD.getAdjacentNum()+1; ++ad0)
					{
						const int T0 = GridD.getType(ad0);
						//const int I0 = GridD.getNatom(ad0);
						//const int iat0 = ucell.itia2iat(T0, I0);
						//const int start0 = ucell.itiaiw2iwt(T0, I0, 0);

						tau0 = GridD.getAdjacentTau(ad0);
						dtau1 = tau0 - tau1;
						double distance1 = dtau1.norm() * ucell.lat0;
						double rcut1 = ORB.Phi[T1].getRcut() + ORB.Beta[T0].get_rcut_max();

						dtau2 = tau0 - tau2;
						double distance2 = dtau2.norm() * ucell.lat0;
						double rcut2 = ORB.Phi[T2].getRcut() + ORB.Beta[T0].get_rcut_max();

						if( distance1 < rcut1 && distance2 < rcut2 )
						{
							is_adj = true;
							break;
						} // dis1, dis2
					}
				}



				if(is_adj)
				{
					++na_each[iat];

					for(int ii=0; ii<atom1->nw * GlobalV::NPOL; ++ii)
					{
						// the index of orbitals in this processor
						const int iw1_all = start1 + ii;
						const int mu = ParaO.trace_loc_row[iw1_all];
						if(mu<0)continue;

						for(int jj=0; jj<ucell.atoms[T2].nw * GlobalV::NPOL; ++jj)
						{
							const int iw2_all = start2 + jj;
							const int nu = ParaO.trace_loc_col[iw2_all];
							if(nu<0)continue;
							
							++irr;
						}
					}
				}//end is_adj
			}//end ad
			++iat;
		}//end I1
	}//end T1

 	//xiaohui add "GlobalV::OUT_LEVEL", 2015-09-16
	if(GlobalV::OUT_LEVEL != "m") OUT(GlobalV::ofs_running,"irr",irr);
	if(GlobalV::OUT_LEVEL != "m") OUT(GlobalV::ofs_running,"LNNR.nnr",LNNR.nnr);



	//------------------------------------------------
	// info will identify each atom in each unitcell.
	//------------------------------------------------
	this->info = new int**[na_proc];
	for(int i=0; i<na_proc; i++)
	{
//		GlobalV::ofs_running << " atom" << setw(5) << i << setw(10) << na_each[i] << endl;
		if( na_each[i] > 0)
		{
			info[i] = new int*[ na_each[i] ];
			for(int j=0; j<na_each[i]; j++)
			{
				// (Rx, Ry, Rz, T, I)
				info[i][j] = new int[5];
				ZEROS(info[i][j],5);
			}
		}
	}

	iat = 0;
	for (int T1 = 0; T1 < ucell.ntype; ++T1)
	{
		Atom* atom1 = &ucell.atoms[T1];
		for (int I1 = 0; I1 < atom1->na; ++I1)
		{
			tau1 = atom1->tau[I1];
			//GridD.Find_atom( tau1 );
			GridD.Find_atom(ucell,  tau1 ,T1, I1);

			// (2) search among all adjacent atoms.
			int cb = 0;
			for (int ad = 0; ad < GridD.getAdjacentNum()+1; ++ad)
			{
				const int T2 = GridD.getType(ad);
				const int I2 = GridD.getNatom(ad);
				tau2 = GridD.getAdjacentTau(ad);
				dtau = tau2 - tau1;
				double distance = dtau.norm() * ucell.lat0;
				double rcut = ORB.Phi[T1].getRcut() + ORB.Phi[T2].getRcut();


				bool is_adj = false;
				if(distance < rcut) is_adj = true;
				else if(distance >= rcut) 
				{
					for (int ad0 = 0; ad0 < GridD.getAdjacentNum()+1; ++ad0)
					{
						const int T0 = GridD.getType(ad0);
						//const int I0 = GridD.getNatom(ad0);
						//const int iat0 = ucell.itia2iat(T0, I0);
						//const int start0 = ucell.itiaiw2iwt(T0, I0, 0);

						tau0 = GridD.getAdjacentTau(ad0);
						dtau1 = tau0 - tau1;
						double distance1 = dtau1.norm() * ucell.lat0;
						double rcut1 = ORB.Phi[T1].getRcut() + ORB.Beta[T0].get_rcut_max();

						dtau2 = tau0 - tau2;
						double distance2 = dtau2.norm() * ucell.lat0;
						double rcut2 = ORB.Phi[T2].getRcut() + ORB.Beta[T0].get_rcut_max();

						if( distance1 < rcut1 && distance2 < rcut2 )
						{
							is_adj = true;
							break;
						} // dis1, dis2
					}
				}

				if(is_adj)
				{
					info[iat][cb][0] = GridD.getBox(ad).x; 
					info[iat][cb][1] = GridD.getBox(ad).y; 
					info[iat][cb][2] = GridD.getBox(ad).z; 
					info[iat][cb][3] = T2;
					info[iat][cb][4] = I2;
					++cb;
				}
			}//end ad
//			GlobalV::ofs_running << " nadj = " << cb << endl;
			++iat;
		}//end I1
	}//end T1
	timer::tick("Record_adj","for_2d");

	return;
}


//--------------------------------------------
// This will record the orbitals according to
// grid division (cut along z direction) 
//--------------------------------------------
void Record_adj::for_grid(const Grid_Technique &gt)
{
    TITLE("Record_adj","for_grid");
	timer::tick("Record_adj","for_grid");
	
	Vector3<double> tau1, tau2, dtau;
	Vector3<double> tau0, dtau1, dtau2;
	
	this->na_proc = 0;
	for(int T1=0; T1<ucell.ntype; ++T1)
	{
		for(int I1=0; I1<ucell.atoms[T1].na; ++I1)
		{
			const int iat = ucell.itia2iat(T1,I1);
			if(gt.in_this_processor[iat])
			{
				++na_proc;
			}
		}
	}

	// number of adjacents for each atom.	
	this->na_each = new int[na_proc];
	ZEROS(na_each, na_proc);

	int ca = 0;
	for(int T1=0; T1<ucell.ntype; ++T1)
	{
		Atom* atom1 = &ucell.atoms[T1];
		for(int I1=0; I1<atom1->na; ++I1)
		{
			const int iat = ucell.itia2iat(T1,I1);
			// key in this function
			if(gt.in_this_processor[iat])
			{	
				tau1 = atom1->tau[I1];
				//GridD.Find_atom(tau1);
				GridD.Find_atom(ucell, tau1, T1, I1);
				for (int ad = 0; ad < GridD.getAdjacentNum()+1; ad++)
				{
					const int T2 = GridD.getType(ad);
					const int I2 = GridD.getNatom(ad);
					const int iat2 = ucell.itia2iat(T2, I2);
					if(gt.in_this_processor[iat2])
					{
						//Atom* atom2 = &ucell.atoms[T2];
						tau2 = GridD.getAdjacentTau(ad);
						dtau = tau2 - tau1;
						double distance = dtau.norm() * ucell.lat0;
						double rcut = ORB.Phi[T1].getRcut() + ORB.Phi[T2].getRcut();


						bool is_adj = false;
						if(distance < rcut) is_adj = true;
						/*
						else if(distance >= rcut)
						{
                            for (int ad0 = 0; ad0 < GridD.getAdjacentNum()+1; ++ad0)
                            {
                                const int T0 = GridD.getType(ad0);
                                const int I0 = GridD.getNatom(ad0);
                                const int iat0 = ucell.itia2iat(T0, I0);
                                const int start0 = ucell.itiaiw2iwt(T0, I0, 0);

                                tau0 = GridD.getAdjacentTau(ad0);
                                dtau1 = tau0 - tau1;
                                dtau2 = tau0 - tau2;

                                double distance1 = dtau1.norm() * ucell.lat0;
                                double distance2 = dtau2.norm() * ucell.lat0;

                                double rcut1 = ORB.Phi[T1].getRcut() + ORB.Beta[T0].get_rcut_max();
                                double rcut2 = ORB.Phi[T2].getRcut() + ORB.Beta[T0].get_rcut_max();

                                if( distance1 < rcut1 && distance2 < rcut2 )
                                {
                                   	is_adj = true; 
                                    break;
								} // dis1, dis2
							}
						}
						*/

						// check the distance
						if(is_adj)
						{
							++na_each[ca];
						}
					}//end judge 2
				}//end ad
				++ca;
			}//end judge 1
		}//end I1
	}//end T1


	this->info = new int**[na_proc];
	for(int i=0; i<na_proc; i++)
	{
		assert(na_each[i]>0);
		info[i] = new int*[ na_each[i] ];
		for(int j=0; j<na_each[i]; j++)
		{
			// (Rx, Ry, Rz, T, I)
			info[i][j] = new int[5];
			ZEROS(info[i][j],5);
		}
	}

	ca = 0;
	for(int T1=0; T1<ucell.ntype; ++T1)
	{
		Atom* atom1 = &ucell.atoms[T1];
		for(int I1=0; I1 < atom1->na; ++I1)
		{
			const int iat = ucell.itia2iat(T1,I1);

			// key of this function
			if(gt.in_this_processor[iat])
			{
				tau1 = atom1->tau[I1];
				//GridD.Find_atom(tau1);
				GridD.Find_atom(ucell, tau1, T1, I1);

				int cb = 0;
				for (int ad = 0; ad < GridD.getAdjacentNum()+1; ad++)
				{
					const int T2 = GridD.getType(ad);
					const int I2 = GridD.getNatom(ad);
					const int iat2 = ucell.itia2iat(T2, I2);

					// key of this function
					if(gt.in_this_processor[iat2])
					{
						//Atom* atom2 = &ucell.atoms[T2];
						tau2 = GridD.getAdjacentTau(ad);
						dtau = tau2 - tau1;
						double distance = dtau.norm() * ucell.lat0;
						double rcut = ORB.Phi[T1].getRcut() + ORB.Phi[T2].getRcut();

						// check the distance
						if(distance < rcut)
						{
							info[ca][cb][0] = GridD.getBox(ad).x; 
							info[ca][cb][1] = GridD.getBox(ad).y; 
							info[ca][cb][2] = GridD.getBox(ad).z; 
							info[ca][cb][3] = T2;
							info[ca][cb][4] = I2;
							++cb;
						}
						/*
                        else if(distance >= rcut)
                        {
                            for (int ad0 = 0; ad0 < GridD.getAdjacentNum()+1; ++ad0)
                            {
                                const int T0 = GridD.getType(ad0);
                                const int I0 = GridD.getNatom(ad0);
                                const int iat0 = ucell.itia2iat(T0, I0);
                                const int start0 = ucell.itiaiw2iwt(T0, I0, 0);

                                tau0 = GridD.getAdjacentTau(ad0);
                                dtau1 = tau0 - tau1;
                                dtau2 = tau0 - tau2;

                                double distance1 = dtau1.norm() * ucell.lat0;
                                double distance2 = dtau2.norm() * ucell.lat0;

                                double rcut1 = ORB.Phi[T1].getRcut() + ORB.Beta[T0].get_rcut_max();
                                double rcut2 = ORB.Phi[T2].getRcut() + ORB.Beta[T0].get_rcut_max();

                                if( distance1 < rcut1 && distance2 < rcut2 )
                                {
									info[ca][cb][0] = GridD.getBox(ad).x; 
									info[ca][cb][1] = GridD.getBox(ad).y; 
									info[ca][cb][2] = GridD.getBox(ad).z; 
									info[ca][cb][3] = T2;
									info[ca][cb][4] = I2;
									++cb;
									break;
								} // dis1, dis2
							}
                        }
						*/
					}
				}// end ad

				assert(cb == na_each[ca]);
				++ca;
			}
		}
	}
	assert(ca==na_proc);
	timer::tick("Record_adj","for_grid");

//	cout << " after for_grid" << endl;
	return;
}

