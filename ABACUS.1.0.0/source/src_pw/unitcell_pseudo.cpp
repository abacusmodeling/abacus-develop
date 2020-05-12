//==========================================================
// AUTHOR : mohan
// DATE : 2008-11-10
//==========================================================
#include "unitcell_pseudo.h"
#include "global.h"
#include <cstring>		// Peize Lin fix bug about strcmp 2016-08-02

UnitCell_pseudo::UnitCell_pseudo()
{
	if(test_pseudo_cell) TITLE("unitcell_pseudo","Constructor");
	set_atom_flag = false;
}

UnitCell_pseudo::~UnitCell_pseudo()
{
	if(set_atom_flag)
	{
		delete[] atoms;
	}
}

//==============================================================
//Calculate various lattice related quantities for given latvec
//==============================================================
void UnitCell_pseudo::setup_cell(
		const string &s_pseudopot_dir ,  
		const string &fn , ofstream &log)
{
	TITLE("UnitCell_pseudo","setup_cell");	
	// (1) init mag (global class)
	assert(ntype>0);
	delete[] mag.start_magnetization;
	mag.start_magnetization = new double[this->ntype];

	// (2) init *Atom class array.
	this->atoms = new Atom[this->ntype]; // atom species.
	this->set_atom_flag = true;


	bool ok = true;
	bool ok2 = true;

	// (3) read in atom information
	if(MY_RANK == 0)
	{
		// open "atom_unitcell" file.
		ifstream ifa(fn.c_str(), ios::in);
		if (!ifa)
		{
			ofs_warning << fn;
			ok = false;
		}

		if(ok)
		{

			ofs_running << "\n\n\n\n";
			ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
			ofs_running << " |                                                                    |" << endl;
			ofs_running << " | Reading atom information in unitcell:                              |" << endl;
			ofs_running << " | From the input file and the structure file we know the number of   |" << endl;
			ofs_running << " | different elments in this unitcell, then we list the detail        |" << endl;
			ofs_running << " | information for each element, especially the zeta and polar atomic |" << endl;
			ofs_running << " | orbital number for each element. The total atom number is counted. |" << endl;
			ofs_running << " | We calculate the nearest atom distance for each atom and show the  |" << endl;
			ofs_running << " | Cartesian and Direct coordinates for each atom. We list the file   |" << endl;
			ofs_running << " | address for atomic orbitals and nonlocal projectors. The volume    |" << endl;
			ofs_running << " | and the lattice vectors in real and reciprocal space is also shown.|" << endl;
			ofs_running << " |                                                                    |" << endl;
			ofs_running << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
			ofs_running << "\n\n\n\n";



			ofs_running << " READING UNITCELL INFORMATION" << endl;
			//========================
			// call read_atom_species
			//========================
			this->read_atom_species(ifa);

			//==========================
			// call read_atom_positions
			//==========================
			ok2 = this->read_atom_positions(ifa);

			if(ok2)
			{
				for(int i=0;i<this->ntype;i++)
				{
					Global_File::make_dir_atom( this->atoms[i].label );
				}
			}
		}
	}
#ifdef __MPI
	Parallel_Common::bcast_bool(ok);
	Parallel_Common::bcast_bool(ok2);
	if(NONCOLIN) Parallel_Common::bcast_bool(DOMAG);
	if(NONCOLIN) Parallel_Common::bcast_bool(DOMAG_Z);
#endif
	if(!ok)
	{
		WARNING_QUIT("UnitCell_pseudo::setup_cell","Can not find the file containing atom positions.!");
	}
	if(!ok2)
	{
		WARNING_QUIT("UnitCell_pseudo::setup_cell","Something wrong during read_atom_positions.");
	}

#ifdef __MPI
	this->bcast_unitcell_pseudo();

	// mohan add 2010-09-29
	ORB.bcast_files();
#endif
	
	//========================================================
	// Calculate unit cell volume
	// the reason to calculate volume here is 
	// Firstly, latvec must be read in.
	//========================================================
	assert(lat0 > 0.0);
	this->omega = abs( latvec.Det() ) * this->lat0 * lat0 * lat0 ;
	if(this->omega<=0)
	{
		WARNING_QUIT("setup_cell","omega <= 0 .");
	}
	else
	{
		ofs_running << endl;
		OUT(ofs_running,"Volume (Bohr^3)", this->omega);
		OUT(ofs_running,"Volume (A^3)", this->omega * pow(BOHR_TO_A, 3));
	}
		
	//==========================================================
	// Calculate recip. lattice vectors and dot products
	// latvec has the unit of lat0, but G has the unit 2Pi/lat0
	//==========================================================
	this->GT = latvec.Inverse();
	this->G  = GT.Transpose();
	this->GGT = G * GT;
	this->invGGT = GGT.Inverse();
        //LiuXh add 20180515
        this->GT0 = latvec.Inverse();
        this->G0  = GT.Transpose();
        this->GGT0 = G * GT;
        this->invGGT0 = GGT.Inverse();

	ofs_running << endl;
	out.printM3(ofs_running,"Lattice vectors: (Cartesian coordinate: in unit of a_0)",latvec); 
	out.printM3(ofs_running,"Reciprocal vectors: (Cartesian coordinate: in unit of 2 pi/a_0)",G);
//	OUT(ofs_running,"lattice center x",latcenter.x);
//	OUT(ofs_running,"lattice center y",latcenter.y);
//	OUT(ofs_running,"lattice center z",latcenter.z);

	// read in non-local pseudopotential and ouput the projectors.

	ofs_running << "\n\n\n\n";
	ofs_running << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;
	ofs_running << " |                                                                    |" << endl;
	ofs_running << " | Reading pseudopotentials files:                                    |" << endl;
	ofs_running << " | The pseudopotential file is in UPF format. The 'NC' indicates that |" << endl;
	ofs_running << " | the type of pseudopotential is 'norm conserving'. Functional of    |" << endl;
	ofs_running << " | exchange and correlation is decided by 4 given parameters in UPF   |" << endl;
	ofs_running << " | file.  We also read in the 'core correction' if there exists.      |" << endl;
	ofs_running << " | Also we can read the valence electrons number and the maximal      |" << endl;
	ofs_running << " | angular momentum used in this pseudopotential. We also read in the |" << endl;
	ofs_running << " | trail wave function, trail atomic density and local-pseudopotential|" << endl;
	ofs_running << " | on logrithmic grid. The non-local pseudopotential projector is also|" << endl;
	ofs_running << " | read in if there is any.                                           |" << endl;
	ofs_running << " |                                                                    |" << endl;
	ofs_running << " <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
	ofs_running << "\n\n\n\n";


	this->read_pseudopot(s_pseudopot_dir);
	
	if(MY_RANK == 0) 
	{
		for(int it=0; it<this->ntype; it++)
		{
			Atom* atom = &atoms[it];
			stringstream ss;
			ss << global_out_dir << atom->label 
				<< "/" << atom->label
				<< ".NONLOCAL";
			ofstream ofs(ss.str().c_str());

			ofs << "<HEADER>" << endl;
			ofs << setw(10) << atom->label << "\t" << "label" << endl;
			ofs << setw(10) << atom->pp_type << "\t" << "pseudopotential type" << endl;
			ofs << setw(10) << atom->lmax << "\t" << "lmax" << endl;
			ofs << "</HEADER>" << endl;

			ofs << "\n<DIJ>" << endl;
			ofs << setw(10) << atom->nbeta << "\t" << "nummber of projectors." << endl;
			
			for(int ib=0; ib<atom->nbeta; ib++)
			{
				for(int ib2=0; ib2<atom->nbeta; ib2++)
				{
					ofs<<setw(10) << atom->lll[ib] 
						<< " " << atom->lll[ib2]
						<< " " << atom->dion(ib,ib2)<<endl;
				}
			}
			ofs << "</DIJ>" << endl;
			
			for(int i=0; i<atom->nbeta; i++)
			{
				ofs << "<PP_BETA>" << endl;
				ofs << setw(10) << i << "\t" << "the index of projectors." <<endl;
				ofs << setw(10) << atom->lll[i] << "\t" << "the angular momentum." <<endl;

				// mohan add
				// only keep the nonzero part.
				int cut_mesh = atom->mesh; 
				for(int j=atom->mesh-1; j>=0; --j)
				{
					if( abs( atom->betar(i,j) ) > 1.0e-10 )
					{
						cut_mesh = j; 
						break;
					}
				}
				if(cut_mesh %2 == 0) ++cut_mesh;

				ofs << setw(10) << cut_mesh << "\t" << "the number of mesh points." << endl;


				for(int j=0; j<cut_mesh; ++j)
				{
					ofs << setw(15) << atom->r[j]
						<< setw(15) << atom->betar(i, j)
						<< setw(15) << atom->rab[j] << endl;
				}
				ofs << "</PP_BETA>" << endl;
			}

			ofs.close();
		}
	}

#ifdef __MPI
	this->bcast_unitcell_pseudo2();
#endif	

	for(int it=0; it<ntype; it++)
	{
		for(int j=0; j<4; j++)
		{
			if(atoms[0].dft[j]!=atoms[it].dft[j])
			{
				ofs_warning << "\n type " << atoms[0].label << " functional is " 
				<< atoms[0].dft[0] << " " << atoms[0].dft[1] << " "
				<< atoms[0].dft[2] << " " << atoms[0].dft[3];
				
				ofs_warning << "\n type " << atoms[it].label << " functional is " 
				<< atoms[it].dft[0] << " " << atoms[it].dft[1] << " "
				<< atoms[it].dft[2] << " " << atoms[it].dft[3] << endl;
				
				WARNING_QUIT("setup_cell","All DFT functional must consistent.");
			}
		}
	}

	// mohan add 2010-09-06
	// because the number of element type
	// will easily be ignored, so here
	// I warn the user again for each type.
	for(int it=0; it<ntype; it++)
	{
		xcf.which_dft(atoms[it].dft);
	}

	this->cal_nelec();
	this->cal_natomwfc();
	this->cal_nwfc();
	this->cal_meshx();

//	stringstream ss;
//	ss << global_out_dir << "unitcell_pp.log";
//	print_unitcell_pseudo( ss.str() );
	return;
}

void UnitCell_pseudo::read_atom_species(ifstream &ifa)
{
	TITLE("UnitCell_pseudo","read_atom_species");

	this->atom_mass  = new double[ntype]; //atom masses
	this->atom_label = new string[ntype]; //atom labels
	this->pseudo_fn  = new string[ntype]; //file name of pseudopotential

	string word;
	//==========================================
	// read in information of each type of atom
	//==========================================
	if( SCAN_BEGIN(ifa, "ATOMIC_SPECIES") )
	{
		OUT(ofs_running,"ntype",ntype);
		for (int i = 0;i < ntype;i++)
		{
			ifa >> atom_label[i] >> atom_mass[i];
			
			stringstream ss;
			ss << "atom label for species " << i+1;
			OUT(ofs_running,ss.str(),atom_label[i]);	
			READ_VALUE(ifa, pseudo_fn[i]);
			if(test_pseudo_cell==2) 
			{
				ofs_running << "\n" << setw(6) << atom_label[i] 
						<< setw(12) << atom_mass[i] 
						<< setw(18) << pseudo_fn[i];
			}
		}
	}

#ifdef __FP
	if( SCAN_BEGIN(ifa, "NUMERICAL_ORBITAL") )
	{
		ORB.read_in_flag = true;
		for(int i=0; i<ntype; i++)
		{
			string ofile;

			//-----------------------------------
			// Turn off the read in NONLOCAL file
			// function since 2013-08-02 by mohan
			//-----------------------------------
			//string nfile;

			ifa >> ofile;
			//-----------------------------------
			// Turn off the read in NONLOCAL file
			// function since 2013-08-02 by mohan
			//-----------------------------------
			//READ_VALUE(ifa, nfile);
			
			ORB.orbital_file.push_back(ofile);

			//-----------------------------------
			// Turn off the read in NONLOCAL file
			// function since 2013-08-02 by mohan
			//-----------------------------------
			//ORB.nonlocal_file.push_back(nfile);

//			ofs_running << " For atom type " << i + 1 << endl;
//		    ofs_running << " Read in numerical orbitals from file " << ofile << endl;
//		    ofs_running << " Read in nonlocal projectors from file " << nfile << endl;
			
		}
	}	
#endif

	//==========================
	// read in lattice constant
	//==========================
	if( SCAN_BEGIN(ifa, "LATTICE_CONSTANT") )
	{
		READ_VALUE(ifa, lat0);
		if(lat0<=0.0)
		{
			WARNING_QUIT("read_atom_species","lat0<=0.0");
		}
		lat0_angstrom = lat0 * 0.529177 ;
		OUT(ofs_running,"lattice constant (Bohr)",lat0);
		OUT(ofs_running,"lattice constant (Angstrom)",lat0_angstrom);
		this->tpiba  = TWO_PI / lat0;
		this->tpiba2 = tpiba * tpiba;
	}

	//===========================
	// Read in latticies vector
	//===========================
	if( SCAN_BEGIN(ifa, "LATTICE_VECTORS") )
	{
		// Reading lattice vectors. notice
		// here that only one cpu read these
		// parameters.
		ifa >> latvec.e11 >> latvec.e12;
		READ_VALUE(ifa, latvec.e13);
		ifa >> latvec.e21 >> latvec.e22;
		READ_VALUE(ifa, latvec.e23);
		ifa >> latvec.e31 >> latvec.e32;
		READ_VALUE(ifa, latvec.e33);

		// lattice vectors in another form.
		a1.x = latvec.e11;
		a1.y = latvec.e12;
		a1.z = latvec.e13;

		a2.x = latvec.e21;
		a2.y = latvec.e22;
		a2.z = latvec.e23;

		a3.x = latvec.e31;
		a3.y = latvec.e32;
		a3.z = latvec.e33;
		
	}
	return;
}

// Read atomic positions
// return 1: no problem.
// return 0: some problems.
bool UnitCell_pseudo::read_atom_positions(ifstream &ifpos)
{
	TITLE("UnitCell_pseudo","read_atom_positions");

	bool use_xyz = false;
        /* LiuXh 20171109
	if( SCAN_BEGIN(ifpos, "USE_XYZ") )
	{
		use_xyz = true;
	}
        */

	if( SCAN_BEGIN(ifpos, "ATOMIC_POSITIONS"))
	{
		READ_VALUE( ifpos, Coordinate);
		if(Coordinate != "Cartesian" 
			&& Coordinate != "Direct" 
			&& Coordinate != "Cartesian_angstrom"
			&& Coordinate != "Cartesian_au"
			&& Coordinate != "Cartesian_angstrom_center_xy"
			&& Coordinate != "Cartesian_angstrom_center_xz"
			&& Coordinate != "Cartesian_angstrom_center_yz"
			&& Coordinate != "Cartesian_angstrom_center_xyz"
			)
		{
			WARNING("read_atom_position","Cartesian or Direct?");
			ofs_warning << " There are several options for you:" << endl;
			ofs_warning << " Direct" << endl;
			ofs_warning << " Cartesian_angstrom" << endl;
			ofs_warning << " Cartesian_au" << endl;
			ofs_warning << " Cartesian_angstrom_center_xy" << endl;
			ofs_warning << " Cartesian_angstrom_center_xz" << endl;
			ofs_warning << " Cartesian_angstrom_center_yz" << endl;
			ofs_warning << " Cartesian_angstrom_center_xyz" << endl;
			return 0;
		}

		Vector3<double> v;
		Vector3<int> mv;
		int na = 0;
		this->nat = 0;

		//======================================
		// calculate total number of atoms
		// and adjust the order of atom species
		//======================================
		assert(ntype>0);
		for (int it = 0;it < ntype; it++)
		{
			ofs_running << "\n READING ATOM TYPE " << it+1 << endl;
			
			//=======================================
			// (1) read in atom label
			// start magnetization
			//=======================================
			READ_VALUE(ifpos, atoms[it].label);
			bool found = false;
			for(int it2=0; it2<ntype; it2++)
			{
				if( this->atoms[it].label == this->atom_label[it] )
				found = true;
			}
			if(!found)
			{
				ofs_warning << " Label read from ATOMIC_POSITIONS is " << this->atoms[it].label << endl; 
				ofs_warning << " Lable from ATOMIC_SPECIES is " << this->atom_label[it] << endl;
				return 0;
			}
			READ_VALUE(ifpos, mag.start_magnetization[it] );

			OUT(ofs_running, "atom label",atoms[it].label);

			if(NSPIN==4)//added by zhengdy-soc
			{
				if(NONCOLIN && LSPINORB){
					if(fabs(mag.start_magnetization[it])>1e-6) DOMAG_Z = true;
				}
				if(DOMAG)
				{
					soc.m_loc[it].x = mag.start_magnetization[it] *
							sin(soc.angle1[it]) * cos(soc.angle2[it]);
					soc.m_loc[it].y = mag.start_magnetization[it] *
							sin(soc.angle1[it]) * sin(soc.angle2[it]);
					soc.m_loc[it].z = mag.start_magnetization[it] *
							cos(soc.angle1[it]);
				}
				else if(DOMAG_Z)
				{
					soc.m_loc[it].x = 0;
					soc.m_loc[it].y = 0;
					soc.m_loc[it].z = mag.start_magnetization[it];
				}

				OUT(ofs_running, "noncollinear magnetization_x",soc.m_loc[it].x);
				OUT(ofs_running, "noncollinear magnetization_y",soc.m_loc[it].y);
				OUT(ofs_running, "noncollinear magnetization_z",soc.m_loc[it].z);

				ZEROS(soc.ux ,3);
			}
			else if(NSPIN==2)
			{
				soc.m_loc[it].x = mag.start_magnetization[it];
				OUT(ofs_running, "start magnetization",mag.start_magnetization[it]);
			}
			else if(NSPIN==1)
			{
				OUT(ofs_running, "start magnetization","FALSE");
			}

			//===========================================
			// (2) read in numerical orbital information
			// int atoms[it].nwl
			// int* atoms[it].l_nchi;
			//===========================================
			/*READ_VALUE(ifpos, this->atoms[it].nwl);
			assert(this->atoms[it].nwl<10);
			OUT(ofs_running,"L max for local orbitals",atoms[it].nwl);

			delete[] this->atoms[it].l_nchi;
			this->atoms[it].l_nchi = new int[ this->atoms[it].nwl+1];
			this->atoms[it].nw = 0;
			for(int L=0; L<atoms[it].nwl+1; L++)
			{
				READ_VALUE(ifpos, this->atoms[it].l_nchi[L]);

				// calculate the number of local basis(3D)
				this->atoms[it].nw += (2*L + 1) * this->atoms[it].l_nchi[L];
				
				stringstream ss;
				ss << "L=" << L << ", number of zeta";
				OUT(ofs_running,ss.str(),atoms[it].l_nchi[L]);
			}*/
                    if (BASIS_TYPE == "lcao" || BASIS_TYPE == "lcao_in_pw")
                    {    
                        ifstream ifs(ORB.orbital_file[it].c_str(), ios::in);  // pengfei 2014-10-13
                        if (!ifs)
                        {
                            cout << " Can't find the ORBITAL file." << endl;
                        }

                        char word[80];
                        delete[] this->atoms[it].l_nchi;
                        this->atoms[it].l_nchi = new int[ this->atoms[it].nwl+1];
                        this->atoms[it].nw = 0;
                        int L =0;

                        while (ifs.good())
                        {
                             ifs >> word;
                             if (strcmp("Lmax", word) == 0)
                             {
                                 READ_VALUE(ifs, this->atoms[it].nwl);
                             }
                             assert(this->atoms[it].nwl<10);

                             if (strcmp("Cutoff(a.u.)", word) == 0)         // pengfei Li 16-2-29
                             {
                                 READ_VALUE(ifs, this->atoms[it].Rcut);
                             }

                             //cout << "atoms[it].nwl = "<<atoms[it].nwl <<endl;
                             if (strcmp("Sorbital-->", word) == 0)
                             {
                                 READ_VALUE(ifs, this->atoms[it].l_nchi[L]);
                                this->atoms[it].nw += (2*L + 1) * this->atoms[it].l_nchi[L];
                                stringstream ss;
                                ss << "L=" << L << ", number of zeta";
                                OUT(ofs_running,ss.str(),atoms[it].l_nchi[L]);
                                L++;
                             }
                             if (strcmp("Porbital-->", word) == 0)
                             {
                                 READ_VALUE(ifs, this->atoms[it].l_nchi[L]);
                                this->atoms[it].nw += (2*L + 1) * this->atoms[it].l_nchi[L];
                                stringstream ss;
                                ss << "L=" << L << ", number of zeta";
                                OUT(ofs_running,ss.str(),atoms[it].l_nchi[L]);
                                L++;
                             }
                             if (strcmp("Dorbital-->", word) == 0)
                             {
                                 READ_VALUE(ifs, this->atoms[it].l_nchi[L]);
                                this->atoms[it].nw += (2*L + 1) * this->atoms[it].l_nchi[L];
                                stringstream ss;
                                ss << "L=" << L << ", number of zeta";
                                OUT(ofs_running,ss.str(),atoms[it].l_nchi[L]);
                                L++;
                             }
                             if (strcmp("Forbital-->", word) == 0)
                             {
                                 READ_VALUE(ifs, this->atoms[it].l_nchi[L]);
                                this->atoms[it].nw += (2*L + 1) * this->atoms[it].l_nchi[L];
                                stringstream ss;
                                ss << "L=" << L << ", number of zeta";
                                OUT(ofs_running,ss.str(),atoms[it].l_nchi[L]);
                                L++;
                             }
                             if (strcmp("Gorbital-->", word) == 0)
                             {
                                 READ_VALUE(ifs, this->atoms[it].l_nchi[L]);
                                this->atoms[it].nw += (2*L + 1) * this->atoms[it].l_nchi[L];
                                stringstream ss;
                                ss << "L=" << L << ", number of zeta";
                                OUT(ofs_running,ss.str(),atoms[it].l_nchi[L]);
                                L++;
                             }
                             //cout <<" atoms[it].nw = "<<atoms[it].nw<<endl;


                        }
                             ifs.close();
                    }
                    else
                    {
                        delete[] this->atoms[it].l_nchi;
                        this->atoms[it].l_nchi = new int[ this->atoms[it].nwl+1];
                        this->atoms[it].nw = 0;

                             this->atoms[it].nwl = 2;
                             //cout << INPUT.lmaxmax << endl;
                             if ( INPUT.lmaxmax != 2 )
                             {
                                 this->atoms[it].nwl = INPUT.lmaxmax;
                             }
                        for(int L=0; L<atoms[it].nwl+1; L++)
                        {
                                this->atoms[it].l_nchi[L] = 1;
                                // calculate the number of local basis(3D)
                                this->atoms[it].nw += (2*L + 1) * this->atoms[it].l_nchi[L];

                                stringstream ss;
                                ss << "L=" << L << ", number of zeta";
                                OUT(ofs_running,ss.str(),atoms[it].l_nchi[L]);
                        }
    
                         /*   this->atoms[it].l_nchi[0] = 1;
                                this->atoms[it].nw += (2*0 + 1) * this->atoms[it].l_nchi[0];
                             //   stringstream s0;
                             //   s0 << "L=" << 0 << ", number of zeta";
                             //   OUT(ofs_running,s0.str(),atoms[it].l_nchi[0]);

                             this->atoms[it].l_nchi[1] = 1;
                                this->atoms[it].nw += (2*1 + 1) * this->atoms[it].l_nchi[1];
                             //   stringstream s1;
                             //   s1 << "L=" << 1 << ", number of zeta";
                             //   OUT(ofs_running,s1.str(),atoms[it].l_nchi[1]);

                             this->atoms[it].l_nchi[2] = 1;
                                this->atoms[it].nw += (2*2 + 1) * this->atoms[it].l_nchi[2];*/
                             //   stringstream s2;
                             //   s2 << "L=" << 2 << ", number of zeta";
                             //  OUT(ofs_running,s2.str(),atoms[it].l_nchi[2]);*/
                         
                    }

			//OUT(ofs_running,"Total number of local orbitals",atoms[it].nw);

			//=========================
			// (3) read in atom number
			//=========================
			READ_VALUE(ifpos, na);
			this->atoms[it].na = na;

			OUT(ofs_running,"number of atom for this type",na);

			this->nat += na;
			if (na <= 0) 
			{
				WARNING("read_atom_positions"," atom number < 0.");
				return 0;
			}
			if (na > 0)
			{
       			delete[] atoms[it].tau;
				delete[] atoms[it].taud;
       			delete[] atoms[it].mbl;
				delete[] atoms[it].mag;
       			atoms[it].tau = new Vector3<double>[na];
       			atoms[it].taud = new Vector3<double>[na];
       			atoms[it].mbl = new Vector3<int>[na];
				atoms[it].mag = new double[na];
				atoms[it].mass = this->atom_mass[it]; //mohan add 2011-11-07 
				ZEROS(atoms[it].mag,na);
				for (int ia = 0;ia < na; ia++)
				{
					if(use_xyz)
					{
						string tmpid;
						ifpos >> tmpid;
						if(tmpid != atoms[it].label)
						{
							WARNING("read_atom_positions","atom label not match");	
							return 0;
						} 
						else
						{
							ifpos >> v.x >> v.y >> v.z;

							mv.x = true;
							mv.y = true;
							mv.z = true;
							
							string mags;
							std::getline( ifpos, mags );
							// change string to double.
							atoms[it].mag[ia] = std::atof(mags.c_str());
						}
					}
					else
					{
						ifpos >> v.x >> v.y >> v.z
							>> mv.x >> mv.y >> mv.z;
						string mags;
						std::getline( ifpos, mags );
						// change string to double.
						atoms[it].mag[ia] = std::atof(mags.c_str());
					}	

					if(Coordinate=="Direct")
					{
						// change v from direct to cartesian,
						// the unit is pw.lat0
						atoms[it].taud[ia] = v;
						atoms[it].tau[ia] = v * latvec;
					}
					else if(Coordinate=="Cartesian")
					{
						atoms[it].tau[ia] = v ;// in unit lat0
						//cout << " T=" << it << " I=" << ia << " tau=" << atoms[it].tau[ia].x << " " << 
						//atoms[it].tau[ia].y << " " << atoms[it].tau[ia].z << endl;
					}
					else if(Coordinate=="Cartesian_angstrom")
					{
						atoms[it].tau[ia] = v / 0.529177 / lat0;
					}	
					else if(Coordinate=="Cartesian_angstrom_center_xy")
					{
						// calculate lattice center 
						latcenter.x = (latvec.e11 + latvec.e21 + latvec.e31)/2.0;
						latcenter.y = (latvec.e12 + latvec.e22 + latvec.e32)/2.0;
						latcenter.z = 0.0;
						atoms[it].tau[ia] = v / 0.529177 / lat0 + latcenter; 
					}
					else if(Coordinate=="Cartesian_angstrom_center_xz")
					{
						// calculate lattice center 
						latcenter.x = (latvec.e11 + latvec.e21 + latvec.e31)/2.0;
						latcenter.y = 0.0; 
						latcenter.z = (latvec.e13 + latvec.e23 + latvec.e33)/2.0;	
						atoms[it].tau[ia] = v / 0.529177 / lat0 + latcenter; 
					}
					else if(Coordinate=="Cartesian_angstrom_center_yz")
					{
						// calculate lattice center 
						latcenter.x = 0.0; 
						latcenter.y = (latvec.e12 + latvec.e22 + latvec.e32)/2.0;
						latcenter.z = (latvec.e13 + latvec.e23 + latvec.e33)/2.0;	
						atoms[it].tau[ia] = v / 0.529177 / lat0 + latcenter; 
					}
					else if(Coordinate=="Cartesian_angstrom_center_xyz")
					{
						// calculate lattice center 
						latcenter.x = (latvec.e11 + latvec.e21 + latvec.e31)/2.0;
						latcenter.y = (latvec.e12 + latvec.e22 + latvec.e32)/2.0;
						latcenter.z = (latvec.e13 + latvec.e23 + latvec.e33)/2.0;	
						atoms[it].tau[ia] = v / 0.529177 / lat0 + latcenter; 
					}
					else if(Coordinate=="Cartesian_au")
					{
						atoms[it].tau[ia] = v / lat0;
					}

					if(Coordinate=="Cartesian" || 
						Coordinate=="Cartesian_angstrom" || 
						Coordinate=="Cartesian_angstrom_center_xy" || 
						Coordinate=="Cartesian_angstrom_center_xz" || 
						Coordinate=="Cartesian_angstrom_center_yz" || 
						Coordinate=="Cartesian_angstrom_center_xyz" || 
						Coordinate=="Cartesian_au")
					{
						double dx,dy,dz;
						Mathzone::Cartesian_to_Direct(atoms[it].tau[ia].x, atoms[it].tau[ia].y, atoms[it].tau[ia].z,
						latvec.e11, latvec.e12, latvec.e13,
						latvec.e21, latvec.e22, latvec.e23,
						latvec.e31, latvec.e32, latvec.e33,
						dx,dy,dz);
					
						atoms[it].taud[ia].x = dx;
						atoms[it].taud[ia].y = dy;
						atoms[it].taud[ia].z = dz;

					}
					
					atoms[it].mbl[ia] = mv;
				}//endj
			}// end na
		}//end for ntype
	}// end scan_begin


	ofs_running << endl;
	OUT(ofs_running,"TOTAL ATOM NUMBER",nat);

	// mohan add 2010-06-30	
	//xiaohui modify 2015-03-15, cancel outputfile "STRU_READIN.xyz"
	//this->print_cell_xyz("STRU_READIN.xyz");
	this->check_dtau();

	if ( this->check_tau() )
	{

	}
	else
	{
		return 0;
	}
	this->print_tau();
	//xiaohui modify 2015-03-15, cancel outputfile "STRU_READIN.xyz"
	//this->print_cell_xyz("STRU_READIN_ADJUST.xyz");
	this->print_cell_cif("STRU_READIN_ADJUST.cif");

	return 1;
}//end read_atom_positions

bool UnitCell_pseudo::check_tau(void)const
{
	TITLE("UnitCell_pseudo","check_tau");
	timer::tick("UnitCell_pseudo","check_tau");
	
	Vector3<double> diff = 0.0;
	double norm = 0.0;
	double tolerence_bohr = 1.0e-3;

	ofs_running << "\n Output nearest atom not considering periodic boundary condition" << endl;
	ofs_running << " " << setw(5) << "TYPE" << setw(6) << "INDEX" 
	<< setw(20) << "NEAREST(Bohr)" 
	<< setw(20) << "NEAREST(Angstrom)" << endl; 
	for(int T1=0; T1< this->ntype; T1++)
	{
		for(int I1=0; I1< this->atoms[T1].na; I1++)
		{	
			double shortest_norm = 10000.0; // a large number
			//int nearest_atom_type = 0;
			//int nearest_atom_index = 0;
			for(int T2=0; T2<this->ntype; T2++)
			{
				for(int I2=0; I2<this->atoms[T2].na; I2++)
				{
					if(T1==T2 && I1==I2)
					{
						shortest_norm = 0.0;
						//nearest_atom_type = T1;
						//nearest_atom_index = I2;
						// self atom
					}
					else
					{
						diff = atoms[T1].tau[I1] - atoms[T2].tau[I2];
						norm = diff.norm() * ucell.lat0;
						if( shortest_norm > norm )
						{
							shortest_norm = norm;
							//nearest_atom_type = T2;
							//nearest_atom_index = I2;
						}
						if( norm < tolerence_bohr ) // unit is Bohr
						{	
							ofs_warning << " two atoms are too close!" << endl;
							ofs_warning << " type:" << this->atoms[T1].label << " atom " << I1 + 1 << endl; 
							ofs_warning << " type:" << this->atoms[T2].label << " atom " << I2 + 1 << endl; 
							ofs_warning << " distance = " << norm << " Bohr" << endl;
							return 0;
						}
					}
				}
			}
			ofs_running << " " << setw(5) << atoms[T1].label << setw(6) << I1+1 
			<< setw(20) << shortest_norm  
			<< setw(20) << shortest_norm * BOHR_TO_A << endl;
		}
	}

	timer::tick("UnitCell_pseudo","check_tau");
	return 1;
}

void UnitCell_pseudo::print_stru_file(const string &fn, const int &type)const
{
	TITLE("UnitCell_pseudo","print_stru_file");
	
	if(MY_RANK!=0) return;

	ofstream ofs(fn.c_str());

	ofs << "ATOMIC_SPECIES" << endl;
	ofs << setprecision(12);

	for(int it=0; it<ntype; it++)
	{
                //modified by zhengdy 2015-07-24
		ofs << atom_label[it] << " " << atom_mass[it] << " " << pseudo_fn[it] << endl;
	}
	
#ifdef __FP
	//if(LOCAL_BASIS) xiaohui modify 2013-09-02 //mohan fix bug 2011-05-01
	if(BASIS_TYPE=="lcao" || BASIS_TYPE=="lcao_in_pw") //xiaohui add 2013-09-02. Attention...
	{	
		ofs << "\nNUMERICAL_ORBITAL" << endl;
		for(int it=0; it<ntype; it++)
		{
			//-----------------------------------
			// Turn off the read in NONLOCAL file
			// function since 2013-08-02 by mohan
			//-----------------------------------
//			ofs << ORB.orbital_file[it] << " " << ORB.nonlocal_file[it] << " #local_orbital; non-local projector" << endl;
			//modified by zhengdy 2015-07-24
                        ofs << ORB.orbital_file[it] << endl;
		}
	}
#endif

	ofs << "\nLATTICE_CONSTANT" << endl;
        //modified by zhengdy 2015-07-24
	ofs << lat0 << endl;

	ofs << "\nLATTICE_VECTORS" << endl;
	ofs << latvec.e11 << " " << latvec.e12 << " " << latvec.e13 << " #latvec1" << endl; 
	ofs << latvec.e21 << " " << latvec.e22 << " " << latvec.e23 << " #latvec2" << endl;
	ofs << latvec.e31 << " " << latvec.e32 << " " << latvec.e33 << " #latvec3" << endl;
	
	ofs << "\nATOMIC_POSITIONS" << endl;

	if(type==1)
	{
		ofs << "Cartesian" << endl;
		for(int it=0; it<ntype; it++)
		{
			ofs << endl;
			ofs << atoms[it].label << " #label" << endl;
			ofs << mag.start_magnetization[it] << " #magnetism" << endl;
			//2015-05-07, modify
			//ofs << atoms[it].nwl << " #max angular momentum" << endl;
			//xiaohui modify 2015-03-15
			//for(int l=0; l<=atoms[it].nwl; l++)
			//{
			//	ofs << atoms[it].l_nchi[l] << " #number of zeta for l=" << l << endl;
			//}
			ofs << atoms[it].na << " #number of atoms" << endl;
			for(int ia=0; ia<atoms[it].na; ia++)
			{
				ofs << atoms[it].tau[ia].x << " " << atoms[it].tau[ia].y << " " << atoms[it].tau[ia].z << " " 
					<< atoms[it].mbl[ia].x << " " << atoms[it].mbl[ia].y << " " << atoms[it].mbl[ia].z << endl;
			}
		}
	}
	else if(type==2)
	{
		ofs << "Direct" << endl;
		for(int it=0; it<ntype; it++)
		{
			ofs << endl;
			ofs << atoms[it].label << " #label" << endl;
			ofs << mag.start_magnetization[it] << " #magnetism" << endl;
			//ofs << atoms[it].nwl << " #max angular momentum" << endl;
			//xiaohui modify 2015-03-15
			//for(int l=0; l<=atoms[it].nwl; l++)
			//{
			//	ofs << atoms[it].l_nchi[l] << " #number of zeta for l=" << l << endl;
			//}
			ofs << atoms[it].na << " #number of atoms" << endl;
			for(int ia=0; ia<atoms[it].na; ia++)
			{
				ofs << atoms[it].taud[ia].x << " " << atoms[it].taud[ia].y << " " << atoms[it].taud[ia].z << " " 
					<< atoms[it].mbl[ia].x << " " << atoms[it].mbl[ia].y << " " << atoms[it].mbl[ia].z << endl;
			}
		}
	}

	ofs.close();

	return;
}


void UnitCell_pseudo::print_tau(void)const
{
    TITLE("UnitCell_pseudo","print_tau");
    if(Coordinate == "Cartesian" || Coordinate == "Cartesian_angstrom")
    {
        ofs_running << "\n CARTESIAN COORDINATES ( UNIT = " << ucell.lat0 << " Bohr )." << endl;
        ofs_running << setw(13) << " atom"
        //<< setw(20) << "x" 
        //<< setw(20) << "y" 
        //<< setw(20) << "z" 
        //<< " mag"
        << setw(20) << "x"
        << setw(20) << "y"
        << setw(20) << "z"
        << setw(20) << "mag"
        << endl;
        ofs_running << setprecision(12);

        int iat=0;
        for(int it=0; it<ucell.ntype; it++)
        {
            for (int ia = 0; ia < ucell.atoms[it].na; ia++)
            {
                stringstream ss;
                ss << "tauc_" << ucell.atoms[it].label << ia+1;

                ofs_running << " " << setw(12) << ss.str()
                //<< setw(20) << atoms[it].tau[ia].x 
                //<< setw(20) << atoms[it].tau[ia].y
                //<< setw(20) << atoms[it].tau[ia].z
                //<< " " << atoms[it].mag[ia]
                << setw(20) << atoms[it].tau[ia].x
                << setw(20) << atoms[it].tau[ia].y
                << setw(20) << atoms[it].tau[ia].z
                << setw(20) << atoms[it].mag[ia]
                << endl;

                ++iat;
            }
        }
    }

    if(Coordinate == "Direct")
    {
        ofs_running << "\n DIRECT COORDINATES" << endl;
        ofs_running << setw(13) << " atom"
        //<< setw(20) << "x"
        //<< setw(20) << "y"
        //<< setw(20) << "z"
        //<< " mag"
        << setw(20) << "x"
        << setw(20) << "y"
        << setw(20) << "z"
        << setw(20) << "mag"
        << endl;

        int iat=0;
        for(int it=0; it<ucell.ntype; it++)
        {
            for (int ia = 0; ia < ucell.atoms[it].na; ia++)
            {
                stringstream ss;
                ss << "taud_" << ucell.atoms[it].label << ia+1;

                ofs_running << " " << setw(12) << ss.str()
                //<< setw(20) << atoms[it].taud[ia].x
                //<< setw(20) << atoms[it].taud[ia].y
                //<< setw(20) << atoms[it].taud[ia].z
                //<< " " << atoms[it].mag[ia]
                << setw(20) << atoms[it].taud[ia].x
                << setw(20) << atoms[it].taud[ia].y
                << setw(20) << atoms[it].taud[ia].z
                << setw(20) << atoms[it].mag[ia]
                << endl;

                ++iat;
            }
        }
    }

	ofs_running << endl;
	return;
}	

//==========================================================
// Read pseudopotential according to the dir
//==========================================================
void UnitCell_pseudo::read_pseudopot(const string &pp_dir)
{
	if(test_pseudo_cell) TITLE("UnitCell_pseudo","read_pseudopot");
//----------------------------------------------------------
// EXPLAIN : setup reading log for pseudopot_upf
//----------------------------------------------------------
	stringstream ss;
	ss << global_out_dir << "atom_pseudo.log";
	
//	ofstream ofs;
	
//	if(MY_RANK==0)
//	{
//		ofs.open( ss.str().c_str(), ios::out);
//	}

//----------------------------------------------------------
// EXPLAIN : Read in the atomic pseudo potential
//----------------------------------------------------------
	string pp_address;
	for (int i = 0;i < ntype;i++)
	{
		Pseudopot_upf upf;
	
		// mohan update 2010-09-12	
		int error = 0;
		
		if(MY_RANK==0)
		{
			pp_address = pp_dir + this->pseudo_fn[i];
			//error = upf.read_pseudo_upf( pp_address ); xiaohui modify 2013-06-23
			error = upf.init_pseudo_reader( pp_address ); //xiaohui add 2013-06-23
		}

#ifdef __MPI
		Parallel_Common::bcast_int(error);
#endif

		if(error==1)
		{
			cout << " Pseudopotential directory now is : " << pp_address << endl;
			ofs_warning << " Pseudopotential directory now is : " << pp_address << endl;
			WARNING_QUIT("UnitCell_pseudo::read_pseudopot","Couldn't find pseudopotential file.");
		}
		else if(error==2)
		{
			WARNING_QUIT("UnitCell_pseudo::read_pseudopot","Something in pseudopotential not match.");
		}
		else if(error==3)
		{
			//xiaohui modify 2015-03-25
			//WARNING_QUIT("UnitCell_pseudo::read_pseudopot","Please check the reference states in pseudopotential .vwr file.\n Also the norm of the read in pseudo wave functions\n explicitly please check S, P and D channels.\n If the norm of the wave function is \n unreasonable large (should be near 1.0), MESIA would quit. \n The solution is to turn off the wave functions  \n and the corresponding non-local projectors together\n in .vwr pseudopotential file.");
			WARNING_QUIT("UnitCell_pseudo::read_pseudopot","Please check the reference states in pseudopotential .vwr file.\n Also the norm of the read in pseudo wave functions\n explicitly please check S, P and D channels.\n If the norm of the wave function is \n unreasonable large (should be near 1.0), ABACUS would quit. \n The solution is to turn off the wave functions  \n and the corresponding non-local projectors together\n in .vwr pseudopotential file.");
		}
//		OUT(ofs_running,"PP_ERRROR",error);

//xiaohui add 2015-03-24
#ifdef __MPI
		Parallel_Common::bcast_bool(upf.functional_error);
#endif
		//xiaohui add 2015-03-24
		if(upf.functional_error == 1)
		{
			WARNING_QUIT("Pseudopot_upf::read_pseudo_header","input xc functional does not match that in pseudopot file");
		}

		if(MY_RANK==0)
		{
//			upf.print_pseudo_upf( ofs );
			atoms[i].set_pseudo_us( upf );

			ofs_running << "\n Read in pseudopotential file is " << pseudo_fn[i] << endl;
			OUT(ofs_running,"pseudopotential type",atoms[i].pp_type);
			OUT(ofs_running,"functional Ex", atoms[i].dft[0]);
			OUT(ofs_running,"functional Ec", atoms[i].dft[1]);
			OUT(ofs_running,"functional GCEx", atoms[i].dft[2]);
			OUT(ofs_running,"functional GCEc", atoms[i].dft[3]);
			OUT(ofs_running,"nonlocal core correction", atoms[i].nlcc);
//			OUT(ofs_running,"spin orbital",atoms[i].has_so);
			OUT(ofs_running,"valence electrons", atoms[i].zv);
			OUT(ofs_running,"lmax", atoms[i].lmax);
			OUT(ofs_running,"number of zeta", atoms[i].nchi);
			OUT(ofs_running,"number of projectors", atoms[i].nbeta);
			for(int ib=0; ib<atoms[i].nbeta; ib++)
			{
				OUT(ofs_running,"L of projector", atoms[i].lll[ib]);
			}
//			OUT(ofs_running,"Grid Mesh Number", atoms[i].mesh);
		}
			
		//atoms[i].print_pseudo_us(ofs);
	}

//	if(MY_RANK==0)
//	{
//		ofs.close();
//	}
	return;
}


//=========================================================
// calculate total number of electrons (nelec) and default
// number of bands (NBANDS).
//=========================================================
void UnitCell_pseudo::cal_nelec(void)
{
	if(test_pseudo_cell) TITLE("UnitCell_pseudo","cal_nelec");
	//=======================================================
	// calculate the total number of electrons in the system
	// if nelec <>0; use input number (setup.f90)
	//=======================================================

	ofs_running << "\n SETUP THE ELECTRONS NUMBER" << endl;

	if (nelec == 0)
	{
		for (int it = 0; it < ntype;it++)
		{
			stringstream ss1, ss2;
			ss1 << "electron number of element " << atoms[it].label;
			const int nelec_it = atoms[it].zv * atoms[it].na;
			nelec += nelec_it;
			ss2 << "total electron number of element " << atoms[it].label; 
			
			OUT(ofs_running,ss1.str(),atoms[it].zv);
			OUT(ofs_running,ss2.str(),nelec_it);
		}
	}

//	OUT(ofs_running,"Total nelec",nelec);

	//=======================================
	// calculate number of bands (setup.f90)
	//=======================================
	double occupied_bands = static_cast<double>(nelec/DEGSPIN);	

	if( (occupied_bands - std::floor(occupied_bands)) > 0.0 )
	{
		occupied_bands = std::floor(occupied_bands) + 1.0; //mohan fix 2012-04-16
	}

	OUT(ofs_running,"occupied bands",occupied_bands);
	
	// mohan add 2010-09-04
        //cout << "nbands(ucell) = " <<NBANDS <<endl;
	if(NBANDS==occupied_bands)
	{
		if( Occupy::gauss() || Occupy::tetra() )
		{
			WARNING_QUIT("UnitCell_pseudo::cal_nelec","If you use smearing, the number of bands should be larger than occupied bands.");
		}
	}
	
	
	if(NBANDS == 0)
	{
		// mohan add 2011-03-31
		NBANDS = int( Mathzone::Max3( occupied_bands , mag.get_nelup() , mag.get_neldw() )  + 1.0e-6);
		if( Occupy::gauss() || Occupy::tetra() )
		{
			//========================================================
			// max function can only be used for bool,double,int,float
			// otherwise it's a bug
			//=========================================================
                         NBANDS = static_cast<int>(1.2*occupied_bands); // pengfei 2014-10-13
                         if (NSPIN == 2)
                         {
                             NBANDS = static_cast<int>(1.2*occupied_bands); //zhengdy 2020-05-12
                             //if (NBANDS > NLOCAL)
                             //    NBANDS = NLOCAL;
                         }
                         NBANDS = NBANDS>6 ? NBANDS : 6 ; 

	//		NBANDS = static_cast<int>( Mathzone::Max3( 1.2*occupied_bands, 
	//					1.2*mag.get_nelup() , 
	//					1.2*mag.get_neldw() ) );
		}
		if(NSPIN == 4) NBANDS = NBANDS * 2;//added by zhengdy-soc
		AUTO_SET("NBANDS",NBANDS);
	}
	else if ( CALCULATION=="scf" || CALCULATION=="md" || CALCULATION=="relax") //pengfei 2014-10-13
	{
		if(NBANDS < occupied_bands) WARNING_QUIT("unitcell","Too few bands!");
		if(NBANDS < mag.get_nelup() ) 
		{
			OUT(ofs_running,"nelup",mag.get_nelup());
			WARNING_QUIT("unitcell","Too few spin up bands!");
		}
		if(NBANDS < mag.get_neldw() ) 
		{
			WARNING_QUIT("unitcell","Too few spin down bands!");
		}
	}

	OUT(ofs_running,"NBANDS",NBANDS);
	return;
}

//===========================================
// calculate the total number of local basis
// Target : nwfc, lmax,
// 			atoms[].stapos_wf
// 			NBANDS
//===========================================
void UnitCell_pseudo::cal_nwfc()
{
	if(test_pseudo_cell) TITLE("UnitCell_pseudo","cal_nwfc");
	assert(ntype>0);
	assert(nat>0);

	//===========================
	// (1) set iw2l, iw2n, iw2m
	//===========================
	for(int it=0; it<ntype; it++)
	{
		this->atoms[it].set_index();
	}

	//===========================
	// (2) set namax and nwmax
	//===========================
	this->namax = 0;
	this->nwmax = 0;
	for(int it=0; it<ntype; it++)
	{
		this->namax = std::max( atoms[it].na, namax );
		this->nwmax = std::max( atoms[it].nw, nwmax );
	}
	assert(namax>0);
//	OUT(ofs_running,"max input atom number",namax);
//	OUT(ofs_running,"max wave function number",nwmax);	

	//===========================
	// (3) set nwfc and stapos_wf
	//===========================
	NLOCAL = 0;
	for(int it=0; it<ntype; it++)
	{
		atoms[it].stapos_wf = NLOCAL;
		const int nlocal_it = atoms[it].nw * atoms[it].na;
                if(!NONCOLIN) NLOCAL += nlocal_it;
		else NLOCAL += nlocal_it * 2;//zhengdy-soc
//		stringstream ss1;
//		ss1 << "number of local orbitals for species " << it;
//		if(LOCAL_BASIS)OUT(ofs_running,ss1.str(),nlocal_it);
//		OUT(ofs_running,"start position of local orbitals",atoms[it].stapos_wf);
	}
	
	//OUT(ofs_running,"NLOCAL",NLOCAL);
            ofs_running << " " << setw(40) << "NLOCAL" << " = " << NLOCAL <<endl;
	//========================================================
	// (4) set index for iat2it, iat2ia, itia2iat, itiaiw2iwt
	//========================================================

	this->set_iat2it();
	
	delete[] iat2ia;
	this->iat2ia = new int[nat];// bug fix 2009-3-8

	// mohan add 2010-09-26
	assert(NLOCAL>0);
	delete[] iwt2iat;
	this->iwt2iat = new int[NLOCAL];

	this->itia2iat.create(ntype, namax);
	this->itiaiw2iwt.create(ntype, namax, nwmax*NPOL);
	int iat=0;
	int iwt=0;
	for(int it = 0;it < ntype;it++)
	{
		for(int ia=0; ia<atoms[it].na; ia++)
		{
			this->itia2iat(it, ia) = iat;
			this->iat2ia[iat] = ia;
			for(int iw=0; iw<atoms[it].nw * NPOL; iw++)
			{
				this->itiaiw2iwt(it, ia, iw) = iwt;
				this->iwt2iat[iwt] = iat;
				++iwt;
			}
			++iat;
		}	
	}

	//========================
	// (5) set lmax and nmax
	//========================
	this->lmax = 0;	
	this->nmax = 0;
	for(int it=0; it<ntype; it++)
	{
		lmax = std::max( lmax, atoms[it].nwl);
		for(int l=0; l<atoms[it].nwl+1; l++)
		{
			nmax = std::max( nmax, atoms[it].l_nchi[l] );
		}

		int nchi = 0;
		for(int l=0; l<atoms[it].nwl+1; l++)
		{
			nchi += atoms[it].l_nchi[l];
		}
		this->nmax_total = std::max(nmax_total, nchi);
	}

	//=======================
	// (6) set lmax_ppwf
	//=======================
	this->lmax_ppwf = 0;
	for(int it=0; it< ntype; it++)
	{
		for(int ic=0; ic<ucell.atoms[it].nchi; ic++)
		{
			if( lmax_ppwf < ucell.atoms[it].lchi[ic] )
			{
				this->lmax_ppwf = ucell.atoms[it].lchi[ic]; 
			}
		}
	}

	/*
	for(int it=0; it< ntype; it++)
	{
		cout << " label=" << it << " nbeta=" << ucell.atoms[it].nbeta << endl;
		for(int ic=0; ic<ucell.atoms[it].nbeta; ic++)
		{
			cout << " " << ucell.atoms[it].lll[ic] << endl;
		}
	}
	*/

//	OUT(ofs_running,"lmax between L(pseudopotential)",lmax_ppwf);

	//=====================
	// Use localized basis
	//=====================
	//if(LOCAL_BASIS) xiaohui modify 2013-09-02
	if(BASIS_TYPE=="lcao" || BASIS_TYPE=="lcao_in_pw") //xiaohui add 2013-09-02
	{
		//if(winput::b_recon)
		//{
		//	NBANDS = NLOCAL;
		//}
		AUTO_SET("NBANDS",NBANDS);
	}
	else // plane wave basis
	{
		//if(winput::after_iter && winput::sph_proj)
		//{
		//	if(NBANDS < NLOCAL)
		//	{
		//		WARNING_QUIT("cal_nwfc","NBANDS must > NLOCAL !");
		//	}
		//}
	}

	return;
}

//======================
// Target : meshx
// Demand : atoms[].msh
//======================
void UnitCell_pseudo::cal_meshx()
{
	if(test_pseudo_cell) TITLE("UnitCell_pseudo","cal_meshx");
	this->meshx = 0;
	for (int it = 0;it < this->ntype;it++)
	{
		const int mesh = this->atoms[it].msh;
		if (mesh > this->meshx)
		{
			this->meshx = mesh;
		}
	}
	return;
}

//=========================
// Target : natomwfc
// Demand : atoms[].nchi
// 			atoms[].lchi
// 			atoms[].oc
// 			atoms[].na
//=========================
void UnitCell_pseudo::cal_natomwfc(void)
{
	if(test_pseudo_cell) TITLE("UnitCell_pseudo","cal_natomwfc");
	this->natomwfc = 0;
	for (int it = 0;it < ntype;it++)
	{
		//============================
		// Use pseudo-atomic orbitals
		//============================
		int tmp=0;
		for (int l = 0;l < atoms[it].nchi;l++)
		{
			if (atoms[it].oc[l] >= 0)
			{
				if(NONCOLIN)
				{
					if(atoms[it].has_so)
					{
						tmp += 2 * atoms[it].lchi[l];
						if(fabs(atoms[it].jchi[l] - atoms[it].lchi[l] - 0.5)<1e-6)
						tmp += 2 ;
					}
					else
					{
						tmp += 2 * (2 * atoms[it].lchi[l] + 1);
					}
				}
				else
					tmp += 2 * atoms[it].lchi[l] + 1;
			}
		}
		natomwfc += tmp * atoms[it].na;
	}
	OUT(ofs_running,"initial pseudo atomic orbital number",natomwfc);
	return;
}

void UnitCell_pseudo::print_unitcell_pseudo(const string &fn)
{
	if(test_pseudo_cell) TITLE("UnitCell_pseudo","print_unitcell_pseudo");
	ofstream ofs( fn.c_str() );

	this->print_cell(ofs);
	for (int i = 0;i < ntype;i++)
	{
		atoms[i].print_Atom(ofs);
	}

	ofs.close();
	return;
}


int UnitCell_pseudo::find_type(const string &label)
{
	if(test_pseudo_cell) TITLE("UnitCell_pseudo","find_type");
	assert(ntype>0);
	for(int it=0;it<ntype;it++)
	{
		if(atoms[it].label == label)
		{
			return it;
		}
	}
	WARNING_QUIT("UnitCell_pseudo::find_type","Can not find the atom type!");
	return -1;
}

#ifdef __MPI
void UnitCell_pseudo::bcast_unitcell_pseudo(void)
{
	Parallel_Common::bcast_int( meshx );
	Parallel_Common::bcast_int( natomwfc );
	Parallel_Common::bcast_int( lmax );
	Parallel_Common::bcast_int( lmax_ppwf );
	Parallel_Common::bcast_double( nelec );

	bcast_unitcell();
}

void UnitCell_pseudo::bcast_unitcell_pseudo2(void)
{
	bcast_unitcell2();
}
#endif

void UnitCell_pseudo::check_dtau(void)
{
	for(int it=0; it<ntype; it++)
	{
		Atom* atom1 = &atoms[it];
		for(int ia=0; ia<atoms[it].na; ia++)
		{
			double dx2 = (atom1->taud[ia].x+10000) - int(atom1->taud[ia].x+10000);
			double dy2 = (atom1->taud[ia].y+10000) - int(atom1->taud[ia].y+10000);
			double dz2 = (atom1->taud[ia].z+10000) - int(atom1->taud[ia].z+10000);

			// mohan add 2011-04-07			
			while(dx2 >= 1) 
			{
				ofs_warning << " dx2 is >=1 " << endl;
				dx2 -= 1.0;
			}
			while(dy2 >= 1) 
			{
				ofs_warning << " dy2 is >=1 " << endl;
				dy2 -= 1.0;
			}
			while(dz2 >= 1) 
			{
				ofs_warning << " dz2 is >=1 " << endl;
				dz2 -= 1.0;
			}
			// mohan add 2011-04-07			
			while(dx2<0) 
			{
				ofs_warning << " dx2 is <0 " << endl;
				dx2 += 1.0;
			}
			while(dy2<0) 
			{
				ofs_warning << " dy2 is <0 " << endl;
				dy2 += 1.0;
			}
			while(dz2<0) 
			{
				ofs_warning << " dz2 is <0 " << endl;
				dz2 += 1.0;
			}

			atom1->taud[ia].x = dx2;
			atom1->taud[ia].y = dy2;
			atom1->taud[ia].z = dz2;

			double cx2, cy2, cz2;

			Mathzone::Direct_to_Cartesian(
			atom1->taud[ia].x, atom1->taud[ia].y, atom1->taud[ia].z,
			latvec.e11, latvec.e12, latvec.e13,
			latvec.e21, latvec.e22, latvec.e23,
			latvec.e31, latvec.e32, latvec.e33,
			cx2, cy2, cz2);

			atom1->tau[ia].x = cx2;
			atom1->tau[ia].y = cy2;
			atom1->tau[ia].z = cz2;

	//		cout << setw(15) << dx2 << setw(15) << dy2 << setw(15) << dz2 
	//		<< setw(15) << cx2 << setw(15) << cy2 << setw(15) << cz2
	//		<< endl;
			
		}
	}
	return;
}

//LiuXh add a new function here,
//20180515
void UnitCell_pseudo::setup_cell_after_vc(
        const string &s_pseudopot_dir,
        const string &fn, ofstream &log)
{
    if(MY_RANK == 0)
    {
        //ifstream ifa(fn.c_str(), ios::in);
        //this->read_atom_species_after_vc(ifa);
    }

    assert(lat0 > 0.0);
    this->omega = abs(latvec.Det()) * this->lat0 * lat0 * lat0;
    if(this->omega <= 0)
    {
        WARNING_QUIT("setup_cell_after_vc", "omega <= 0 .");
    }
    else
    {
        ofs_running << endl;
        OUT(ofs_running, "Volume (Bohr^3)", this->omega);
        OUT(ofs_running, "Volume (A^3))", this->omega * pow(BOHR_TO_A, 3));
    }

    //==========================================================
    // Calculate recip. lattice vectors and dot products
    // latvec has the unit of lat0, but G has the unit 2Pi/lat0
    //==========================================================
    this->GT = latvec.Inverse();
    this->G  = GT.Transpose();
    this->GGT = G * GT;
    this->invGGT = GGT.Inverse();

    for(int it=0; it<ucell.ntype; it++)
    {
        Atom* atom = &ucell.atoms[it];
        for(int ia =0;ia< atom->na;ia++)
        {
            atom->tau[ia] = atom->taud[ia] * ucell.latvec;
/*
#ifdef __MPI
Parallel_Common::bcast_double( atom->tau[ia].x );
Parallel_Common::bcast_double( atom->tau[ia].y );
Parallel_Common::bcast_double( atom->tau[ia].z );
Parallel_Common::bcast_double( atom->taud[ia].x );
Parallel_Common::bcast_double( atom->taud[ia].y );
Parallel_Common::bcast_double( atom->taud[ia].z );
#endif
*/
        }
    }
#ifdef __MPI
    MPI_Barrier(MPI_COMM_WORLD);
    for (int i=0;i<ucell.ntype;i++)
    {
        ucell.atoms[i].bcast_atom(); // bcast tau array
    }
#endif

    ofs_running << endl;
    out.printM3(ofs_running,"Lattice vectors: (Cartesian coordinate: in unit of a_0)",latvec);
    out.printM3(ofs_running,"Reciprocal vectors: (Cartesian coordinate: in unit of 2 pi/a_0)",G);

    return;
}
