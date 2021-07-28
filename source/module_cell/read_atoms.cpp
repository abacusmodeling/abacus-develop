#include "unitcell_pseudo.h"
#ifdef __LCAO
#include "../module_orbital/ORB_read.h" // to use 'ORB' -- mohan 2021-01-30
#endif

#ifndef __CELL
#include "../src_pw/global.h"
#endif
#include <cstring>		// Peize Lin fix bug about strcmp 2016-08-02

void UnitCell_pseudo::read_atom_species(ifstream &ifa, ofstream &ofs_running)
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
			if(GlobalV::test_pseudo_cell==2) 
			{
				ofs_running << "\n" << setw(6) << atom_label[i] 
						<< setw(12) << atom_mass[i] 
						<< setw(18) << pseudo_fn[i];
			}

			// Peize Lin test for bsse 2021.04.07
			const string bsse_label = "empty";
			this->atoms[i].flag_empty_element = 
				(search( atom_label[i].begin(), atom_label[i].end(), bsse_label.begin(), bsse_label.end() ) != atom_label[i].end())
				? true : false;
		}
	}
#ifdef __LCAO
	if(GlobalV::BASIS_TYPE=="lcao" || GlobalV::BASIS_TYPE=="lcao_in_pw")
	{
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

//				GlobalV::ofs_running << " For atom type " << i + 1 << endl;
//			    GlobalV::ofs_running << " Read in numerical orbitals from file " << ofile << endl;
//			    GlobalV::ofs_running << " Read in nonlocal projectors from file " << nfile << endl;
				
			}
		}	
		// caoyu add 2021-03-16
		if (SCAN_BEGIN(ifa, "NUMERICAL_DESCRIPTOR")) {
			ifa >> ORB.descriptor_file;
		}
	}

	// Peize Lin add 2016-09-23
	if( Exx_Global::Hybrid_Type::HF   == exx_lcao.info.hybrid_type || 
	    Exx_Global::Hybrid_Type::PBE0 == exx_lcao.info.hybrid_type || 
		Exx_Global::Hybrid_Type::HSE  == exx_lcao.info.hybrid_type )
	{
		if( SCAN_BEGIN(ifa, "ABFS_ORBITAL") )
		{
			for(int i=0; i<ntype; i++)
			{
				string ofile;
				ifa >> ofile;
				exx_lcao.info.files_abfs.push_back(ofile);
			}
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
	if(latName=="test"){	
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
		}
		if( SCAN_BEGIN(ifa, "LATTICE_PARAMETERS") )
		{
			WARNING_QUIT("UnitCell_pseudo::read_atom_species","do not use LATTICE_PARAMETERS without explicit specification of lattice type");
		}
	}//supply lattice vectors
	else{
		if( SCAN_BEGIN(ifa, "LATTICE_VECTORS") )
		{
			WARNING_QUIT("UnitCell_pseudo::read_atom_species","do not use LATTICE_VECTORS along with explicit specification of lattice type");
		}
		if(latName=="sc"){//simple-cubic
			latvec.e11 = 1.0; latvec.e12 = 0.0; latvec.e13 = 0.0;
			latvec.e21 = 0.0; latvec.e22 = 1.0;	latvec.e23 = 0.0;
			latvec.e31 = 0.0; latvec.e32 = 0.0;	latvec.e33 = 1.0;
		}
		else if(latName=="fcc"){//face-centered cubic
			latvec.e11 =-0.5; latvec.e12 = 0.0; latvec.e13 = 0.5;
			latvec.e21 = 0.0; latvec.e22 = 0.5;	latvec.e23 = 0.5;
			latvec.e31 =-0.5; latvec.e32 = 0.5;	latvec.e33 = 0.0;
		}
		else if(latName=="bcc"){//body-centered cubic
			latvec.e11 = 0.5; latvec.e12 = 0.5; latvec.e13 = 0.5;
			latvec.e21 =-0.5; latvec.e22 = 0.5;	latvec.e23 = 0.5;
			latvec.e31 =-0.5; latvec.e32 =-0.5;	latvec.e33 = 0.5;
		}
		else if(latName=="hexagonal"){//hexagonal
			double e22 = sqrt(3.0) / 2.0;
			latvec.e11 = 1.0; latvec.e12 = 0.0; latvec.e13 = 0.0;
			latvec.e21 =-0.5; latvec.e22 = e22; latvec.e23 = 0.0;
			latvec.e31 = 0.0; latvec.e32 = 0.0;	latvec.e33 = 0.0;
			if( SCAN_BEGIN(ifa, "LATTICE_PARAMETERS") )
			{
				READ_VALUE(ifa, latvec.e33);
			}
		}
		else if(latName=="trigonal"){//trigonal
			double t1 = 0.0;
			double t2 = 0.0;
			if( SCAN_BEGIN(ifa, "LATTICE_PARAMETERS") )
			{
				double cosab=0.0;
				READ_VALUE(ifa, cosab);
				t1 = sqrt(1.0 + 2.0*cosab);
				t2 = sqrt(1.0 - cosab);
			}
			double e11 = t2 / sqrt(2.0);
			double e12 = -t2 / sqrt(6.0);
			double e13 = t1 / sqrt(3.0);
			double e22 = sqrt(2.0) * t2 / sqrt(3.0);
		
			latvec.e11 = e11; latvec.e12 = e12; latvec.e13 = e13;
			latvec.e21 = 0.0; latvec.e22 = e22;	latvec.e23 = e13;
			latvec.e31 =-e11; latvec.e32 = e12;	latvec.e33 = e13;
		}
		else if(latName=="st"){//simple tetragonal
			latvec.e11 = 1.0; latvec.e12 = 0.0; latvec.e13 = 0.0;
			latvec.e21 = 0.0; latvec.e22 = 1.0; latvec.e23 = 0.0;
			latvec.e31 = 0.0; latvec.e32 = 0.0;	latvec.e33 = 0.0;
			if( SCAN_BEGIN(ifa, "LATTICE_PARAMETERS") )
			{
				READ_VALUE(ifa, latvec.e33);
			}
		}
		else if(latName=="bct"){//body-centered tetragonal
			double cba = 0.0;
			if( SCAN_BEGIN(ifa, "LATTICE_PARAMETERS") )
			{
				READ_VALUE(ifa, cba);
				cba = cba / 2.0;
			}
			latvec.e11 = 0.5; latvec.e12 =-0.5; latvec.e13 = cba;
			latvec.e21 = 0.5; latvec.e22 = 0.5; latvec.e23 = cba;
			latvec.e31 =-0.5; latvec.e32 =-0.5;	latvec.e33 = cba;
		}
		else if(latName=="so"){//simple orthorhombic
			latvec.e11 = 1.0; latvec.e12 = 0.0; latvec.e13 = 0.0;
			latvec.e21 = 0.0; latvec.e22 = 0.0;	latvec.e23 = 0.0;
			latvec.e31 = 0.0; latvec.e32 = 0.0;	latvec.e33 = 0.0;
			if( SCAN_BEGIN(ifa, "LATTICE_PARAMETERS") )
			{
				ifa >> latvec.e22;
				READ_VALUE(ifa, latvec.e33);
			}
		}
		else if(latName=="baco"){//base-centered orthorhombic
			latvec.e11 = 0.5; latvec.e12 = 0.0; latvec.e13 = 0.0;
			latvec.e21 =-0.5; latvec.e22 = 0.0;	latvec.e23 = 0.0;
			latvec.e31 = 0.0; latvec.e32 = 0.0;	latvec.e33 = 0.0;
			if( SCAN_BEGIN(ifa, "LATTICE_PARAMETERS") )
			{
				ifa >> latvec.e12;
				latvec.e12 = latvec.e12 / 2.0;
				latvec.e22 = latvec.e12;
				READ_VALUE(ifa, latvec.e33);
			}
		}
		else if(latName=="fco"){//face-centered orthorhombic
			double bba = 0.0; double cba = 0.0;
			if( SCAN_BEGIN(ifa, "LATTICE_PARAMETERS") )
			{
				ifa >> bba;
				READ_VALUE(ifa, cba);
				bba = bba / 2.0; cba = cba / 2.0;
			}
			latvec.e11 = 0.5; latvec.e12 = 0.0; latvec.e13 = cba;
			latvec.e21 = 0.5; latvec.e22 = bba;	latvec.e23 = 0.0;
			latvec.e31 = 0.0; latvec.e32 = bba;	latvec.e33 = cba;
		}
		else if(latName=="bco"){//body-centered orthorhombic
			double bba = 0.0; double cba = 0.0;
			if( SCAN_BEGIN(ifa, "LATTICE_PARAMETERS") )
			{
				ifa >> bba;
				READ_VALUE(ifa, cba);
				bba = bba / 2.0; cba = cba / 2.0;
			}
			latvec.e11 = 0.5; latvec.e12 = bba; latvec.e13 = cba;
			latvec.e21 =-0.5; latvec.e22 = bba;	latvec.e23 = cba;
			latvec.e31 =-0.5; latvec.e32 =-bba;	latvec.e33 = cba;
		}
		else if(latName=="sm"){//simple monoclinic
			double bba = 0.0; double cba = 0.0;
			double cosab = 0.0;
			double e21 = 0.0; double e22 = 0.0;
			if( SCAN_BEGIN(ifa, "LATTICE_PARAMETERS") )
			{
				ifa >> bba >> cba;
				READ_VALUE(ifa, cosab);
				e21 = bba * cosab;
				e22 = bba * sqrt(1.0-cosab*cosab);
			}
			latvec.e11 = 1.0; latvec.e12 = 0.0; latvec.e13 = 0.0;
			latvec.e21 = e21; latvec.e22 = e22;	latvec.e23 = 0.0;
			latvec.e31 = 0.0; latvec.e32 = 0.0;	latvec.e33 = cba;
		}
		else if(latName=="bacm"){//base-centered monoclinic
			double bba = 0.0; double cba = 0.0;
			double cosab = 0.0;
			double e21 = 0.0; double e22 = 0.0;
			if( SCAN_BEGIN(ifa, "LATTICE_PARAMETERS") )
			{
				ifa >> bba >> cba;
				READ_VALUE(ifa, cosab);
				e21 = bba * cosab;
				e22 = bba * sqrt(1.0-cosab*cosab);
				cba = cba / 2.0;
			}
			latvec.e11 = 0.5; latvec.e12 = 0.0; latvec.e13 =-cba;
			latvec.e21 = e21; latvec.e22 = e22;	latvec.e23 = 0.0;
			latvec.e31 = 0.5; latvec.e32 = 0.0;	latvec.e33 = cba;
		}
		else if(latName=="triclinic"){//triclinic
			double bba = 0.0; double cba = 0.0;
			double cosab = 0.0; double cosac = 0.0;
			double cosbc = 0.0; double sinab = 0.0;
			double term = 0.0;
			if( SCAN_BEGIN(ifa, "LATTICE_PARAMETERS") )
			{
				ifa >> bba >> cba >> cosab >> cosac;
				READ_VALUE(ifa, cosbc);
				sinab = sqrt(1.0-cosab*cosab);
			}
			latvec.e11 = 1.0; latvec.e12 = 0.0; latvec.e13 = 0.0;
			latvec.e21 = bba * cosab;
			latvec.e22 = bba * sinab;
			latvec.e23 = 0.0;
			latvec.e31 = cba * cosac;
			latvec.e32 = cba * (cosbc - cosac*cosab/sinab);
			term = 1.0 + 2.0 * cosab*cosac*cosbc - cosab*cosab - cosac*cosac - cosbc*cosbc;
			term = sqrt(term)/sinab;
			latvec.e33 = cba * term;
		}
		else{ 
			cout << "latname is : " << latName << endl;
			WARNING_QUIT("UnitCell_pseudo::read_atom_species","latname not supported!");
		}
	}

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
	return;
}

// Read atomic positions
// return 1: no problem.
// return 0: some problems.
bool UnitCell_pseudo::read_atom_positions(ifstream &ifpos, ofstream &ofs_running, ofstream &ofs_warning)
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
			return 0; // means something wrong
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
				{
					found = true;
				}
			}
			if(!found)
			{
				ofs_warning << " Label read from ATOMIC_POSITIONS is " << this->atoms[it].label << endl; 
				ofs_warning << " Lable from ATOMIC_SPECIES is " << this->atom_label[it] << endl;
				return 0;
			}
			READ_VALUE(ifpos, magnet.start_magnetization[it] );

			OUT(ofs_running, "atom label",atoms[it].label);
			#ifndef __SYMMETRY
			if(GlobalV::NSPIN==4)//added by zhengdy-soc
			{
				if(GlobalV::NONCOLIN)
				{
					magnet.m_loc_[it].x = magnet.start_magnetization[it] *
							sin(magnet.angle1_[it]) * cos(magnet.angle2_[it]);
					magnet.m_loc_[it].y = magnet.start_magnetization[it] *
							sin(magnet.angle1_[it]) * sin(magnet.angle2_[it]);
					magnet.m_loc_[it].z = magnet.start_magnetization[it] *
							cos(magnet.angle1_[it]);
				}
				else
				{
					magnet.m_loc_[it].x = 0;
					magnet.m_loc_[it].y = 0;
					magnet.m_loc_[it].z = magnet.start_magnetization[it];
				}

				OUT(ofs_running, "noncollinear magnetization_x",magnet.m_loc_[it].x);
				OUT(ofs_running, "noncollinear magnetization_y",magnet.m_loc_[it].y);
				OUT(ofs_running, "noncollinear magnetization_z",magnet.m_loc_[it].z);

				ZEROS(magnet.ux_ ,3);
			}
			else if(GlobalV::NSPIN==2)
			{
				magnet.m_loc_[it].x = magnet.start_magnetization[it];
				OUT(ofs_running, "start magnetization",magnet.start_magnetization[it]);
			}
			else if(GlobalV::NSPIN==1)
			{
				OUT(ofs_running, "start magnetization","FALSE");
			}
			#endif

			//===========================================
			// (2) read in numerical orbital information
			// int atoms[it].nwl
			// int* atoms[it].l_nchi;
			//===========================================
#ifndef __SYMMETRY
#ifdef __LCAO
			if (GlobalV::BASIS_TYPE == "lcao" || GlobalV::BASIS_TYPE == "lcao_in_pw")
			{    
				ifstream ifs(ORB.orbital_file[it].c_str(), ios::in);  // pengfei 2014-10-13

				// mohan add return 2021-04-26
				if (!ifs)
				{
					cout << " Element index " << it+1 << endl;
					cout << " orbital file: " << ORB.orbital_file[it] << endl;
					WARNING("read_atom_positions","ABACUS Cannot find the ORBITAL file (basis sets)");
					return 0; // means something wrong
				}

				char word[80];
				this->atoms[it].nw = 0;
				int L =0;

				while (ifs.good())
				{
					ifs >> word;
					if (strcmp("Lmax", word) == 0)
					{
						READ_VALUE(ifs, this->atoms[it].nwl);
						delete[] this->atoms[it].l_nchi;
						this->atoms[it].l_nchi = new int[ this->atoms[it].nwl+1];
					}
					assert(this->atoms[it].nwl<10);

					if (strcmp("Cutoff(a.u.)", word) == 0)         // pengfei Li 16-2-29
					{
						READ_VALUE(ifs, this->atoms[it].Rcut);
					}

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
				}
				ifs.close();
			}
			else
#else
			if(GlobalV::BASIS_TYPE == "pw")
#endif
			{
				this->atoms[it].nw = 0;

				this->atoms[it].nwl = 2;
				//cout << lmaxmax << endl;
				if ( lmaxmax != 2 )
				{
					this->atoms[it].nwl = lmaxmax;
				}
				delete[] this->atoms[it].l_nchi;
				this->atoms[it].l_nchi = new int[ this->atoms[it].nwl+1];
				for(int L=0; L<atoms[it].nwl+1; L++)
				{
					this->atoms[it].l_nchi[L] = 1;
					// calculate the number of local basis(3D)
					this->atoms[it].nw += (2*L + 1) * this->atoms[it].l_nchi[L];
					stringstream ss;
					ss << "L=" << L << ", number of zeta";
					OUT(ofs_running,ss.str(),atoms[it].l_nchi[L]);
				}
			} // end basis type
#endif


			//OUT(GlobalV::ofs_running,"Total number of local orbitals",atoms[it].nw);

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
							//atoms[it].mag[ia] = std::atof(mags.c_str());
							atoms[it].mag[ia] = 0.0;
						}
					}
					else
					{
						ifpos >> v.x >> v.y >> v.z
							>> mv.x >> mv.y >> mv.z;
						string mags;
						std::getline( ifpos, mags );
						// change string to double.
						//atoms[it].mag[ia] = std::atof(mags.c_str());
						atoms[it].mag[ia] = 0.0;
					}	

					if(Coordinate=="Direct")
					{
						// change v from direct to cartesian,
						// the unit is GlobalC::pw.lat0
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

	//GlobalV::ofs_running << "\n Output nearest atom not considering periodic boundary condition" << endl;
	//GlobalV::ofs_running << " " << setw(5) << "TYPE" << setw(6) << "INDEX" 
	//<< setw(20) << "NEAREST(Bohr)" 
	//<< setw(20) << "NEAREST(Angstrom)" << endl; 
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
						norm = diff.norm() * lat0;
						if( shortest_norm > norm )
						{
							shortest_norm = norm;
							//nearest_atom_type = T2;
							//nearest_atom_index = I2;
						}
						if( norm < tolerence_bohr ) // unit is Bohr
						{	
							GlobalV::ofs_warning << " two atoms are too close!" << endl;
							GlobalV::ofs_warning << " type:" << this->atoms[T1].label << " atom " << I1 + 1 << endl; 
							GlobalV::ofs_warning << " type:" << this->atoms[T2].label << " atom " << I2 + 1 << endl; 
							GlobalV::ofs_warning << " distance = " << norm << " Bohr" << endl;
							return 0;
						}
					}
				}
			}
			//GlobalV::ofs_running << " " << setw(5) << atoms[T1].label << setw(6) << I1+1 
			//<< setw(20) << shortest_norm  
			//<< setw(20) << shortest_norm * BOHR_TO_A << endl;
		}
	}

	timer::tick("UnitCell_pseudo","check_tau");
	return 1;
}

void UnitCell_pseudo::print_stru_file(const string &fn, const int &type)const
{
	TITLE("UnitCell_pseudo","print_stru_file");
	
	if(GlobalV::MY_RANK!=0) return;

	ofstream ofs(fn.c_str());

	ofs << "ATOMIC_SPECIES" << endl;
	ofs << setprecision(12);

	for(int it=0; it<ntype; it++)
	{
                //modified by zhengdy 2015-07-24
		ofs << atom_label[it] << " " << atom_mass[it] << " " << pseudo_fn[it] << endl;
	}

#ifdef __LCAO
	if(GlobalV::BASIS_TYPE=="lcao" || GlobalV::BASIS_TYPE=="lcao_in_pw") //xiaohui add 2013-09-02. Attention...
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
			ofs << magnet.start_magnetization[it] << " #magnetism" << endl;
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
			ofs << magnet.start_magnetization[it] << " #magnetism" << endl;
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
        GlobalV::ofs_running << "\n CARTESIAN COORDINATES ( UNIT = " << lat0 << " Bohr )." << endl;
        GlobalV::ofs_running << setw(13) << " atom"
        //<< setw(20) << "x" 
        //<< setw(20) << "y" 
        //<< setw(20) << "z" 
        //<< " mag"
        << setw(20) << "x"
        << setw(20) << "y"
        << setw(20) << "z"
        << setw(20) << "mag"
        << endl;
        GlobalV::ofs_running << setprecision(12);

        int iat=0;
        for(int it=0; it<ntype; it++)
        {
            for (int ia = 0; ia < atoms[it].na; ia++)
            {
                stringstream ss;
                ss << "tauc_" << atoms[it].label << ia+1;

                GlobalV::ofs_running << " " << setw(12) << ss.str()
                //<< setw(20) << atoms[it].tau[ia].x 
                //<< setw(20) << atoms[it].tau[ia].y
                //<< setw(20) << atoms[it].tau[ia].z
                //<< " " << atoms[it].mag[ia]
                << setw(20) << atoms[it].tau[ia].x
                << setw(20) << atoms[it].tau[ia].y
                << setw(20) << atoms[it].tau[ia].z
                //<< setw(20) << atoms[it].mag[ia]
                << setw(20) << magnet.start_magnetization[it]
                << endl;

                ++iat;
            }
        }
    }

    if(Coordinate == "Direct")
    {
        GlobalV::ofs_running << "\n DIRECT COORDINATES" << endl;
        GlobalV::ofs_running << setw(13) << " atom"
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
        for(int it=0; it<ntype; it++)
        {
            for (int ia = 0; ia < atoms[it].na; ia++)
            {
                stringstream ss;
                ss << "taud_" << atoms[it].label << ia+1;

                GlobalV::ofs_running << " " << setw(12) << ss.str()
                //<< setw(20) << atoms[it].taud[ia].x
                //<< setw(20) << atoms[it].taud[ia].y
                //<< setw(20) << atoms[it].taud[ia].z
                //<< " " << atoms[it].mag[ia]
                << setw(20) << atoms[it].taud[ia].x
                << setw(20) << atoms[it].taud[ia].y
                << setw(20) << atoms[it].taud[ia].z
                //<< setw(20) << atoms[it].mag[ia]
                << setw(20) << magnet.start_magnetization[it]
                << endl;

                ++iat;
            }
        }
    }

	GlobalV::ofs_running << endl;
	return;
}	


int UnitCell_pseudo::find_type(const string &label)
{
	if(GlobalV::test_pseudo_cell) TITLE("UnitCell_pseudo","find_type");
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
				GlobalV::ofs_warning << " dx2 is >=1 " << endl;
				dx2 -= 1.0;
			}
			while(dy2 >= 1) 
			{
				GlobalV::ofs_warning << " dy2 is >=1 " << endl;
				dy2 -= 1.0;
			}
			while(dz2 >= 1) 
			{
				GlobalV::ofs_warning << " dz2 is >=1 " << endl;
				dz2 -= 1.0;
			}
			// mohan add 2011-04-07			
			while(dx2<0) 
			{
				GlobalV::ofs_warning << " dx2 is <0 " << endl;
				dx2 += 1.0;
			}
			while(dy2<0) 
			{
				GlobalV::ofs_warning << " dy2 is <0 " << endl;
				dy2 += 1.0;
			}
			while(dz2<0) 
			{
				GlobalV::ofs_warning << " dz2 is <0 " << endl;
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
