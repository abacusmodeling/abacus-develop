#include "lcao_matrix.h"
#include "global_fp.h"

LCAO_Matrix::LCAO_Matrix()
{
	// for gamma_only
	Sloc = new double[1];
	Hloc_fixed = new double[1];
	Hloc = new double[1];

	// for many k points
	Sloc2 = new complex<double>[1];
	Hloc_fixed2 = new complex<double>[1];
	Hloc2 = new complex<double>[1];
}

LCAO_Matrix::~LCAO_Matrix()
{
 	// delete matrix for gamma_only.
    delete[] Sloc;
    delete[] Hloc_fixed;
    delete[] Hloc;

    // delete matrix for many k points
    delete[] Sloc2;
    delete[] Hloc_fixed2;
    delete[] Hloc2;	
}

void LCAO_Matrix::divide_HS_in_frag(void)
{
	TITLE("LCAO_Matrix","divide_HS_in_frag");

	ofs_running << "\n SETUP THE DIVISION OF H/S MATRIX" << endl;
	
	// (1) calculate nrow, ncol, nloc.
	//if (DIAGO_TYPE=="hpseps" || DIAGO_TYPE=="scalpack" || DIAGO_TYPE=="selinv") xiaohui modify 2013-09-02
	if (KS_SOLVER=="hpseps" || KS_SOLVER=="scalpack" || KS_SOLVER=="selinv") //xiaohui add 2013-09-02
	{
		ofs_running << " divide the H&S matrix using 2D block algorithms." << endl;
#ifdef __MPI
		ParaO.divide_HS_2d(DIAG_WORLD);
#else
		WARNING_QUIT("LCAO_Matrix::init","diago_type = 'HPSEPS' is not available for series version.\n You can use 'LAPACK' instead.");
#endif
	}
    	//else if (DIAGO_TYPE == "canonical"
        //    	|| DIAGO_TYPE == "trace_resetting"
        //    	|| DIAGO_TYPE == "trace_correcting") xiaohui modify 2013-09-02
	else if(KS_SOLVER == "canonical"
		|| KS_SOLVER == "trace_resetting"
		|| KS_SOLVER == "trace_correcting") //xiaohui add 2013-09-02
	{
		ofs_running << " divide the H&S matrix according to atoms." << endl;
		// done nothing.
		// the nloc is calculated in "ParaA.set_trace"
	}
	else
	{
		ParaO.nloc = NLOCAL * NLOCAL;
	}


	// (2) set the trace, then we 
	// can calculate the nnr.
	// for 2d: calculate ParaO.nloc first, then trace_loc_row and trace_loc_col
	// for O(N): calculate the three together.
	ParaO.set_trace();



	// (3) allocate matrix.
	if(GAMMA_ONLY_LOCAL)
	{
		if(BFIELD)
		{
			// if b field is chosen,
			// the hamiltonian is complex<double>.
			allocate_HS_k(ParaO.nloc);	
		}
		else
		{
			allocate_HS_gamma(ParaO.nloc);
		}
	}
	else
	{
		allocate_HS_k(ParaO.nloc);
	}


	return;
}

void LCAO_Matrix::allocate_HS_gamma(const int &nloc)
{
	TITLE("LCAO_Matrix","allocate_HS_gamma");

	OUT(ofs_running,"nloc",nloc);
	if(nloc==0) return; //mohan fix bug 2012-05-25

	delete[] Sloc;
	delete[] Hloc_fixed;
	delete[] Hloc;

	this->Sloc = new double[nloc];
	this->Hloc_fixed = new double[nloc];
	this->Hloc = new double[nloc];

	ZEROS(Sloc,nloc);
	ZEROS(Hloc_fixed,nloc);
	ZEROS(Hloc,nloc);

	return;
}


void LCAO_Matrix::allocate_HS_k(const int &nloc)
{
	TITLE("LCAO_Matrix","allocate_HS_k");

	OUT(ofs_running,"nloc",nloc);
	if(nloc==0) return; //mohan fix bug 2012-05-25

	delete[] Sloc2;
	delete[] Hloc_fixed2;
	delete[] Hloc2;

	this->Sloc2 = new complex<double>[nloc];
	this->Hloc_fixed2 = new complex<double>[nloc];
	this->Hloc2 = new complex<double>[nloc];

	ZEROS(Sloc2,nloc);
	ZEROS(Hloc_fixed2,nloc);
	ZEROS(Hloc2,nloc);
	
	return;
}

void LCAO_Matrix::allocate_HS_R(const int &nnR)
{
	delete[] HlocR;
	delete[] SlocR;
	delete[] Hloc_fixedR;	

	this->HlocR = new double[nnR];
	this->SlocR = new double[nnR];
	this->Hloc_fixedR = new double[nnR];

	ZEROS(HlocR, nnR);
	ZEROS(SlocR, nnR);
	ZEROS(Hloc_fixedR, nnR);

	return;
}

void LCAO_Matrix::set_HSgamma(const int &iw1_all, const int &iw2_all, const double &v, const char &dtype)
{
    // use iw1_all and iw2_all to set Hloc
    // becareful! The ir and ic may < 0!!!!!!!!!!!!!!!!
    const int ir = ParaO.trace_loc_row[ iw1_all ];
    const int ic = ParaO.trace_loc_col[ iw2_all ];
    const int index = ir * ParaO.ncol + ic;
   
   	if( index >= ParaO.nloc)
	{
		cout << " iw1_all = " << iw1_all << endl;
		cout << " iw2_all = " << iw2_all << endl;
		cout << " ir = " << ir << endl;
		cout << " ic = " << ic << endl;
		cout << " index = " << index << endl;
		cout << " ParaO.nloc = " << ParaO.nloc << endl;
		WARNING_QUIT("LCAO_Matrix","set_HSgamma");
	}	 

	// S : S matrix element.
	// T : T matrix element.
	// N : nonlocal H matrix element.
	// L : local H matrix element.

    if (dtype=='S')
    {
        this->Sloc[index] += v;
    }
    else if (dtype=='T' || dtype=='N')
    {
        this->Hloc_fixed[index] += v;
    }
    else if (dtype=='L')
    {
        this->Hloc[index] += v;
    }

    return;
}

void LCAO_Matrix::set_HSk(const int &iw1_all, const int &iw2_all, const complex<double> &v, const char &dtype)
{
    // use iw1_all and iw2_all to set Hloc
    // becareful! The ir and ic may < 0!!!!!!!!!!!!!!!!
    const int ir = ParaO.trace_loc_row[ iw1_all ];
    const int ic = ParaO.trace_loc_col[ iw2_all ];
    const int index = ir * ParaO.ncol + ic;
    assert(index < ParaO.nloc);

    if (dtype=='S')//overlap Hamiltonian.
    {
        this->Sloc2[index] += v;
    }
    else if (dtype=='T' || dtype=='N')// kinetic and nonlocal Hamiltonian.
    {
        this->Hloc_fixed2[index] += v;
    }
    else if (dtype=='L') // Local potential Hamiltonian.
    {
        this->Hloc2[index] += v;
    }
	else
	{
		WARNING_QUIT("LCAO_Matrix","set_HSk");
	}

    return;
}

void LCAO_Matrix::set_force
(
    const int &iw1_all,
    const int &iw2_all,
    const double& vx,
    const double& vy,
    const double& vz,
    const char &dtype)
{
    // use iw1_all and iw2_all to set Hloc
    // becareful! The ir and ic may < 0!!!!!!!!!!!!!!!!
    const int ir = ParaO.trace_loc_row[ iw1_all ];
    const int ic = ParaO.trace_loc_col[ iw2_all ];
    const int index = ir * ParaO.ncol + ic;
	
	if( index >= ParaO.nloc)
	{
		cout << " iw1_all = " << iw1_all << endl;
		cout << " iw2_all = " << iw2_all << endl;
		cout << " ir = " << ir << endl;
		cout << " ic = " << ic << endl;
		cout << " index = " << index << endl;
		cout << " ParaO.nloc = " << ParaO.nloc << endl;
		WARNING_QUIT("LCAO_Matrix","set_force");
	}	 

    if (dtype == 'S')
    {
        this->DSloc_x[index] += vx;
        this->DSloc_y[index] += vy;
        this->DSloc_z[index] += vz;
    }
    else if (dtype == 'T')
    {
		// notice, the sign is '-', minus.
        this->DHloc_fixed_x[index] -= vx;
        this->DHloc_fixed_y[index] -= vy;
        this->DHloc_fixed_z[index] -= vz;
    }
    else if (dtype == 'N')
	{
		this->DHloc_fixed_x[index] += vx;
		this->DHloc_fixed_y[index] += vy;
		this->DHloc_fixed_z[index] += vz;
	}

	return;
}


void LCAO_Matrix::zeros_HSgamma(const char &mtype)
{
	if (mtype=='S') ZEROS(Sloc,ParaO.nloc);
	else if (mtype=='T') ZEROS(Hloc_fixed,ParaO.nloc);
	else if (mtype=='H') ZEROS(Hloc,ParaO.nloc);
	return;
}

void LCAO_Matrix::zeros_HSk(const char &mtype)
{
	if (mtype=='S') ZEROS(Sloc2,ParaO.nloc);
	else if (mtype=='T') ZEROS(Hloc_fixed2,ParaO.nloc);
	else if (mtype=='H') ZEROS(Hloc2,ParaO.nloc);
	return;
}

void LCAO_Matrix::zeros_HSR(const char &mtype, const int &nnr)
{
	if (mtype=='S') ZEROS(SlocR, nnr);
	else if (mtype=='T') ZEROS(Hloc_fixedR, nnr);
	else if (mtype=='H') ZEROS(HlocR, nnr);
	return;
}

void LCAO_Matrix::print_HSk(const char &mtype, const char &vtype, const double &accuracy)
{
	TITLE("LCAO_Matrix","print_HSk");
	if(mtype=='S') cout << "Sloc2 matrix" << endl;
	else if(mtype=='T') cout << "Hloc_fixed2 matrix" << endl;
	else if(mtype=='H') cout << "Hloc matrix" << endl;
	else
	{
		WARNING_QUIT("LCAO_Matrix::print_HSk","Check input parameter: mtype.");
	}

	if(vtype=='C') cout << " Output norm."  << endl;
	else if(vtype=='R') cout << " Output real part."  << endl;
	else if(vtype=='I') cout << " Output imag part."  << endl;


	cout << setprecision(8) << endl;
	for(int i=0; i<ParaO.nrow; i++)
	{
		cout << " " ;
		for(int j=0; j<ParaO.ncol; j++)
		{
			const int index = i * ParaO.ncol + j;
			double v=-888.888;//wrong number
			if(vtype=='R')
			{
				if(mtype=='S') v = Sloc2[index].real();
				else if(mtype=='T') v = Hloc_fixed2[index].real();
				else if(mtype=='H') v = Hloc2[index].real();
			}
			else if(vtype=='C')
			{
				if(mtype=='S') v = sqrt( norm ( Sloc2[index] ) );
				else if(mtype=='T') v = sqrt( norm ( Hloc_fixed2[index] ) );
				else if(mtype=='H') v = sqrt( norm ( Hloc2[index] ) );
			}
			else if(vtype=='I')
			{
				if(mtype=='S') v = Sloc2[index].imag();
				else if(mtype=='T') v = Hloc_fixed2[index].imag();
				else if(mtype=='H') v = Hloc2[index].imag();
			}

			if( abs(v) > accuracy )
			{
				cout << setw(15) << v;
			}
			else
			{
				cout << setw(15) << "0"; 
			}
		}
		cout << endl;
	}
	cout << endl;
	cout << setprecision(6) << endl;
	return;
}


void LCAO_Matrix::print_HSgamma(const char &mtype)
{
	TITLE("Parallel_Orbitals","print_HSgamma");

	ofs_running << " " << mtype << " matrix" << endl;
	ofs_running << " nrow=" << ParaO.nrow << endl;
	ofs_running << " ncol=" << ParaO.ncol << endl;
	ofs_running << " element number = " << ParaO.ncol << endl;
	if (mtype=='S')
	{
		if(!BFIELD)
		{
			cout << setprecision(8);
			cout << " print Sloc" << endl;
			for(int i=0; i<NLOCAL; ++i)
			{
				for(int j=0; j<NLOCAL; ++j)
				{
					double v = Sloc[i*ParaO.ncol+j];
					if( abs(v) > 1.0e-8)
					{
						cout << setw(15) << v;
					}
					else
					{
						cout << setw(15) << "0";
					}
				}//end j
				cout << endl;
			}//end i
		}

		/*
		   ofs_running << " " << setw(10) << "LocalRow" << setw(10) << "LocalCol"
		   << setw(10) << "GlobalRow" << setw(10) << "GloablCol"
		   << setw(10) << "Sloc" << endl;
		   for (int i=0; i<this->nrow; i++)
		   {
		   for (int j=0; j<this->ncol; j++)
		   {
		   if ( abs(Sloc[i * this->ncol + j]) > 1.0e-5 )
		   ofs_running << " " << setw(10) << i << setw(10) << j
		   << setw(10) << MatrixInfo.row_set[i] << setw(10) << MatrixInfo.col_set[j]
		   << setw(10) << Sloc[i * this->ncol + j] << endl;
		   }
		   }
		cout << "\n Smatrix" << endl;
		//ofs_running << setprecision(5) << endl;

		for(int i=0; i<ParaO.nrow; i++)
		{
			for(int j=0; j<ParaO.ncol; j++)
			{
				double s = Sloc[ i * ParaO.ncol + j];
				double h = Hloc[ i * ParaO.ncol + j];
				
				if( abs( h - s ) > 1.0e-3 )
				{
				//	ofs_running << setw(5) << i+1 << setw(5) << j+1 << setw(12) << h << setw(12) << s << endl;
				}
				
				
				if( abs(s) > 1.0e-5 )
				{
					ofs_running << setw(10) << s;
				//	ofs_running << setw(5) << i << setw(5) << j << setw(15) << s << endl;

				}
				else
				{
					ofs_running << setw(10) << "0";
				}
			}
			ofs_running << endl;
		}
		 */
	}
	if (mtype=='T')
	{


		/*
		ofs_running << " " << setw(10) << "LocalRow" << setw(10) << "LocalCol"
			<< setw(10) << "GlobalRow" << setw(10) << "GloablCol"
			<< setw(10) << "Hloc_fixed" << endl;
		for (int i=0; i<ParaO.nrow; i++)
		{
			for (int j=0; j<ParaO.ncol; j++)
			{
//				if ( abs(Hloc_fixed[i * ParaO.ncol + j]) > 1.0e-5 )
				{
					ofs_running << " " << setw(10) << i+1 << setw(10) << j+1;

					//mohan fix bug 2012-02-22
					if(ATOM_DISTRIBUTION!=1)
					{
						ofs_running << setw(10) << ParaO.MatrixInfo.row_set[i] << setw(10) << ParaO.MatrixInfo.col_set[j];
					}
					else
					{
						ofs_running << setw(10) << "0" << setw(10) << "0";
					}
					ofs_running << setw(10) << Hloc_fixed[i * ParaO.ncol + j] << endl;
				}
			}
		}
		for(int i=0; i<ParaO.nrow; i++)
		{
			for(int j=0; j<ParaO.ncol; j++)
			{
				double v = Hloc_fixed[i*ParaO.ncol+j];
				if( abs(v) > 1.0e-5 )
				{
					ofs_running << setw(15) << v;
				}
				else
				{
					ofs_running << setw(15) << "0";
				}
			}
			ofs_running << endl;
		}
		*/

		if(!BFIELD)
		{
			cout << " print Hloc_fixed" << endl;
			for(int i=0; i<NLOCAL; ++i)
			{
				for(int j=0; j<NLOCAL; ++j)
				{
					double v = Hloc_fixed[i*ParaO.ncol+j];
					if( abs(v) > 1.0e-8)
					{
						cout << setw(15) << v;
					}
					else
					{
						cout << setw(15) << "0";
					}
				}//end j
				cout << endl;
			}//end i
		}
	}
	if (mtype=='H')
	{

		if(!BFIELD)
		{
			cout << " print Hloc" << endl;
			for(int i=0; i<NLOCAL; ++i)
			{
				for(int j=0; j<NLOCAL; ++j)
				{
					double v = Hloc[i*ParaO.ncol+j];
					if( abs(v) > 1.0e-8)
					{
						cout << setw(15) << v;
					}
					else
					{
						cout << setw(15) << "0";
					}
				}//end j
				cout << endl;
			}//end i
		}

	/*
		ofs_running << " " << setw(10) << "LocalRow" << setw(10) << "LocalCol"
			<< setw(10) << "GlobalRow" << setw(10) << "GloablCol"
			<< setw(10) << "Hloc" << endl;

		ofs_running << "\n Hmatrix" << endl;
		ofs_running << setprecision(8) << endl;
		for (int i=0; i<ParaO.nrow; i++)
		{
			ofs_running << " ";
			for (int j=0; j<ParaO.ncol; j++)
			{
				double a = Hloc[ i * ParaO.ncol + j];
				if ( abs( a ) > 1.0e-5 )
				{
				//	ofs_running << " " << setw(10) << i << setw(10) << j
				//		<< setw(10) << ParaO.MatrixInfo.row_set[i] << setw(10) << ParaO.MatrixInfo.col_set[j]
				//		<< setw(10) << a << endl;
					
					//
					ofs_running << setw(15) << a;
					//ofs_running << setw(5) << i << setw(5) << j << setw(15) << a << endl;

				}
				else
				{
					//
					ofs_running << setw(15) << "0";
				}
			}
			//
			ofs_running << endl;
		}
		*/
	}

	return;
}

// becareful! Update Hloc, we add new members to it.
void LCAO_Matrix::update_Hloc(void)
{
	for (int i=0; i<ParaO.nloc; i++)
	{
		Hloc[i] += Hloc_fixed[i];
	}
	return;
}

void LCAO_Matrix::update_Hloc2(void)
{
	for (int i=0; i<ParaO.nloc; i++)
	{
		Hloc2[i] += Hloc_fixed2[i];
	}
	return;
}


void LCAO_Matrix::output_HSk(const char &mtype, string &fn)
{
	TITLE("LCAO_Matrix","output_HSk");
	stringstream ss;
	ss << global_out_dir << fn;
	ofstream ofs(ss.str().c_str());
	ofs << NLOCAL << endl;
	for(int i=0; i<NLOCAL; i++)
	{
		for(int j=0; j<NLOCAL; j++)
		{	
			const int index = i * NLOCAL + j;
			if(mtype=='S') ofs << Sloc2[index].real() << " " << Sloc2[index].imag() << endl;
			else if(mtype=='T') ofs << Hloc_fixed2[index].real() << " " << Hloc_fixed2[index].imag() << endl;
			else if(mtype=='H') ofs << Hloc2[index].real() << " " << Hloc2[index].imag() << endl;
		}
	}
	ofs.close();
	return;
}
