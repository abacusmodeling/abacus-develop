#include "LCAO_matrix.h"
#include "global_fp.h"

LCAO_Matrix::LCAO_Matrix()
{
	// for gamma_only
	Sloc = new double[1];
	Hloc_fixed = new double[1];
	Hloc = new double[1];
	Sdiag = new double[1];

	// for many k points
	Sloc2 = new std::complex<double>[1];
	Hloc_fixed2 = new std::complex<double>[1];
	Hloc2 = new std::complex<double>[1];
	Sdiag2 = new std::complex<double>[1];
}

LCAO_Matrix::~LCAO_Matrix()
{
 	// delete matrix for gamma_only.
    delete[] Sloc;
    delete[] Hloc_fixed;
    delete[] Hloc;
    delete[] Sdiag;

    // delete matrix for many k points
    delete[] Sloc2;
    delete[] Hloc_fixed2;
    delete[] Hloc2;	
    delete[] Sdiag2;
}


void LCAO_Matrix::divide_HS_in_frag(const bool isGamma, Parallel_Orbitals &po)
{
	TITLE("LCAO_Matrix","divide_HS_in_frag");

	GlobalV::ofs_running << "\n SETUP THE DIVISION OF H/S MATRIX" << std::endl;
	
	// (1) calculate nrow, ncol, nloc.
	if (GlobalV::KS_SOLVER=="genelpa" || GlobalV::KS_SOLVER=="hpseps" || GlobalV::KS_SOLVER=="scalpack" 
		|| GlobalV::KS_SOLVER=="selinv" || GlobalV::KS_SOLVER=="scalapack_gvx")
	{
		GlobalV::ofs_running << " divide the H&S matrix using 2D block algorithms." << std::endl;
#ifdef __MPI
		// storage form of H and S matrices on each processor
		// is determined in 'divide_HS_2d' subroutine
		po.divide_HS_2d(DIAG_WORLD);
#else
		WARNING_QUIT("LCAO_Matrix::init","diago method is not ready.");
#endif
	}
	else
	{
		// the full matrix
		po.nloc = GlobalV::NLOCAL * GlobalV::NLOCAL;
	}

	// (2) std::set the trace, then we can calculate the nnr.
	// for 2d: calculate po.nloc first, then trace_loc_row and trace_loc_col
	// for O(N): calculate the three together.
	po.set_trace();

	// (3) allocate for S, H_fixed, H, and S_diag
	if(isGamma)
	{
		allocate_HS_gamma(po.nloc);
	}
	else
	{
		allocate_HS_k(po.nloc);
	}

	return;
}


void LCAO_Matrix::allocate_HS_gamma(const long &nloc)
{
	TITLE("LCAO_Matrix","allocate_HS_gamma");

	OUT(GlobalV::ofs_running,"nloc",nloc);
	if(nloc==0) return; //mohan fix bug 2012-05-25

	// because we initilize in the constructor function
	// with dimension '1', so here we reconstruct these
	// matrices
	delete[] Sloc;
	delete[] Hloc_fixed;
	delete[] Hloc;
	delete[] Sdiag;

	this->Sloc = new double[nloc];
	this->Hloc_fixed = new double[nloc];
	this->Hloc = new double[nloc];
	this->Sdiag = new double[nloc];

	ZEROS(Sloc,nloc);
	ZEROS(Hloc_fixed,nloc);
	ZEROS(Hloc,nloc);
	ZEROS(Sdiag,nloc); // mohan add 2021-01-30

	return;
}


void LCAO_Matrix::allocate_HS_k(const long &nloc)
{
	TITLE("LCAO_Matrix","allocate_HS_k");

	OUT(GlobalV::ofs_running,"nloc",nloc);
	if(nloc==0) return; //mohan fix bug 2012-05-25

	// because we initilize in the constructor function
	// with dimension '1', so here we reconstruct these
	// matrices
	delete[] Sloc2;
	delete[] Hloc_fixed2;
	delete[] Hloc2;
	delete[] Sdiag2;

	this->Sloc2 = new std::complex<double>[nloc];
	this->Hloc_fixed2 = new std::complex<double>[nloc];
	this->Hloc2 = new std::complex<double>[nloc];
	this->Sdiag2 = new std::complex<double>[nloc];

	ZEROS(Sloc2,nloc);
	ZEROS(Hloc_fixed2,nloc);
	ZEROS(Hloc2,nloc);
	
	return;
}

void LCAO_Matrix::allocate_HS_R(const int &nnR)
{
	if(GlobalV::NSPIN!=4)
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
	}
	else
	{
		delete[] HlocR_soc;
		delete[] SlocR_soc;
		delete[] Hloc_fixedR_soc;

		this->HlocR_soc = new std::complex<double>[nnR];
		this->SlocR_soc = new std::complex<double>[nnR];
		this->Hloc_fixedR_soc = new std::complex<double>[nnR];
		
		ZEROS(HlocR_soc, nnR);
		ZEROS(SlocR_soc, nnR);
		ZEROS(Hloc_fixedR_soc, nnR);
		
	}

	return;
}

void LCAO_Matrix::set_HSgamma(const int &iw1_all, const int &iw2_all, const double &v, const char &dtype)
{
    // use iw1_all and iw2_all to std::set Hloc
    // becareful! The ir and ic may < 0!!!!!!!!!!!!!!!!
    const int ir = GlobalC::ParaO.trace_loc_row[ iw1_all ];
    const int ic = GlobalC::ParaO.trace_loc_col[ iw2_all ];
    //const int index = ir * GlobalC::ParaO.ncol + ic;
	long index;
	if(GlobalV::KS_SOLVER=="genelpa" || GlobalV::KS_SOLVER=="scalapack_gvx")  // save the matrix as column major format
	{
		index=ic*GlobalC::ParaO.nrow+ir;
	}
	else
	{
		index=ir*GlobalC::ParaO.ncol+ic;
  	}
   
   	if( index >= GlobalC::ParaO.nloc)
	{
		std::cout << " iw1_all = " << iw1_all << std::endl;
		std::cout << " iw2_all = " << iw2_all << std::endl;
		std::cout << " ir = " << ir << std::endl;
		std::cout << " ic = " << ic << std::endl;
		std::cout << " index = " << index << std::endl;
		std::cout << " GlobalC::ParaO.nloc = " << GlobalC::ParaO.nloc << std::endl;
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

void LCAO_Matrix::set_HSk(const int &iw1_all, const int &iw2_all, const std::complex<double> &v, const char &dtype, const int spin)
{
    // use iw1_all and iw2_all to std::set Hloc
    // becareful! The ir and ic may < 0!!!!!!!!!!!!!!!!
    const int ir = GlobalC::ParaO.trace_loc_row[ iw1_all ];
    const int ic = GlobalC::ParaO.trace_loc_col[ iw2_all ];
    //const int index = ir * GlobalC::ParaO.ncol + ic;
	long index;
	if(GlobalV::KS_SOLVER=="genelpa" || GlobalV::KS_SOLVER=="scalapack_gvx")  // save the matrix as column major format
	{
		index=ic*GlobalC::ParaO.nrow+ir;
	}
	else
	{
		index=ir*GlobalC::ParaO.ncol+ic;
  	}
    assert(index < GlobalC::ParaO.nloc);
	if (dtype=='S')//overlap Hamiltonian.
	{
		this->Sloc2[index] += v;
	}
	else if (dtype=='T' || dtype=='N')// kinetic and nonlocal Hamiltonian.
	{
		this->Hloc_fixed2[index] += v; // because kinetic and nonlocal Hamiltonian matrices are already block-cycle staraged after caculated in lcao_nnr.cpp
                                      // this statement will not be used.
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
    // use iw1_all and iw2_all to std::set Hloc
    // becareful! The ir and ic may < 0!!!!!!!!!!!!!!!!
    const int ir = GlobalC::ParaO.trace_loc_row[ iw1_all ];
    const int ic = GlobalC::ParaO.trace_loc_col[ iw2_all ];
    const long index = ir * GlobalC::ParaO.ncol + ic;
	
	if( index >= GlobalC::ParaO.nloc)
	{
		std::cout << " iw1_all = " << iw1_all << std::endl;
		std::cout << " iw2_all = " << iw2_all << std::endl;
		std::cout << " ir = " << ir << std::endl;
		std::cout << " ic = " << ic << std::endl;
		std::cout << " index = " << index << std::endl;
		std::cout << " GlobalC::ParaO.nloc = " << GlobalC::ParaO.nloc << std::endl;
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

void LCAO_Matrix::set_stress
(
    const int &iw1_all,
    const int &iw2_all,
    const double& vx,
    const double& vy,
    const double& vz,
    const char &dtype,
    const Vector3<double> &dtau)
{
    // use iw1_all and iw2_all to std::set Hloc
    // becareful! The ir and ic may < 0!!!!!!!!!!!!!!!!
    const int ir = GlobalC::ParaO.trace_loc_row[ iw1_all ];
    const int ic = GlobalC::ParaO.trace_loc_col[ iw2_all ];
    const long index = ir * GlobalC::ParaO.ncol + ic;

	if( index >= GlobalC::ParaO.nloc)
	{
		std::cout << " iw1_all = " << iw1_all << std::endl;
		std::cout << " iw2_all = " << iw2_all << std::endl;
		std::cout << " ir = " << ir << std::endl;
		std::cout << " ic = " << ic << std::endl;
		std::cout << " index = " << index << std::endl;
		std::cout << " GlobalC::ParaO.nloc = " << GlobalC::ParaO.nloc << std::endl;
		WARNING_QUIT("LCAO_Matrix","set_stress");
	}

	if (dtype == 'S')
	{
		this->DSloc_11[index] += vx * dtau.x;
		this->DSloc_12[index] += vx * dtau.y;
		this->DSloc_13[index] += vx * dtau.z;
		this->DSloc_22[index] += vy * dtau.y;
		this->DSloc_23[index] += vy * dtau.z;
		this->DSloc_33[index] += vz * dtau.z;
	}
	else if (dtype == 'T')
	{
		// notice, the sign is '-', minus.
		this->DHloc_fixed_11[index] -= vx * dtau.x;
		this->DHloc_fixed_12[index] -= vx * dtau.y;
		this->DHloc_fixed_13[index] -= vx * dtau.z;
		this->DHloc_fixed_22[index] -= vy * dtau.y;
		this->DHloc_fixed_23[index] -= vy * dtau.z;
		this->DHloc_fixed_33[index] -= vz * dtau.z;
	}
	else if (dtype == 'N')
	{
		this->DHloc_fixed_11[index] += vx * dtau.x;
		this->DHloc_fixed_12[index] += vx * dtau.y;
		this->DHloc_fixed_13[index] += vx * dtau.z;
		this->DHloc_fixed_22[index] += vy * dtau.y;
		this->DHloc_fixed_23[index] += vy * dtau.z;
		this->DHloc_fixed_33[index] += vz * dtau.z;
	}

	return;
}

void LCAO_Matrix::zeros_HSgamma(const char &mtype)
{
	if (mtype=='S') ZEROS(Sloc,GlobalC::ParaO.nloc);
	else if (mtype=='T') ZEROS(Hloc_fixed,GlobalC::ParaO.nloc);
	else if (mtype=='H') ZEROS(Hloc,GlobalC::ParaO.nloc);
	return;
}

void LCAO_Matrix::zeros_HSk(const char &mtype)
{
	if (mtype=='S') ZEROS(Sloc2,GlobalC::ParaO.nloc);
	else if (mtype=='T') ZEROS(Hloc_fixed2,GlobalC::ParaO.nloc);
	else if (mtype=='H') ZEROS(Hloc2,GlobalC::ParaO.nloc);
	return;
}

void LCAO_Matrix::zeros_HSR(const char &mtype, const int &nnr)
{
	if(GlobalV::NSPIN!=4)
	{
		if (mtype=='S') ZEROS(SlocR, nnr);
		else if (mtype=='T') ZEROS(Hloc_fixedR, nnr);
		else if (mtype=='H') ZEROS(HlocR, nnr);
	}
	else
	{
		if (mtype=='H') ZEROS(this->HlocR_soc, nnr);
		else if (mtype=='S') ZEROS(this->SlocR_soc, nnr);
		else if (mtype=='T') ZEROS(this->Hloc_fixedR_soc, nnr);
	}
	return;
}

// Peize Lin add vtype='A' 2018-11-30
void LCAO_Matrix::print_HSk(const char &mtype, const char &vtype, const double &accuracy, std::ostream &os)
{
	TITLE("LCAO_Matrix","print_HSk");
	if(mtype=='S') os << "Sloc2 matrix" << std::endl;
	else if(mtype=='T') os << "Hloc_fixed2 matrix" << std::endl;
	else if(mtype=='H') os << "Hloc2 matrix" << std::endl;
	else
	{
		WARNING_QUIT("LCAO_Matrix::print_HSk","Check input parameter: mtype.");
	}

	if(vtype=='C') os << " Output norm."  << std::endl;
	else if(vtype=='R') os << " Output real part."  << std::endl;
	else if(vtype=='I') os << " Output imag part."  << std::endl;
	else if(vtype=='A') os << " Output std::complex." << std::endl;


	os << setprecision(8) << std::endl;
	for(int i=0; i<GlobalC::ParaO.nrow; i++)
	{
		os << " " ;
		for(int j=0; j<GlobalC::ParaO.ncol; j++)
		{
			const int index = i * GlobalC::ParaO.ncol + j;
			if(vtype=='A')
			{
				std::complex<double> v;
				if(mtype=='S')	v = Sloc2[index];
				else if(mtype=='T') v = Hloc_fixed2[index];
				else if(mtype=='H') v = Hloc2[index];
				auto threshold = [accuracy]( const double v ){ return abs(v)>accuracy ? v : 0.0; };
				os << '(' << threshold(v.real()) << ',' << threshold(v.imag()) << "\t";
			}
			else
			{
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
	//				os << setw(15) << v;
					os << v << "\t";
				}
				else
				{
	//				os << setw(15) << "0"; 
					os << "0" << "\t"; 
				}
			}
		}
		os << std::endl;
	}
	os << std::endl;
	os << setprecision(6) << std::endl;
	return;
}


void LCAO_Matrix::print_HSgamma(const char &mtype, std::ostream &os)
{
	TITLE("Parallel_Orbitals","print_HSgamma");

	GlobalV::ofs_running << " " << mtype << " matrix" << std::endl;
	GlobalV::ofs_running << " nrow=" << GlobalC::ParaO.nrow << std::endl;
	GlobalV::ofs_running << " ncol=" << GlobalC::ParaO.ncol << std::endl;
	GlobalV::ofs_running << " element number = " << GlobalC::ParaO.ncol << std::endl;

	if (mtype=='S')
	{
		os << setprecision(8);
		os << " print Sloc" << std::endl;
		for(int i=0; i<GlobalV::NLOCAL; ++i)
		{
			for(int j=0; j<GlobalV::NLOCAL; ++j)
			{
				double v = Sloc[i*GlobalC::ParaO.ncol+j];
				if( abs(v) > 1.0e-8)
				{
					os << setw(15) << v;
				}
				else
				{
					os << setw(15) << "0";
				}
			}//end j
			os << std::endl;
		}//end i
	}
	if (mtype=='T')
	{
		os << " print Hloc_fixed" << std::endl;
		for(int i=0; i<GlobalV::NLOCAL; ++i)
		{
			for(int j=0; j<GlobalV::NLOCAL; ++j)
			{
				double v = Hloc_fixed[i*GlobalC::ParaO.ncol+j];
				if( abs(v) > 1.0e-8)
				{
					os << setw(15) << v;
				}
				else
				{
					os << setw(15) << "0";
				}
			}//end j
			os << std::endl;
		}//end i
	}
	if (mtype=='H')
	{
		os << " print Hloc" << std::endl;
		for(int i=0; i<GlobalV::NLOCAL; ++i)
		{
			for(int j=0; j<GlobalV::NLOCAL; ++j)
			{
				double v = Hloc[i*GlobalC::ParaO.ncol+j];
				if( abs(v) > 1.0e-8)
				{
					os << setw(15) << v;
				}
				else
				{
					os << setw(15) << "0";
				}
			}//end j
			os << std::endl;
		}//end i
	}

	return;
}

// becareful! Update Hloc, we add new members to it.
void LCAO_Matrix::update_Hloc(void)
{
	for (long i=0; i<GlobalC::ParaO.nloc; i++)
	{
		Hloc[i] += Hloc_fixed[i];
	}
	return;
}

void LCAO_Matrix::update_Hloc2(void)
{
	for (long i=0; i<GlobalC::ParaO.nloc; i++)
	{
		Hloc2[i] += Hloc_fixed2[i];
	}
	return;
}


void LCAO_Matrix::output_HSk(const char &mtype, std::string &fn)
{
	TITLE("LCAO_Matrix","output_HSk");
	std::stringstream ss;
	ss << GlobalV::global_out_dir << fn;
	std::ofstream ofs(ss.str().c_str());
	ofs << GlobalV::NLOCAL << std::endl;
	for(int i=0; i<GlobalV::NLOCAL; i++)
	{
		for(int j=0; j<GlobalV::NLOCAL; j++)
		{	
			const int index = i * GlobalV::NLOCAL + j;
			if(mtype=='S') ofs << Sloc2[index].real() << " " << Sloc2[index].imag() << std::endl;
			else if(mtype=='T') ofs << Hloc_fixed2[index].real() << " " << Hloc_fixed2[index].imag() << std::endl;
			else if(mtype=='H') ofs << Hloc2[index].real() << " " << Hloc2[index].imag() << std::endl;
		}
	}
	ofs.close();
	return;
}

void LCAO_Matrix::allocate_Hloc_fixedR_tr(void)
{
    TITLE("LCAO_Matrix","allocate_Hloc_fixedR_tr");

    //int R_x = 10;
    //int R_y = 10;
    //int R_z = 10;
    int R_x = GlobalC::GridD.getCellX();
    int R_y = GlobalC::GridD.getCellY();
    int R_z = GlobalC::GridD.getCellZ();

    if(GlobalV::NSPIN!=4)
    {
        Hloc_fixedR_tr = new double***[R_x];
        //HR_tr = new double***[R_x];
        //SlocR_tr = new double***[R_x];
        for(int ix=0; ix<R_x; ix++)
        {
            Hloc_fixedR_tr[ix] = new double**[R_y];
            //HR_tr[ix] = new double**[R_y];
            //SlocR_tr[ix] = new double**[R_y];
            for(int iy=0; iy<R_y; iy++)
            {
                Hloc_fixedR_tr[ix][iy] = new double*[R_z];
                //HR_tr[ix][iy] = new double*[R_z];
                //SlocR_tr[ix][iy] = new double*[R_z];
                for(int iz=0; iz<R_z; iz++)
                {
                    Hloc_fixedR_tr[ix][iy][iz] = new double[GlobalC::ParaO.nloc];
                    //HR_tr[ix][iy][iz] = new double[GlobalC::ParaO.nloc];
                    //SlocR_tr[ix][iy][iz] = new double[GlobalC::ParaO.nloc];
                    ZEROS(Hloc_fixedR_tr[ix][iy][iz], GlobalC::ParaO.nloc);
                    //ZEROS(HR_tr[ix][iy][iz], GlobalC::ParaO.nloc);
                    //ZEROS(SlocR_tr[ix][iy][iz], GlobalC::ParaO.nloc);
                }
            }
        }
    }
    else
    {
        Hloc_fixedR_tr_soc = new std::complex<double>***[R_x];
        //HR_tr = new double***[R_x];
        //SlocR_tr = new double***[R_x];
        for(int ix=0; ix<R_x; ix++)
        {
            Hloc_fixedR_tr_soc[ix] = new std::complex<double>**[R_y];
            //HR_tr[ix] = new double**[R_y];
            //SlocR_tr[ix] = new double**[R_y];
            for(int iy=0; iy<R_y; iy++)
            {
                Hloc_fixedR_tr_soc[ix][iy] = new std::complex<double>*[R_z];
                //HR_tr[ix][iy] = new double*[R_z];
                //SlocR_tr[ix][iy] = new double*[R_z];
                for(int iz=0; iz<R_z; iz++)
                {
                    Hloc_fixedR_tr_soc[ix][iy][iz] = new std::complex<double>[GlobalC::ParaO.nloc];
                    //HR_tr[ix][iy][iz] = new double[GlobalC::ParaO.nloc];
                    //SlocR_tr[ix][iy][iz] = new double[GlobalC::ParaO.nloc];
                    ZEROS(Hloc_fixedR_tr_soc[ix][iy][iz], GlobalC::ParaO.nloc);
                    //ZEROS(HR_tr[ix][iy][iz], GlobalC::ParaO.nloc);
                    //ZEROS(SlocR_tr[ix][iy][iz], GlobalC::ParaO.nloc);
                }
            }
        }
    }
//std::cout<<"R_x: "<<R_x<<std::endl;
//std::cout<<"R_y: "<<R_y<<std::endl;
//std::cout<<"R_z: "<<R_z<<std::endl;
//std::cout<<"GlobalC::ParaO.nloc: "<<GlobalC::ParaO.nloc<<std::endl;
//std::cout<<"SlocR_tr 1-3-3-27: "<<SlocR_tr[1][3][3][27]<<std::endl;
//std::cout<<"Hloc_fixedR_tr 1-3-3-27: "<<Hloc_fixedR_tr[1][3][3][27]<<std::endl;

    return;
}

void LCAO_Matrix::allocate_HR_tr(void)
{
    TITLE("LCAO_Matrix","allocate_HR_tr");

    //int R_x = 10;
    //int R_y = 10;
    //int R_z = 10;
    int R_x = GlobalC::GridD.getCellX();
    int R_y = GlobalC::GridD.getCellY();
    int R_z = GlobalC::GridD.getCellZ();

    if(GlobalV::NSPIN!=4)
    {
        HR_tr = new double***[R_x];
        for(int ix=0; ix<R_x; ix++)
        {
            HR_tr[ix] = new double**[R_y];
            for(int iy=0; iy<R_y; iy++)
            {
                HR_tr[ix][iy] = new double*[R_z];
                for(int iz=0; iz<R_z; iz++)
                {
                    HR_tr[ix][iy][iz] = new double[GlobalC::ParaO.nloc];
                    ZEROS(HR_tr[ix][iy][iz], GlobalC::ParaO.nloc);
                }
            }
        }
    }
    else
    {
        HR_tr_soc = new std::complex<double>***[R_x];
        for(int ix=0; ix<R_x; ix++)
        {
            HR_tr_soc[ix] = new std::complex<double>**[R_y];
            for(int iy=0; iy<R_y; iy++)
            {
                HR_tr_soc[ix][iy] = new std::complex<double>*[R_z];
                for(int iz=0; iz<R_z; iz++)
                {
                    HR_tr_soc[ix][iy][iz] = new std::complex<double>[GlobalC::ParaO.nloc];
                    ZEROS(HR_tr_soc[ix][iy][iz], GlobalC::ParaO.nloc);
                }
            }
        }
    }

    return;
}

void LCAO_Matrix::allocate_SlocR_tr(void)
{
    TITLE("LCAO_Matrix","allocate_SlocR_tr");

    //int R_x = 10;
    //int R_y = 10;
    //int R_z = 10;
    int R_x = GlobalC::GridD.getCellX();
    int R_y = GlobalC::GridD.getCellY();
    int R_z = GlobalC::GridD.getCellZ();

    if(GlobalV::NSPIN!=4)
    {
        SlocR_tr = new double***[R_x];
        for(int ix=0; ix<R_x; ix++)
        {
            SlocR_tr[ix] = new double**[R_y];
            for(int iy=0; iy<R_y; iy++)
            {
                SlocR_tr[ix][iy] = new double*[R_z];
                for(int iz=0; iz<R_z; iz++)
                {
                    SlocR_tr[ix][iy][iz] = new double[GlobalC::ParaO.nloc];
                    ZEROS(SlocR_tr[ix][iy][iz], GlobalC::ParaO.nloc);
                }
            }
        }
    }
    else
    {
        SlocR_tr_soc = new std::complex<double>***[R_x];
        for(int ix=0; ix<R_x; ix++)
        {
            SlocR_tr_soc[ix] = new std::complex<double>**[R_y];
            for(int iy=0; iy<R_y; iy++)
            {
                SlocR_tr_soc[ix][iy] = new std::complex<double>*[R_z];
                for(int iz=0; iz<R_z; iz++)
                {
                    SlocR_tr_soc[ix][iy][iz] = new std::complex<double>[GlobalC::ParaO.nloc];
                    ZEROS(SlocR_tr_soc[ix][iy][iz], GlobalC::ParaO.nloc);
                }
            }
        }
    }

    return;
}

void LCAO_Matrix::destroy_Hloc_fixedR_tr(void)
{
    TITLE("LCAO_Matrix","destroy_Hloc_fixed2_R");

    //int R_x = 10;
    //int R_y = 10;
    //int R_z = 10;
    int R_x = GlobalC::GridD.getCellX();
    int R_y = GlobalC::GridD.getCellY();
    int R_z = GlobalC::GridD.getCellZ();

    if(GlobalV::NSPIN!=4)
    {
        for(int ix=0; ix<R_x; ix++)
        {
            for(int iy=0; iy<R_y; iy++)
            {
                for(int iz=0; iz<R_z; iz++)
                {
                    delete[] Hloc_fixedR_tr[ix][iy][iz];
                    delete[] HR_tr[ix][iy][iz];
                    delete[] SlocR_tr[ix][iy][iz];
                }
                delete[] Hloc_fixedR_tr[ix][iy];
                delete[] HR_tr[ix][iy];
                delete[] SlocR_tr[ix][iy];
            }
            delete[] Hloc_fixedR_tr[ix];
            delete[] HR_tr[ix];
            delete[] SlocR_tr[ix];
        }
        delete[] Hloc_fixedR_tr;
        delete[] HR_tr;
        delete[] SlocR_tr;
    }
    else
    {
        for(int ix=0; ix<R_x; ix++)
        {
            for(int iy=0; iy<R_y; iy++)
            {
                for(int iz=0; iz<R_z; iz++)
                {
                    delete[] Hloc_fixedR_tr_soc[ix][iy][iz];
                    delete[] HR_tr_soc[ix][iy][iz];
                    delete[] SlocR_tr_soc[ix][iy][iz];
                }
                delete[] Hloc_fixedR_tr_soc[ix][iy];
                delete[] HR_tr_soc[ix][iy];
                delete[] SlocR_tr_soc[ix][iy];
            }
            delete[] Hloc_fixedR_tr_soc[ix];
            delete[] HR_tr_soc[ix];
            delete[] SlocR_tr_soc[ix];
        }
        delete[] Hloc_fixedR_tr_soc;
        delete[] HR_tr_soc;
        delete[] SlocR_tr_soc;
    }

    return;
}

void LCAO_Matrix::set_HR_tr(const int &Rx, const int &Ry, const int &Rz, const int &iw1_all, const int &iw2_all, const double &v)
{
    const int ir = GlobalC::ParaO.trace_loc_row[ iw1_all ];
    const int ic = GlobalC::ParaO.trace_loc_col[ iw2_all ];

//std::cout<<"ir: "<<ir<<std::endl;
//std::cout<<"ic: "<<ic<<std::endl;
    long index;
    if(GlobalV::KS_SOLVER=="genelpa" || GlobalV::KS_SOLVER=="scalapack_gvx")
    {
        index=ic*GlobalC::ParaO.nrow+ir;
//std::cout<<"index: "<<index<<std::endl;
    }
    else
    {
        index=ir*GlobalC::ParaO.ncol+ic;
//std::cout<<"index: "<<index<<std::endl;
    }

//std::cout<<"GlobalC::ParaO.nloc: "<<GlobalC::ParaO.nloc<<std::endl;
    assert(index < GlobalC::ParaO.nloc);
//std::cout<<"Rx: "<<Rx<<std::endl;
//std::cout<<"Ry: "<<Ry<<std::endl;
//std::cout<<"Rz: "<<Rz<<std::endl;
//std::cout<<"Hloc_fixedR_tr: "<<Hloc_fixedR_tr[Rx][Ry][Rz][index]<<std::endl;
//std::cout<<"v: "<<v<<std::endl;
    HR_tr[Rx][Ry][Rz][index] = Hloc_fixedR_tr[Rx][Ry][Rz][index] + v; 
    //HR_tr[Rx][Ry][Rz][index] = Hloc_fixedR_tr[Rx][Ry][Rz][index]; 
    //HR_tr[Rx][Ry][Rz][index] = v; 
    //HR_tr[Rx][Ry][Rz][index] = index; 

    return;
}

//LiuXh add 2019-07-16
void LCAO_Matrix::set_HR_tr_soc(const int &Rx, const int &Ry, const int &Rz, const int &iw1_all, const int &iw2_all, const std::complex<double> &v)
{
    const int ir = GlobalC::ParaO.trace_loc_row[ iw1_all ];
    const int ic = GlobalC::ParaO.trace_loc_col[ iw2_all ];

//std::cout<<"ir: "<<ir<<std::endl;
//std::cout<<"ic: "<<ic<<std::endl;
    long index;
    if(GlobalV::KS_SOLVER=="genelpa" || GlobalV::KS_SOLVER=="scalapack_gvx")
    {
        index=ic*GlobalC::ParaO.nrow+ir;
//std::cout<<"index: "<<index<<std::endl;
    }
    else
    {
        index=ir*GlobalC::ParaO.ncol+ic;
//std::cout<<"index: "<<index<<std::endl;
    }

//std::cout<<"GlobalC::ParaO.nloc: "<<GlobalC::ParaO.nloc<<std::endl;
    assert(index < GlobalC::ParaO.nloc);
//std::cout<<"Rx: "<<Rx<<std::endl;
//std::cout<<"Ry: "<<Ry<<std::endl;
//std::cout<<"Rz: "<<Rz<<std::endl;
//std::cout<<"Hloc_fixedR_tr: "<<Hloc_fixedR_tr[Rx][Ry][Rz][index]<<std::endl;
//std::cout<<"v: "<<v<<std::endl;
    HR_tr_soc[Rx][Ry][Rz][index] = Hloc_fixedR_tr_soc[Rx][Ry][Rz][index] + v; 
    //HR_tr[Rx][Ry][Rz][index] = Hloc_fixedR_tr[Rx][Ry][Rz][index]; 
    //HR_tr[Rx][Ry][Rz][index] = v; 
    //HR_tr[Rx][Ry][Rz][index] = index; 

    return;
}

void LCAO_Matrix::allocate_HS_R_sparse(void)
{
	TITLE("LCAO_Matrix","allocate_HS_R_sparse");

	int R_x = GlobalC::GridD.getCellX();
    int R_y = GlobalC::GridD.getCellY();
    int R_z = GlobalC::GridD.getCellZ();

	if (GlobalV::NSPIN != 4)
	{
		HR_sparse = new std::map<size_t, std::map<size_t, double>> **[R_x];
		SR_sparse = new std::map<size_t, std::map<size_t, double>> **[R_x];
		for (int ix = 0; ix < R_x; ++ix)
		{
			HR_sparse[ix] = new std::map<size_t, std::map<size_t, double>> *[R_y];
			SR_sparse[ix] = new std::map<size_t, std::map<size_t, double>> *[R_y];
			for (int iy = 0; iy < R_y; ++iy)
			{
				HR_sparse[ix][iy] = new std::map<size_t, std::map<size_t, double>>[R_z];
				SR_sparse[ix][iy] = new std::map<size_t, std::map<size_t, double>>[R_z];
			}
		}
	}
	else
	{
		HR_soc_sparse = new std::map<size_t, std::map<size_t, std::complex<double>>> **[R_x];
		SR_soc_sparse = new std::map<size_t, std::map<size_t, std::complex<double>>> **[R_x];
		for (int ix = 0; ix < R_x; ++ix)
		{
			HR_soc_sparse[ix] = new std::map<size_t, std::map<size_t, std::complex<double>>> *[R_y];
			SR_soc_sparse[ix] = new std::map<size_t, std::map<size_t, std::complex<double>>> *[R_y];
			for (int iy = 0; iy < R_y; ++iy)
			{
				HR_soc_sparse[ix][iy] = new std::map<size_t, std::map<size_t, std::complex<double>>>[R_z];
				SR_soc_sparse[ix][iy] = new std::map<size_t, std::map<size_t, std::complex<double>>>[R_z];
			}
		}
	}

	return;
}

void LCAO_Matrix::destroy_HS_R_sparse(void)
{
	TITLE("LCAO_Matrix","destroy_HS_R_sparse");

	int R_x = GlobalC::GridD.getCellX();
    int R_y = GlobalC::GridD.getCellY();

	if (GlobalV::NSPIN != 4)
	{
		for (int ix = 0; ix < R_x; ++ix)
		{
			for (int iy = 0; iy < R_y; ++iy)
			{
				delete[] HR_sparse[ix][iy];
				delete[] SR_sparse[ix][iy];
			}
			delete[] HR_sparse[ix];
			delete[] SR_sparse[ix];
		}
		delete[] HR_sparse;
		delete[] SR_sparse;
		HR_sparse = nullptr;
		SR_sparse = nullptr;
	}
	else
	{
		for (int ix = 0; ix < R_x; ++ix)
		{
			for (int iy = 0; iy < R_y; ++iy)
			{
				delete[] HR_soc_sparse[ix][iy];
				delete[] SR_soc_sparse[ix][iy];
			}
			delete[] HR_soc_sparse[ix];
			delete[] SR_soc_sparse[ix];
		}
		delete[] HR_soc_sparse;
		delete[] SR_soc_sparse;
		HR_soc_sparse = nullptr;
		SR_soc_sparse = nullptr;
	}

	return;
}