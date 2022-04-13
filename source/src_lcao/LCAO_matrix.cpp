#include "LCAO_matrix.h"
#include "global_fp.h"
#include "../src_pw/global.h"
#ifdef __DEEPKS
#include "../module_deepks/LCAO_deepks.h"
#endif

LCAO_Matrix::LCAO_Matrix()
{
}

LCAO_Matrix::~LCAO_Matrix()
{
}


void LCAO_Matrix::divide_HS_in_frag(const bool isGamma, Parallel_Orbitals &pv)
{
    ModuleBase::TITLE("LCAO_Matrix","divide_HS_in_frag");

    //(1), (2): set up matrix division have been moved into ORB_control
    //just pass `ParaV` as pointer is enough
    this->ParaV = &pv;
	// (3) allocate for S, H_fixed, H, and S_diag
	if(isGamma)
	{
		allocate_HS_gamma(this->ParaV->nloc);
	}
	else
	{
		allocate_HS_k(this->ParaV->nloc);
	}
#ifdef __DEEPKS
	//wenfei 2021-12-19
    //preparation for DeePKS

	if (GlobalV::deepks_out_labels || GlobalV::deepks_scf)
	{
        //allocate relevant data structures for calculating descriptors
        std::vector<int> na;
        na.resize(GlobalC::ucell.ntype);
        for(int it=0;it<GlobalC::ucell.ntype;it++)
        {
            na[it] = GlobalC::ucell.atoms[it].na;
        }
		GlobalC::ld.init(GlobalC::ORB,
            GlobalC::ucell.nat,
            GlobalC::ucell.ntype,
            na);
        if(GlobalV::deepks_scf)
        {
            if(isGamma)
            {
                GlobalC::ld.allocate_V_delta(GlobalC::ucell.nat, pv.nloc);
            }
            else
            {
                GlobalC::ld.allocate_V_delta(GlobalC::ucell.nat, pv.nloc,GlobalC::kv.nks);
            }
        }
	}
#endif
	return;
}


void LCAO_Matrix::allocate_HS_gamma(const long &nloc)
{
    ModuleBase::TITLE("LCAO_Matrix","allocate_HS_gamma");

    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"nloc",nloc);
    if(nloc==0) return; //mohan fix bug 2012-05-25

    // because we initilize in the constructor function
    // with dimension '1', so here we reconstruct these
    // matrices

    this->Sloc.resize(nloc);
    this->Hloc_fixed.resize(nloc);
    this->Hloc.resize(nloc);
    this->Sdiag.resize(nloc);

    ModuleBase::GlobalFunc::ZEROS(Sloc.data(),nloc);
    ModuleBase::GlobalFunc::ZEROS(Hloc_fixed.data(),nloc);
    ModuleBase::GlobalFunc::ZEROS(Hloc.data(),nloc);
    ModuleBase::GlobalFunc::ZEROS(Sdiag.data(),nloc); // mohan add 2021-01-30

    return;
}


void LCAO_Matrix::allocate_HS_k(const long &nloc)
{
    ModuleBase::TITLE("LCAO_Matrix","allocate_HS_k");

    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"nloc",nloc);
    if(nloc==0) return; //mohan fix bug 2012-05-25

    // because we initilize in the constructor function
    // with dimension '1', so here we reconstruct these
    // matrices
    this->Sloc2.resize(nloc);
    this->Hloc_fixed2.resize(nloc);
    this->Hloc2.resize(nloc);
    this->Sdiag2.resize(nloc);

    ModuleBase::GlobalFunc::ZEROS(Sloc2.data(),nloc);
    ModuleBase::GlobalFunc::ZEROS(Hloc_fixed2.data(),nloc);
    ModuleBase::GlobalFunc::ZEROS(Hloc2.data(),nloc);
    
    return;
}

void LCAO_Matrix::allocate_HS_R(const int &nnR)
{
    if(GlobalV::NSPIN!=4)
    {	
        this->SlocR.resize(nnR);
        this->Hloc_fixedR.resize(nnR);

        ModuleBase::GlobalFunc::ZEROS(SlocR.data(), nnR);
        ModuleBase::GlobalFunc::ZEROS(Hloc_fixedR.data(), nnR);
    }
    else
    {
        this->SlocR_soc.resize(nnR);
        this->Hloc_fixedR_soc.resize(nnR);
        
        ModuleBase::GlobalFunc::ZEROS(SlocR_soc.data(), nnR);
        ModuleBase::GlobalFunc::ZEROS(Hloc_fixedR_soc.data(), nnR);
        
    }

    return;
}

//------------------------------------------------------
// DESCRIPTION:
// set 'dtype' matrix element (iw1_all, iw2_all) with 
// an input value 'v'
//------------------------------------------------------
void LCAO_Matrix::set_HSgamma(
    const int &iw1_all, // index i for atomic orbital (row)
    const int &iw2_all, // index j for atomic orbital (column)
    const double &v, // value for matrix element (i,j) 
    const char &dtype, // type of the matrix
    double* HSloc) //input pointer for store the matrix
{
    // use iw1_all and iw2_all to set Hloc
    // becareful! The ir and ic may be < 0 !!!
    const int ir = this->ParaV->trace_loc_row[ iw1_all ];
    const int ic = this->ParaV->trace_loc_col[ iw2_all ];

    //const int index = ir * ParaO.ncol + ic;
    long index=0;

    // save the matrix as column major format
    if(GlobalV::KS_SOLVER=="genelpa" || GlobalV::KS_SOLVER=="scalapack_gvx")
    {
        index=ic*this->ParaV->nrow+ir;
    }
    else
    {
        index=ir*this->ParaV->ncol+ic;
      }
   
       if( index >= this->ParaV->nloc)
    {
        std::cout << " iw1_all = " << iw1_all << std::endl;
        std::cout << " iw2_all = " << iw2_all << std::endl;
        std::cout << " ir = " << ir << std::endl;
        std::cout << " ic = " << ic << std::endl;
        std::cout << " index = " << index << std::endl;
        std::cout << " this->ParaV->nloc = " << this->ParaV->nloc << std::endl;
        ModuleBase::WARNING_QUIT("LCAO_Matrix","set_HSgamma");
    }	 

    //-----------------------------------
    // dtype: type of the matrix.
    // S : S matrix element.
    // T : T matrix element.
    // N : nonlocal H matrix element.
    // L : local H matrix element.
    //-----------------------------------
    /*if (dtype=='S')
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
    }*/
    //using input pointer HSloc
    HSloc[index] += v;

    return;
}

void LCAO_Matrix::set_HSk(const int &iw1_all, const int &iw2_all, const std::complex<double> &v, const char &dtype, const int spin)
{
    // use iw1_all and iw2_all to set Hloc
    // becareful! The ir and ic may < 0!!!!!!!!!!!!!!!!
    const int ir = this->ParaV->trace_loc_row[ iw1_all ];
    const int ic = this->ParaV->trace_loc_col[ iw2_all ];
    //const int index = ir * this->ParaV->ncol + ic;
    long index;
    if(GlobalV::KS_SOLVER=="genelpa" || GlobalV::KS_SOLVER=="scalapack_gvx")  // save the matrix as column major format
    {
        index=ic*this->ParaV->nrow+ir;
    }
    else
    {
        index=ir*this->ParaV->ncol+ic;
      }
    assert(index < this->ParaV->nloc);
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
        ModuleBase::WARNING_QUIT("LCAO_Matrix","set_HSk");
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
    const int ir = this->ParaV->trace_loc_row[ iw1_all ];
    const int ic = this->ParaV->trace_loc_col[ iw2_all ];
    const long index = ir * this->ParaV->ncol + ic;
    
    if( index >= this->ParaV->nloc)
    {
        std::cout << " iw1_all = " << iw1_all << std::endl;
        std::cout << " iw2_all = " << iw2_all << std::endl;
        std::cout << " ir = " << ir << std::endl;
        std::cout << " ic = " << ic << std::endl;
        std::cout << " index = " << index << std::endl;
        std::cout << " this->ParaV->nloc = " << this->ParaV->nloc << std::endl;
        ModuleBase::WARNING_QUIT("LCAO_Matrix","set_force");
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
    const ModuleBase::Vector3<double> &dtau)
{
    // use iw1_all and iw2_all to set Hloc
    // becareful! The ir and ic may < 0!!!!!!!!!!!!!!!!
    const int ir = this->ParaV->trace_loc_row[ iw1_all ];
    const int ic = this->ParaV->trace_loc_col[ iw2_all ];
    const long index = ir * this->ParaV->ncol + ic;

    if( index >= this->ParaV->nloc)
    {
        std::cout << " iw1_all = " << iw1_all << std::endl;
        std::cout << " iw2_all = " << iw2_all << std::endl;
        std::cout << " ir = " << ir << std::endl;
        std::cout << " ic = " << ic << std::endl;
        std::cout << " index = " << index << std::endl;
        std::cout << " this->ParaV->nloc = " << this->ParaV->nloc << std::endl;
        ModuleBase::WARNING_QUIT("LCAO_Matrix","set_stress");
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
    if (mtype=='S') ModuleBase::GlobalFunc::ZEROS(this->Sloc.data(), this->Sloc.size());
    else if (mtype=='T') ModuleBase::GlobalFunc::ZEROS(this->Hloc_fixed.data(), this->Hloc_fixed.size());
    else if (mtype=='H') ModuleBase::GlobalFunc::ZEROS(this->Hloc.data(), this->Hloc.size());
    return;
}

void LCAO_Matrix::zeros_HSk(const char &mtype)
{
    if (mtype=='S') ModuleBase::GlobalFunc::ZEROS(this->Sloc2.data(), this->Sloc2.size());
    else if (mtype=='T') ModuleBase::GlobalFunc::ZEROS(this->Hloc_fixed2.data(), this->Hloc_fixed2.size());
    else if (mtype=='H') ModuleBase::GlobalFunc::ZEROS(this->Hloc2.data(), this->Hloc2.size());
    return;
}

void LCAO_Matrix::zeros_HSR(const char &mtype)
{
    if(GlobalV::NSPIN!=4)
    {
        if (mtype=='S') ModuleBase::GlobalFunc::ZEROS(this->SlocR.data(), this->SlocR.size());
        else if (mtype=='T') ModuleBase::GlobalFunc::ZEROS(this->Hloc_fixedR.data(), this->Hloc_fixedR.size());
    }
    else
    {
        if (mtype=='S') ModuleBase::GlobalFunc::ZEROS(this->SlocR_soc.data(), this->SlocR_soc.size());
        else if (mtype=='T') ModuleBase::GlobalFunc::ZEROS(this->Hloc_fixedR_soc.data(), this->Hloc_fixedR_soc.size());
    }
    return;
}

// Peize Lin add vtype='A' 2018-11-30
void LCAO_Matrix::print_HSk(const char &mtype, const char &vtype, const double &accuracy, std::ostream &os)
{
    ModuleBase::TITLE("LCAO_Matrix","print_HSk");
    if(mtype=='S') os << "Sloc2 matrix" << std::endl;
    else if(mtype=='T') os << "Hloc_fixed2 matrix" << std::endl;
    else if(mtype=='H') os << "Hloc2 matrix" << std::endl;
    else
    {
        ModuleBase::WARNING_QUIT("LCAO_Matrix::print_HSk","Check input parameter: mtype.");
    }

    if(vtype=='C') os << " Output norm."  << std::endl;
    else if(vtype=='R') os << " Output real part."  << std::endl;
    else if(vtype=='I') os << " Output imag part."  << std::endl;
    else if(vtype=='A') os << " Output std::complex." << std::endl;


    os << std::setprecision(8) << std::endl;
    for(int i=0; i<this->ParaV->nrow; i++)
    {
        os << " " ;
        for(int j=0; j<this->ParaV->ncol; j++)
        {
            const int index = i * this->ParaV->ncol + j;
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
    //				os << std::setw(15) << v;
                    os << v << "\t";
                }
                else
                {
    //				os << std::setw(15) << "0"; 
                    os << "0" << "\t"; 
                }
            }
        }
        os << std::endl;
    }
    os << std::endl;
    os << std::setprecision(6) << std::endl;
    return;
}


void LCAO_Matrix::print_HSgamma(const char &mtype, std::ostream &os)
{
    ModuleBase::TITLE("Parallel_Orbitals","print_HSgamma");

    GlobalV::ofs_running << " " << mtype << " matrix" << std::endl;
    GlobalV::ofs_running << " nrow=" << this->ParaV->nrow << std::endl;
    GlobalV::ofs_running << " ncol=" << this->ParaV->ncol << std::endl;
    GlobalV::ofs_running << " element number = " << this->ParaV->ncol << std::endl;

    if (mtype=='S')
    {
        os << std::setprecision(8);
        os << " print Sloc" << std::endl;
        for(int i=0; i<GlobalV::NLOCAL; ++i)
        {
            for(int j=0; j<GlobalV::NLOCAL; ++j)
            {
                double v = Sloc[i*this->ParaV->ncol+j];
                if( abs(v) > 1.0e-8)
                {
                    os << std::setw(15) << v;
                }
                else
                {
                    os << std::setw(15) << "0";
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
                double v = Hloc_fixed[i*this->ParaV->ncol+j];
                if( abs(v) > 1.0e-8)
                {
                    os << std::setw(15) << v;
                }
                else
                {
                    os << std::setw(15) << "0";
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
                double v = Hloc[i*this->ParaV->ncol+j];
                if( abs(v) > 1.0e-8)
                {
                    os << std::setw(15) << v;
                }
                else
                {
                    os << std::setw(15) << "0";
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
    for (long i=0; i<this->ParaV->nloc; i++)
    {
        Hloc[i] += Hloc_fixed[i];
    }
    return;
}

void LCAO_Matrix::update_Hloc2(const int &ik)
{
	for (long i = 0; i < this->ParaV->nloc; i++)
	{
		Hloc2[i] += Hloc_fixed2[i];
#ifdef __DEEPKS
		if(GlobalV::deepks_scf)
		{
			Hloc2[i] += GlobalC::ld.H_V_delta_k[ik][i];
		}
#endif
	}

	return;
}


void LCAO_Matrix::output_HSk(const char &mtype, std::string &fn)
{
    ModuleBase::TITLE("LCAO_Matrix","output_HSk");
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
    ModuleBase::TITLE("LCAO_Matrix","allocate_Hloc_fixedR_tr");

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
                    Hloc_fixedR_tr[ix][iy][iz] = new double[this->ParaV->nloc];
                    //HR_tr[ix][iy][iz] = new double[this->ParaV->nloc];
                    //SlocR_tr[ix][iy][iz] = new double[this->ParaV->nloc];
                    ModuleBase::GlobalFunc::ZEROS(Hloc_fixedR_tr[ix][iy][iz], this->ParaV->nloc);
                    //ModuleBase::GlobalFunc::ZEROS(HR_tr[ix][iy][iz], this->ParaV->nloc);
                    //ModuleBase::GlobalFunc::ZEROS(SlocR_tr[ix][iy][iz], this->ParaV->nloc);
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
                    Hloc_fixedR_tr_soc[ix][iy][iz] = new std::complex<double>[this->ParaV->nloc];
                    //HR_tr[ix][iy][iz] = new double[this->ParaV->nloc];
                    //SlocR_tr[ix][iy][iz] = new double[this->ParaV->nloc];
                    ModuleBase::GlobalFunc::ZEROS(Hloc_fixedR_tr_soc[ix][iy][iz], this->ParaV->nloc);
                    //ModuleBase::GlobalFunc::ZEROS(HR_tr[ix][iy][iz], this->ParaV->nloc);
                    //ModuleBase::GlobalFunc::ZEROS(SlocR_tr[ix][iy][iz], this->ParaV->nloc);
                }
            }
        }
    }
//std::cout<<"R_x: "<<R_x<<std::endl;
//std::cout<<"R_y: "<<R_y<<std::endl;
//std::cout<<"R_z: "<<R_z<<std::endl;
//std::cout<<"this->ParaV->nloc: "<<this->ParaV->nloc<<std::endl;
//std::cout<<"SlocR_tr 1-3-3-27: "<<SlocR_tr[1][3][3][27]<<std::endl;
//std::cout<<"Hloc_fixedR_tr 1-3-3-27: "<<Hloc_fixedR_tr[1][3][3][27]<<std::endl;

    return;
}

void LCAO_Matrix::allocate_HR_tr(void)
{
    ModuleBase::TITLE("LCAO_Matrix","allocate_HR_tr");

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
                    HR_tr[ix][iy][iz] = new double[this->ParaV->nloc];
                    ModuleBase::GlobalFunc::ZEROS(HR_tr[ix][iy][iz], this->ParaV->nloc);
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
                    HR_tr_soc[ix][iy][iz] = new std::complex<double>[this->ParaV->nloc];
                    ModuleBase::GlobalFunc::ZEROS(HR_tr_soc[ix][iy][iz], this->ParaV->nloc);
                }
            }
        }
    }

    return;
}

void LCAO_Matrix::allocate_SlocR_tr(void)
{
    ModuleBase::TITLE("LCAO_Matrix","allocate_SlocR_tr");

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
                    SlocR_tr[ix][iy][iz] = new double[this->ParaV->nloc];
                    ModuleBase::GlobalFunc::ZEROS(SlocR_tr[ix][iy][iz], this->ParaV->nloc);
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
                    SlocR_tr_soc[ix][iy][iz] = new std::complex<double>[this->ParaV->nloc];
                    ModuleBase::GlobalFunc::ZEROS(SlocR_tr_soc[ix][iy][iz], this->ParaV->nloc);
                }
            }
        }
    }

    return;
}

void LCAO_Matrix::destroy_Hloc_fixedR_tr(void)
{
    ModuleBase::TITLE("LCAO_Matrix","destroy_Hloc_fixed2_R");

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
    const int ir = this->ParaV->trace_loc_row[ iw1_all ];
    const int ic = this->ParaV->trace_loc_col[ iw2_all ];

//std::cout<<"ir: "<<ir<<std::endl;
//std::cout<<"ic: "<<ic<<std::endl;
    long index;
    if(GlobalV::KS_SOLVER=="genelpa" || GlobalV::KS_SOLVER=="scalapack_gvx")
    {
        index=ic*this->ParaV->nrow+ir;
//std::cout<<"index: "<<index<<std::endl;
    }
    else
    {
        index=ir*this->ParaV->ncol+ic;
//std::cout<<"index: "<<index<<std::endl;
    }

//std::cout<<"this->ParaV->nloc: "<<this->ParaV->nloc<<std::endl;
    assert(index < this->ParaV->nloc);
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
    const int ir = this->ParaV->trace_loc_row[ iw1_all ];
    const int ic = this->ParaV->trace_loc_col[ iw2_all ];

//std::cout<<"ir: "<<ir<<std::endl;
//std::cout<<"ic: "<<ic<<std::endl;
    long index;
    if(GlobalV::KS_SOLVER=="genelpa" || GlobalV::KS_SOLVER=="scalapack_gvx")
    {
        index=ic*this->ParaV->nrow+ir;
//std::cout<<"index: "<<index<<std::endl;
    }
    else
    {
        index=ir*this->ParaV->ncol+ic;
//std::cout<<"index: "<<index<<std::endl;
    }

//std::cout<<"this->ParaV->nloc: "<<this->ParaV->nloc<<std::endl;
    assert(index < this->ParaV->nloc);
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

void LCAO_Matrix::destroy_HS_R_sparse(void)
{
    ModuleBase::TITLE("LCAO_Matrix","destroy_HS_R_sparse");

    if (GlobalV::NSPIN != 4)
    {
        std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> empty_HR_sparse_up;
        std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> empty_HR_sparse_down;
        std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, double>>> empty_SR_sparse;
        HR_sparse[0].swap(empty_HR_sparse_up);
        HR_sparse[1].swap(empty_HR_sparse_down);
        SR_sparse.swap(empty_SR_sparse);
    }
    else
    {
        std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>> empty_HR_soc_sparse;
        std::map<Abfs::Vector3_Order<int>, std::map<size_t, std::map<size_t, std::complex<double>>>> empty_SR_soc_sparse;
        HR_soc_sparse.swap(empty_HR_soc_sparse);
        SR_soc_sparse.swap(empty_SR_soc_sparse);
    }

    // 'all_R_coor' has a small memory requirement and does not need to be deleted.
    // std::set<Abfs::Vector3_Order<int>> empty_all_R_coor;
    // all_R_coor.swap(empty_all_R_coor);

    return;
}
