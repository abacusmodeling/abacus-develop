//caoyu add 2021-03-2
#ifdef __DEEPKS

#include "LCAO_descriptor.h"
#include "LCAO_matrix.h"
#include "../module_base/lapack_connector.h"
#include "../module_base/intarray.h"
#include "../module_base/complexmatrix.h"
#include "global_fp.h"
#include "../src_pw/global.h"
#include "../src_io/winput.h"

#include <torch/script.h>
#include <torch/csrc/autograd/autograd.h>
#include <npy.hpp>

LCAO_Descriptor ld;

LCAO_Descriptor::LCAO_Descriptor()
{
    alpha_index = new IntArray[1];
    inl_index = new IntArray[1];
    inl_l = new int[1];
    d = new double[1];
    H_V_delta = new double[1]; 
}

LCAO_Descriptor::~LCAO_Descriptor()
{
    delete[] alpha_index;
    delete[] inl_index;
    delete[] inl_l;
    delete[] d;
    delete[] H_V_delta;

    //delete S_mu_alpha**
    for (int inl = 0;inl < this->inlmax;inl++)
    {
        delete[] S_mu_alpha[inl];
    }
    delete[] S_mu_alpha;

    //delete DS_mu_alpha**
    for (int inl = 0;inl < this->inlmax;inl++)
    {
        delete[] DS_mu_alpha_x[inl];
        delete[] DS_mu_alpha_y[inl];
        delete[] DS_mu_alpha_z[inl];
    }
    delete[] DS_mu_alpha_x;
    delete[] DS_mu_alpha_y;
    delete[] DS_mu_alpha_z;

    //delete pdm**
    for (int inl = 0;inl < this->inlmax;inl++)
    {
        delete[] pdm[inl];
    }
    delete[] pdm;
    //delete gedm**
    for (int inl = 0;inl < this->inlmax;inl++)
    {
        delete[] gedm[inl];
    }
    delete[] gedm;
}

void LCAO_Descriptor::init(int lm, int nm, int tot_inl)
{
    TITLE("LCAO_Descriptor", "init");
    GlobalV::ofs_running << " Initialize the descriptor index for deepks (lcao line)" << endl;

    assert(lm >= 0);
    assert(nm >= 0);
    assert(tot_inl >= 0);
    
    this->lmaxd = lm;
    this->nmaxd = nm;
    this->inlmax = tot_inl;
    GlobalV::ofs_running << " lmax of descriptor = " << this->lmaxd << endl;
    GlobalV::ofs_running << " nmax of descriptor= " << nmaxd << endl;

    //init S_mu_alpha**
    this->S_mu_alpha = new double* [this->inlmax];    //inl
    for (int inl = 0;inl < this->inlmax;inl++)
    {
        this->S_mu_alpha[inl] = new double[GlobalV::NLOCAL * (2 * this->lmaxd + 1)];     //GlobalV::NLOCAL*nm
        ZEROS(S_mu_alpha[inl], GlobalV::NLOCAL * (2 * this->lmaxd+ 1));
    }

    //init F_delta
    F_delta.create(GlobalC::ucell.nat, 3);
    //init DS_mu_alpha**
    this->DS_mu_alpha_x = new double* [this->inlmax];
    this->DS_mu_alpha_y = new double* [this->inlmax];
    this->DS_mu_alpha_z = new double* [this->inlmax];
    for (int inl = 0;inl < this->inlmax;inl++)
    {
        this->DS_mu_alpha_x[inl] = new double[GlobalV::NLOCAL * (2 * this->lmaxd + 1)];
        this->DS_mu_alpha_y[inl] = new double[GlobalV::NLOCAL * (2 * this->lmaxd + 1)];
        this->DS_mu_alpha_z[inl] = new double[GlobalV::NLOCAL * (2 * this->lmaxd + 1)];
        ZEROS(DS_mu_alpha_x[inl], GlobalV::NLOCAL * (2 * this->lmaxd + 1));
        ZEROS(DS_mu_alpha_y[inl], GlobalV::NLOCAL * (2 * this->lmaxd + 1));
        ZEROS(DS_mu_alpha_z[inl], GlobalV::NLOCAL * (2 * this->lmaxd + 1));
    }
    //init pdm**
    const int PDM_size = (this->lmaxd * 2 + 1) * (this->lmaxd * 2 + 1);
    this->pdm = new double* [this->inlmax];
    for (int inl = 0;inl < this->inlmax;inl++)
    {
        this->pdm[inl] = new double[PDM_size];
        ZEROS(this->pdm[inl], PDM_size);
    }
    //init gedm**
    this->gedm = new double* [this->inlmax];
    for (int inl = 0;inl < this->inlmax;inl++)
    {
        this->gedm[inl] = new double[PDM_size];
        ZEROS(this->gedm[inl], PDM_size);
    }
    // cal n(descriptor) per atom , related to Lmax, nchi(L) and m. (not total_nchi!)
	this->des_per_atom=0; // mohan add 2021-04-21
    for (int l = 0; l <= this->lmaxd; l++)
    {
        this->des_per_atom += GlobalC::ORB.Alpha[0].getNchi(l) * (2 * l + 1);
    }
    this->n_descriptor = GlobalC::ucell.nat * this->des_per_atom;

    this->init_index();
    
    return;
}


void LCAO_Descriptor::init_index(void)
{
    
    delete[] this->alpha_index;
    this->alpha_index = new IntArray[GlobalC::ucell.ntype];
    delete[] this->inl_index;
    this->inl_index = new IntArray[GlobalC::ucell.ntype];
    delete[] this->inl_l;
    this->inl_l = new int[this->inlmax];
    ZEROS(this->inl_l, this->inlmax);

    int inl = 0;
    int alpha = 0;
    for (int it = 0; it < GlobalC::ucell.ntype; it++)
    {
        this->alpha_index[it].create(
            GlobalC::ucell.atoms[it].na,
            this->lmaxd + 1, // l starts from 0
            this->nmaxd,
            2 * this->lmaxd + 1); // m ==> 2*l+1

        this->inl_index[it].create(
            GlobalC::ucell.atoms[it].na,
            this->lmaxd + 1,
            this->nmaxd); 

        GlobalV::ofs_running << "Type " << it + 1
                    << " number_of_atoms " << GlobalC::ucell.atoms[it].na << endl;

        for (int ia = 0; ia < GlobalC::ucell.atoms[it].na; ia++)
        {
            //alpha
            for (int l = 0; l < this->lmaxd + 1; l++)
            {
                for (int n = 0; n < GlobalC::ORB.Alpha[0].getNchi(l); n++)
                {
                    for (int m = 0; m < 2 * l + 1; m++)
                    {
                        this->alpha_index[it](ia, l, n, m) = alpha;
                        alpha++;
                    }
                    this->inl_index[it](ia, l, n) = inl;
                    this->inl_l[inl] = l;
                    inl++;
                }
            }
        }//end ia
    }//end it
    assert(this->n_descriptor == alpha);
    assert(GlobalC::ucell.nat * GlobalC::ORB.Alpha[0].getTotal_nchi() == inl);
    GlobalV::ofs_running << "descriptors_per_atom " << this->des_per_atom << endl;
    GlobalV::ofs_running << "total_descriptors " << this->n_descriptor << endl;
	
	return;
}


void LCAO_Descriptor::build_S_descriptor(const bool& calc_deri)
{
    TITLE("LCAO_Descriptor", "build_S_descriptor");
    //array to store data
    double olm[3] = {0.0, 0.0, 0.0};

    //\sum{T} e**{ikT} <\phi_{ia}|d\phi_{k\beta}(T)>	//???
    Vector3<double> tau1, tau2, dtau;
    Vector3<double> dtau1, dtau2, tau0;
    for (int T1 = 0; T1 < GlobalC::ucell.ntype; ++T1)
    {
        Atom *atom1 = &GlobalC::ucell.atoms[T1];
        for (int I1 = 0; I1 < atom1->na; ++I1)
        {
            tau1 = atom1->tau[I1];
            //GlobalC::GridD.Find_atom(tau1);
            GlobalC::GridD.Find_atom(GlobalC::ucell, tau1, T1, I1);

            for (int ad = 0; ad < GlobalC::GridD.getAdjacentNum() + 1; ++ad)
            {
                const int T2 = GlobalC::GridD.getType(ad);
                const int I2 = GlobalC::GridD.getNatom(ad);
                //Atom *atom2 = &GlobalC::ucell.atoms[T2];
                tau2 = GlobalC::GridD.getAdjacentTau(ad);
                dtau = tau2 - tau1;
                double distance = dtau.norm() * GlobalC::ucell.lat0;
                double rcut = GlobalC::ORB.Phi[T1].getRcut() + GlobalC::ORB.Alpha[0].getRcut(); //Rcut is subject to GlobalC::ORB.Phi to keep dimension of S_mu_alpha same as Sloc
                if (distance < rcut)
                {
                    int iw1_all = GlobalC::ucell.itiaiw2iwt(T1, I1, 0); //iw1_all = combined index (it, ia, iw)

                    for (int jj = 0; jj < atom1->nw * GlobalV::NPOL; ++jj)
                    {
                        const int jj0 = jj / GlobalV::NPOL;
                        const int L1 = atom1->iw2l[jj0];
                        const int N1 = atom1->iw2n[jj0];
                        const int m1 = atom1->iw2m[jj0];

                        for (int L2 = 0; L2 <= GlobalC::ORB.Alpha[0].getLmax(); ++L2)
                        {
                            for (int N2 = 0; N2 < GlobalC::ORB.Alpha[0].getNchi(L2); ++N2)
                            {
                                for (int m2 = 0; m2 < 2 * L2 + 1; ++m2)
                                {
                                    olm[0] = olm[1] = olm[2] = 0.0;
                                    if (!calc_deri)
                                    {
                                        GlobalC::UOT.snap_psipsi(olm, 0, 'D', tau1,
                                                T1, L1, m1, N1, GlobalC::GridD.getAdjacentTau(ad),
                                                T2, L2, m2, N2, GlobalV::NSPIN);
                                        if (GlobalV::GAMMA_ONLY_LOCAL)
                                        {
                                            this->set_S_mu_alpha(
												iw1_all,  // index of numerical atomic orbital mu
												inl_index[T2](I2,L2,N2), // index of alpha (for descriptor)
												m2, 
												olm[0]);// overlap value
                                        }
                                    }
                                    else
                                    {
                                        GlobalC::UOT.snap_psipsi(olm, 1, 'D', tau1,
                                            T1, L1, m1, N1, GlobalC::GridD.getAdjacentTau(ad),
                                            T2, L2, m2, N2, GlobalV::NSPIN);
                                        if (GlobalV::GAMMA_ONLY_LOCAL)
                                        {
                                            this->set_DS_mu_alpha(
												iw1_all, 
												inl_index[T2](I2,L2,N2), 
												m2, 
												olm[0], // derivative for x direction 
												olm[1], // derivative for y direction
												olm[2]);// derivative for z direction
                                        }
                                    }

                                } //m2
                            }     //N2
                        }         //nw2(L2)
                        ++iw1_all;
                    } // nw1
                }     // distance
            }         // ad
        } // I1
    }     // T1
    if (!GlobalV::GAMMA_ONLY_LOCAL)
    {
        WARNING_QUIT("LCAO_Descriptor::build_S_descriptor", 
		"muti-kpoint method for descriptor is not implemented yet! ");
    }
    return;
}


void LCAO_Descriptor::set_S_mu_alpha(const int &iw1_all, const int &inl, const int &im, const double &v)
{
    //const int ir = GlobalC::ParaO.trace_loc_row[iw1_all];
    //const int ic = GlobalC::ParaO.trace_loc_col[iw2_all];
    //no parellel yet
    const int ir = iw1_all;
    const int ic = im;
    //const int index = ir * GlobalC::ParaO.ncol + ic;
    int index=0; // mohan add 2021-07-26
    if (GlobalV::KS_SOLVER == "genelpa" || GlobalV::KS_SOLVER == "scalapack_gvx") // save the matrix as column major format
    {
        index = ic * GlobalV::NLOCAL + ir;
    }
    else
    {
        index = ir * (2*inl_l[inl]+1)  + ic; //row: lcao orbitals; col: descriptor basis
    }
    this->S_mu_alpha[inl][index] += v;
    return;
}


void LCAO_Descriptor::cal_projected_DM(void)
{
    TITLE("LCAO_Descriptor", "cal_projected_DM");
    //step 1: get dm: the coefficient of wfc, not charge density
    double *dm = new double[GlobalV::NLOCAL * GlobalV::NLOCAL];
    ZEROS(dm, GlobalV::NLOCAL * GlobalV::NLOCAL);
    this->getdm(dm);

    //step 2: get S_alpha_mu and S_nu_beta
    double **ss = this->S_mu_alpha;

    //step 3 : multiply: cal ST*DM*S
    
    //init tmp_pdm*
    const int tmp_PDM_size = GlobalV::NLOCAL * (lmaxd*2+1);
    double* tmp_pdm = new double[tmp_PDM_size];
    ZEROS(tmp_pdm, tmp_PDM_size);

    for (int inl = 0;inl < inlmax;inl++)
    {   
        int nm = 2 * inl_l[inl] + 1;   //1,3,5,...
        const char t = 'T';  //transpose
        const char nt = 'N'; //non transpose
        const double alpha = 1;
        const double beta = 0;
        double *a = dm;
        double *b = ss[inl];
        double *c = tmp_pdm;

        dgemm_(&nt, 
			&nt, 
			&GlobalV::NLOCAL, 
			&nm, 
			&GlobalV::NLOCAL, 
			&alpha, 
			a, 
			&GlobalV::NLOCAL, 
			b, 
			&GlobalV::NLOCAL, 
			&beta, 
			c, 
			&GlobalV::NLOCAL); //DM*S

        a = ss[inl];
        b = c;
        c = this->pdm[inl];

        dgemm_(
			&t, 
			&nt, 
			&nm, 
			&nm, 
			&GlobalV::NLOCAL, 
			&alpha, 
			a, 
			&GlobalV::NLOCAL, 
			b, 
			&GlobalV::NLOCAL, 
			&beta, 
			c, 
			&nm); //ST*DM*S
    }
    
    delete[] tmp_pdm;
    delete[] dm;
    return;
}

void LCAO_Descriptor::cal_descriptor(void)
{
    TITLE("LCAO_Descriptor", "cal_descriptor");
    delete[] d;
    d = new double[this->n_descriptor];

    //==========print preparation=============
    GlobalV::ofs_running << " print out each DM_inl" << endl;
    ofstream ofs;
    stringstream ss;
    ss << winput::spillage_outdir << "/"
       << "projected_DM.dat";
    if (GlobalV::MY_RANK == 0)
    {
        ofs.open(ss.str().c_str());
    }
    //==========print preparation=============
    const int lmax = GlobalC::ORB.get_lmax_d();
    int id = 0;
    for (int it = 0; it < GlobalC::ucell.ntype; it++)
    {
        for (int ia = 0; ia < GlobalC::ucell.atoms[it].na; ia++)
        {
            ofs << GlobalC::ucell.atoms[it].label << " atom_index " 
			<< ia + 1 << " n_descriptor " 
			<< this->des_per_atom << endl;

            for (int l = 0; l <= lmax; l++)
            {
                int nmax = GlobalC::ORB.Alpha[0].getNchi(l);
                for (int n = 0; n < nmax; n++)
                {
                    const int dim = 2 * l + 1;
                    const int inl = inl_index[it](ia, l, n);
                    // descriptor for atom (it, ia)
                    ComplexMatrix des(dim, dim);
                    for (int m = 0; m < dim; m++)
                    {
                        for (int m2 = 0; m2 < dim; m2++)
                        {
                            int index = m * dim + m2;
                            complex<double> tmp(this->pdm[inl][index], 0);
                            des(m, m2) += tmp;
                        }
                    }

                    this->print_projected_DM(ofs, des, it, ia, l, n);

                    //GlobalV::ofs_running << "dimension of des is " << 2 * l + 1 << endl;
                    if (l == 0)
                    {
                        this->d[id] = des(0, 0).real();
                        ++id;
                    }
                    else
                    {
                        // diagonalizae
                        // assume des matrix is Hermitian
                        char jobz = 'N'; // eigenvalues only
                        char uplo = 'U'; // upper matrix is stored
                        int ndim = des.nr;
                        double *tmpd = new double[ndim]();
                        const int lwork = 2 * ndim;
                        complex<double> *work = new complex<double>[lwork]();
                        double *rwork = new double[3 * ndim - 2]();
                        int infor = 0;
                        // diag by calling zheev
                        LapackConnector::zheev(jobz, uplo, ndim, des, ndim, tmpd, work, lwork, rwork, &infor);
                        // put the eigenvalues into d (descriptor)
                        for (int idim = 0; idim < ndim; ++idim)
                        {
                            this->d[id] = tmpd[idim];
                            ++id;
                        }
                        delete[] tmpd;
                        delete[] rwork;
                        delete[] work;
                    }
                } //n
            }     //l
        }         //ia
    }             //it
    this->cal_descriptor_tensor();  //use torch::symeig
    this->print_descriptor();
    return;
}



void LCAO_Descriptor::print_projected_DM(
	ofstream& ofs, 
	ComplexMatrix& des, 
	const int& it, 
	const int& ia, 
	const int& l, 
	const int& n)
{
    ofs << "L=" << l << "   N=" << n << endl;
    for (int i = 0; i < 2 * l + 1; i++)
    {
        for (int j = 0; j < 2 * l + 1; j++)
        {
            ofs << des(i, j).real() << " ";
        }
        ofs << endl;
    }
    return;
}

void LCAO_Descriptor::print_descriptor(void)
{
    TITLE("LCAO_Descriptor", "print_descriptor");
    ofstream ofs;
    stringstream ss;
    // the parameter 'winput::spillage_outdir' is read from INPUTw.
    ss << winput::spillage_outdir << "/"
       << "descriptor.dat";
    if (GlobalV::MY_RANK == 0)
    {
        ofs.open(ss.str().c_str());
    }
    for (int it = 0; it < GlobalC::ucell.ntype; it++)
    {
        for (int ia = 0; ia < GlobalC::ucell.atoms[it].na; ia++)
        {
            ofs << GlobalC::ucell.atoms[it].label << " atom_index " 
			<< ia + 1 << " n_descriptor " << this->des_per_atom << endl;

            int id0 = this->alpha_index[it](ia, 0, 0, 0);

            for (int id = id0; id < id0 + this->des_per_atom; ++id)
            {
                if ((id - id0) > 0 && (id - id0) % 8 == 0)
				{
                    ofs << endl;
				}
                ofs << d[id] << " ";
            }
            ofs << endl << endl;
        }
    }
    GlobalV::ofs_running << "descriptors are printed" << endl;

    return;
}

void LCAO_Descriptor::set_DS_mu_alpha(
	const int& iw1_all, 
	const int& inl, 
	const int& im,
    const double& vx, 
	const double& vy, 
	const double& vz)
{
    //const int ir = GlobalC::ParaO.trace_loc_row[iw1_all];
    //const int ic = GlobalC::ParaO.trace_loc_col[iw2_all];
    //no parellel yet
    const int ir = iw1_all;
    const int ic = im;
    //const int index = ir * GlobalC::ParaO.ncol + ic;
    int index;
    if (GlobalV::KS_SOLVER == "genelpa" || GlobalV::KS_SOLVER == "scalapack_gvx") // save the matrix as column major format
    {
        index = ic * GlobalV::NLOCAL + ir;
    }
    else
    {
        index = ir * (2*inl_l[inl]+1)  + ic; //row: lcao orbitals; col: descriptor basis
    }
    this->DS_mu_alpha_x[inl][index] += vx;
    this->DS_mu_alpha_y[inl][index] += vy;
    this->DS_mu_alpha_z[inl][index] += vz;
    return;
}

void LCAO_Descriptor::getdm(double* dm)
{
    for (int i = 0; i < GlobalC::LOC.wfc_dm_2d.dm_gamma[0].nr; i++)
    {
        for (int j = 0; j < GlobalC::LOC.wfc_dm_2d.dm_gamma[0].nc; j++)
        {
            dm[i * GlobalV::NLOCAL + j] = GlobalC::LOC.wfc_dm_2d.dm_gamma[0](i, j); //only consider default GlobalV::NSPIN = 1
        }
    }
}

void LCAO_Descriptor::cal_gdmx(matrix &dm)
{
    TITLE("LCAO_Descriptor", "cal_gdmx");
    //get DS_alpha_mu and S_nu_beta
    double** ss = this->S_mu_alpha;
    double** dsx = this->DS_mu_alpha_x;
    double** dsy = this->DS_mu_alpha_y;
    double** dsz = this->DS_mu_alpha_z;
    for (int inl = 0;inl < inlmax;inl++)
    {
        //dE/dD will be multiplied in cal_f_delta, here only calculate dD/dx_I
        int nm = 2 * inl_l[inl] + 1;   //1,3,5,...
        for (int i =0; i < GlobalV::NLOCAL;++i) 
        {
            const int iat = GlobalC::ucell.iwt2iat[i];//the atom whose force being calculated
            for (int j= 0;j < GlobalV::NLOCAL; ++j)
            {
                const int mu = GlobalC::ParaO.trace_loc_row[j];
                const int nu = GlobalC::ParaO.trace_loc_col[i];
                if (mu >= 0 && nu >= 0)
                {
                    for (int m1 = 0;m1 < nm;++m1)
                    {
                        for (int m2 = 0;m2 < nm;++m2)
                        {
                            for (int is = 0;is < GlobalV::NSPIN;++is)
                            {
                                if (GlobalV::KS_SOLVER == "genelpa" || GlobalV::KS_SOLVER == "scalapack_gvx") // save the matrix as column major format
                                {
                                    gdmx[iat][inl][m1 * nm + m2] += 4 * dsx[inl][m1 * GlobalV::NLOCAL + mu] * dm(mu, nu) * ss[inl][m2 * GlobalV::NLOCAL + nu];
                                    gdmy[iat][inl][m1 * nm + m2] += 4 * dsy[inl][m1 * GlobalV::NLOCAL + mu] *  dm(mu, nu)  * ss[inl][m2 * GlobalV::NLOCAL + nu];
                                    gdmz[iat][inl][m1 * nm + m2] += 4 * dsz[inl][m1 * GlobalV::NLOCAL + mu] * dm(mu, nu)  * ss[inl][m2 * GlobalV::NLOCAL + nu];
                                }
                                else
                                {
                                    gdmx[iat][inl][m1 * nm + m2] += 4 * dsx[inl][mu * nm + m1] * dm(mu, nu) * ss[inl][nu * nm + m2];
                                    gdmy[iat][inl][m1 * nm + m2] += 4 * dsy[inl][mu * nm + m1] * dm(mu, nu) * ss[inl][nu * nm + m2];
                                    gdmz[iat][inl][m1 * nm + m2] += 4 * dsz[inl][mu * nm + m1] *  dm(mu, nu) * ss[inl][nu * nm + m2];
                                }
                            }
                        }//end m2
                    } //end m1
                }//end if
            }//end nu
        }//end mu
    }//end inl
    return;
}

void LCAO_Descriptor::init_gdmx()
{
    this->gdmx = new double** [GlobalC::ucell.nat];
    this->gdmy = new double** [GlobalC::ucell.nat];
    this->gdmz = new double** [GlobalC::ucell.nat];
    for (int iat = 0;iat < GlobalC::ucell.nat;iat++)
    {
        this->gdmx[iat] = new double* [inlmax];
        this->gdmy[iat] = new double* [inlmax];
        this->gdmz[iat] = new double* [inlmax];
        for (int inl = 0;inl < inlmax;inl++)
        {
            this->gdmx[iat][inl] = new double [(2 * lmaxd + 1) * (2 * lmaxd + 1)];
            this->gdmy[iat][inl] = new double [(2 * lmaxd + 1) * (2 * lmaxd + 1)];
            this->gdmz[iat][inl] = new double[(2 * lmaxd + 1) * (2 * lmaxd + 1)];
            ZEROS(gdmx[iat][inl], (2 * lmaxd + 1) * (2 * lmaxd + 1));
            ZEROS(gdmy[iat][inl], (2 * lmaxd + 1) * (2 * lmaxd + 1));
            ZEROS(gdmz[iat][inl], (2 * lmaxd + 1) * (2 * lmaxd + 1));
        }
    }
    return;
}

void LCAO_Descriptor::del_gdmx(void)
{
    for (int iat = 0;iat < GlobalC::ucell.nat;iat++)
    {
        for (int inl = 0;inl < inlmax;inl++)
        {
            delete[] this->gdmx[iat][inl];
            delete[] this->gdmy[iat][inl];
            delete[] this->gdmz[iat][inl];
        }
        delete[] this->gdmx[iat];
        delete[] this->gdmy[iat];
        delete[] this->gdmz[iat];
    }
    delete[] this->gdmx;
    delete[] this->gdmy;
    delete[] this->gdmz;
    return;
}

void LCAO_Descriptor::cal_v_delta(const string& model_file)
{
    TITLE("LCAO_Descriptor", "cal_v_delta");
    //1.  (dE/dD)<alpha_m'|psi_nv>
    this->load_model(model_file);
    this->cal_gedm();

    //2. multiply and sum
    double* tmp_v1 = new double[(2 * lmaxd + 1) * GlobalV::NLOCAL];
    ZEROS(tmp_v1, (2 * lmaxd + 1) * GlobalV::NLOCAL);
    double* tmp_v2 = new double[GlobalV::NLOCAL *GlobalV::NLOCAL];
    ZEROS(tmp_v2, GlobalV::NLOCAL * GlobalV::NLOCAL);
    //init H_V_delta
    delete[] this->H_V_delta;
    this->H_V_delta = new double[GlobalV::NLOCAL * GlobalV::NLOCAL];
    
    for (int inl = 0;inl < inlmax;inl++)
    {
        int nm = 2 * inl_l[inl] + 1;   //1,3,5,...
        const char t = 'T';  //transpose
        const char nt = 'N'; //non transpose
        const double alpha = 1;
        const double beta = 0;
        double* a = this->gedm[inl];//[nm][nm]
        double* b = S_mu_alpha[inl];//[GlobalV::NLOCAL][nm]--trans->[nm][GlobalV::NLOCAL]
        double* c = tmp_v1;
        
        dgemm_(&nt, &t, &nm, &GlobalV::NLOCAL, &nm, &alpha, a, &nm, b, &GlobalV::NLOCAL, &beta, c, &nm);

        //2. <psi_mu|alpha_m>*(dE/dD)*<alpha_m'|psi_nv>
        a = b; //[GlobalV::NLOCAL][nm]
        b = c;//[nm][GlobalV::NLOCAL]
        c = tmp_v2;//[GlobalV::NLOCAL][GlobalV::NLOCAL]
        dgemm_(&nt, &nt, &GlobalV::NLOCAL, &GlobalV::NLOCAL, &nm, &alpha, a, &GlobalV::NLOCAL, b, &nm, &beta, c, &GlobalV::NLOCAL);

        //3. sum of Inl
        for (int i = 0;i < GlobalV::NLOCAL * GlobalV::NLOCAL;++i)
        {
            this->H_V_delta[i] += c[i];
        }
    }
    delete[] tmp_v1;
    delete[] tmp_v2;
    GlobalV::ofs_running << "finish calculating H_V_delta" << endl;
    return;
}

void LCAO_Descriptor::cal_f_delta(matrix& dm)
{
    TITLE("LCAO_Descriptor", "cal_f_delta");
    int iat = 0;    //check if the index same as GlobalC::ucell.iw2iat or not !!
    for (int it = 0;it < GlobalC::ucell.ntype;++it)
    {
        for (int ia = 0;ia < GlobalC::ucell.atoms[it].na;++ia)
        {
            for (int inl = 0;inl < inlmax;++inl)
            {
                int nm = 2 * inl_l[inl] + 1;
                
                //1. cal gdmx
                this->init_gdmx();
                this->cal_gdmx(dm);

                //2.multiply and sum for each atom
                // \sum_{Inl}\sum_{mm'} <gedm, gdmx>_{mm'}
                //notice: sum of multiplied corresponding element(mm') , not matrix multiplication !
                for (int m1 = 0;m1 < nm;++m1)
                {
                    for (int m2 = 0; m2 < nm;++m2)
                    {
                        this->F_delta(iat, 0) += this->gedm[inl][m1 * nm + m2] * gdmx[iat][inl][m1 * nm + m2];
                        this->F_delta(iat, 1) += this->gedm[inl][m1 * nm + m2] * gdmy[iat][inl][m1 * nm + m2];
                        this->F_delta(iat, 2) += this->gedm[inl][m1 * nm + m2] * gdmz[iat][inl][m1 * nm + m2];
                    }
                }
                
            }//end inl
            ++iat;
        }//end ia
    }//end it
    this->del_gdmx();
    return;
}

void LCAO_Descriptor::cal_descriptor_tensor()
{
    TITLE("LCAO_Descriptor", "cal_descriptor_tensor");
    //init pdm_tensor and d_tensor
    torch::Tensor tmp;
    for (int inl = 0;inl < this->inlmax;++inl)
    {
        int nm = 2 * inl_l[inl] + 1;
        tmp = torch::ones({ nm, nm }, torch::TensorOptions().dtype(torch::kFloat64));
        for (int m1 = 0;m1 < nm;++m1)
        {
            for (int m2 = 0;m2 < nm;++m2)
            {
                tmp.index_put_({m1, m2}, this->pdm[inl][m1 * nm + m2]);
            }
        }
        //torch::Tensor tmp = torch::from_blob(this->pdm[inl], { nm, nm }, torch::requires_grad());
        tmp.requires_grad_(true);
        this->pdm_tensor.push_back(tmp);
        this->d_tensor.push_back(torch::ones({ nm }, torch::requires_grad(true)));
    }

    //cal d_tensor
    for (int inl = 0;inl < inlmax;++inl)
    {
        torch::Tensor vd;
        std::tuple<torch::Tensor, torch::Tensor> d_v(this->d_tensor[inl], vd);
        d_v = torch::symeig(pdm_tensor[inl], /*eigenvalues=*/true, /*upper=*/true);
        d_tensor[inl] = std::get<0>(d_v);
    }
    return;
}

void LCAO_Descriptor::load_model(const string& model_file)
{
    try {
        module = torch::jit::load(model_file);
    }
    catch (const c10::Error& e) {
        std::cerr << "error loading the model" << std::endl;
        return;
    }
}
void LCAO_Descriptor::cal_gedm()
{
    TITLE("LCAO_Descriptor", "cal_gedm");
    //forward
    std::vector<torch::jit::IValue> inputs;
    //input_dim:(natom, des_per_atom)
    inputs.push_back(torch::cat(d_tensor, /*dim=*/0).reshape({ GlobalC::ucell.nat, des_per_atom }));
    std::vector<torch::Tensor> ec;
    ec.push_back(module.forward(inputs).toTensor());    //Hartree
    this->E_delta = ec[0].item().toDouble() * 2;//Ry; *2 is for Hartree to Ry
    
    //cal gedm
    std::vector<torch::Tensor> grad_shell;
    grad_shell.push_back(torch::ones_like(ec[0]));
    this->gedm_tensor = torch::autograd::grad(ec, this->pdm_tensor, grad_shell);

    //gedm_tensor(Hartree) to gedm(Ry)
    for (int inl = 0;inl < inlmax;++inl)
    {
        int nm = 2 * inl_l[inl] + 1;
        for (int m1 = 0;m1 < nm;++m1)
        {
            for (int m2 = 0;m2 < nm;++m2)
            {
                int index = m1 * nm + m2;
                //*2 is for Hartree to Ry
                this->gedm[inl][index] = this->gedm_tensor[inl].index({ m1,m2 }).item().toDouble() * 2;
            }
        }
    }
    return;
}

void LCAO_Descriptor::print_H_V_delta()
{
    TITLE("LCAO_Descriptor", "print_H_V_delta");
    ofstream ofs;
    stringstream ss;
    // the parameter 'winput::spillage_outdir' is read from INPUTw.
    ss << winput::spillage_outdir << "/"
       << "H_V_delta.dat";
    if (GlobalV::MY_RANK == 0)
    {
        ofs.open(ss.str().c_str());
    }
    ofs << "E_delta(Ry) from deepks model: " << this->E_delta << endl;
    ofs << "E_delta(eV) from deepks model: " << this->E_delta * Hartree_to_eV << endl;
    ofs << "H_delta(Hartree)(gamma only)) from deepks model: " << endl;
    for (int i = 0;i < GlobalV::NLOCAL;++i)
    {
        for (int j = 0;j < GlobalV::NLOCAL;++j)
        {
            ofs<< setw(12)<< this->H_V_delta[i * GlobalV::NLOCAL + j] << " ";
        }
        ofs << endl;
    }
    ofs << "H_delta(eV)(gamma only)) from deepks model: " << endl;
    for (int i = 0;i < GlobalV::NLOCAL;++i)
    {
        for (int j = 0;j < GlobalV::NLOCAL;++j)
        {
            ofs<< setw(12)<< this->H_V_delta[i * GlobalV::NLOCAL + j] *Hartree_to_eV<< " ";
        }
        ofs << endl;
    }
    GlobalV::ofs_running << "H_delta is printed" << endl;
    return;
}

void LCAO_Descriptor::print_F_delta()
{
    TITLE("LCAO_Descriptor", "print_F_delta");
    ofstream ofs;
    stringstream ss;
    // the parameter 'winput::spillage_outdir' is read from INPUTw.
    ss << winput::spillage_outdir << "/"
       << "F_delta.dat";
    if (GlobalV::MY_RANK == 0)
    {
        ofs.open(ss.str().c_str());
    }
    ofs << "F_delta(Hatree/Bohr) from deepks model: " << endl;
    ofs << setw(12) << "type" << setw(12) << "atom" << setw(15) << "dF_x" << setw(15) << "dF_y" << setw(15) << "dF_z" << endl;
    for (int it = 0;it < GlobalC::ucell.ntype;++it)
    {
        for (int ia = 0;ia < GlobalC::ucell.atoms[it].na;++ia)
        {
            int iat = GlobalC::ucell.itia2iat(it, ia);
            ofs << setw(12) << GlobalC::ucell.atoms[it].label << setw(12) << ia
                << setw(15) << this->F_delta(iat, 0) / 2 << setw(15) << this->F_delta(iat, 1) / 2
                << setw(15) << this->F_delta(iat, 2) / 2 << endl;
        }
    }
    ofs << "F_delta(eV/Angstrom) from deepks model: " << endl;
    ofs << setw(12) << "type" << setw(12) << "atom" << setw(15) << "dF_x" << setw(15) << "dF_y" << setw(15) << "dF_z" << endl;
    for (int it = 0;it < GlobalC::ucell.ntype;++it)
    {
        for (int ia = 0;ia < GlobalC::ucell.atoms[it].na;++ia)
        {
            int iat = GlobalC::ucell.itia2iat(it, ia);
            ofs << setw(12) << GlobalC::ucell.atoms[it].label << setw(12)
                << ia << setw(15) << this->F_delta(iat, 0) * Ry_to_eV/BOHR_TO_A
                << setw(15) << this->F_delta(iat, 1) * Ry_to_eV/BOHR_TO_A
                << setw(15) << this->F_delta(iat, 2) * Ry_to_eV/BOHR_TO_A << endl;
        }
    }
    GlobalV::ofs_running << "F_delta is printed" << endl;
    return;
}

void LCAO_Descriptor::save_npy_d()
{
    //save descriptor in .npy format
    vector<double> npy_des;
    for (int i = 0;i < this->n_descriptor;++i)
    {
        npy_des.push_back(this->d[i]);
    }
    const long unsigned dshape[] = {(long unsigned) GlobalC::ucell.nat, (long unsigned) this->des_per_atom };
    npy::SaveArrayAsNumpy("dm_eig.npy", false, 2, dshape, npy_des);
    return;
}

void LCAO_Descriptor::save_npy_e(double& ebase)
{   
    //save e_base
    const long unsigned eshape[] = { 1 };
    vector<double> npy_ebase;
    npy_ebase.push_back(ebase);
    npy::SaveArrayAsNumpy("e_base.npy", false, 1, eshape, npy_ebase);
    return;
}

void LCAO_Descriptor::save_npy_f(matrix& fbase)
{
    //save f_base
    //caution: unit: Rydberg/Bohr
    const long unsigned fshape[] = {(long unsigned) GlobalC::ucell.nat, 3 };
    vector<double> npy_fbase;
    for (int iat = 0;iat < GlobalC::ucell.nat;++iat)
    {
        for (int i = 0;i < 3;i++)
        {
            npy_fbase.push_back(fbase(iat, i));
        }
    }
    npy::SaveArrayAsNumpy("f_base.npy", false, 2, fshape, npy_fbase);
    return;
}

#endif
