//caoyu add 2021-03-2
#ifdef __DEEPKS

#include "LCAO_descriptor.h"

//===============================
//DeePKS Part 1
//deals with generation of descriptors as well as labels
//===============================

namespace GlobalC
{
    LCAO_Descriptor ld;
}

LCAO_Descriptor::LCAO_Descriptor()
{
    alpha_index = new ModuleBase::IntArray[1];
    inl_index = new ModuleBase::IntArray[1];
    inl_l = new int[1];
    d = new double[1];
    H_V_delta = new double[1];
    dm_double = new double[1];
    DH_V_delta_x = new double[1];
    DH_V_delta_y = new double[1];
    DH_V_delta_z = new double[1];
}

LCAO_Descriptor::~LCAO_Descriptor()
{
    delete[] alpha_index;
    delete[] inl_index;
    delete[] inl_l;
    delete[] d;
    delete[] H_V_delta;
    delete[] dm_double;
    delete[] DH_V_delta_x;
    delete[] DH_V_delta_y;
    delete[] DH_V_delta_z;

    //=======1. "out_descriptor" part==========
    //delete S_mu_alpha**
    for (int inl = 0;inl < this->inlmax;inl++)
    {
        delete[] S_mu_alpha[inl];
    }
    delete[] S_mu_alpha;
    //delete pdm**
    for (int inl = 0;inl < this->inlmax;inl++)
    {
        delete[] pdm[inl];
    }
    delete[] pdm;
    //=======2. "deepks_scf" part==========
    if (GlobalV::deepks_scf)
    {
        //delete gedm**
        for (int inl = 0;inl < this->inlmax;inl++)
        {
            delete[] gedm[inl];
        }
        delete[] gedm;
        if (GlobalV::FORCE)
        {
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
        }
    }
}

void LCAO_Descriptor::init(
	const int lm, // max L for descriptor 
	const int nm, // max n for descriptor
	const int tot_inl) // total number of atomic orbital basis for all of the atoms 
{
    ModuleBase::TITLE("LCAO_Descriptor", "init");

    GlobalV::ofs_running << " Initialize the descriptor index for DeePKS (lcao line)" << std::endl;

    assert(lm >= 0);
    assert(nm >= 0);
    assert(tot_inl >= 0);
    
    this->lmaxd = lm;
    this->nmaxd = nm;
    this->inlmax = tot_inl;
    GlobalV::ofs_running << " lmax of descriptor = " << this->lmaxd << std::endl;
    GlobalV::ofs_running << " nmax of descriptor= " << nmaxd << std::endl;
	GlobalV::ofs_running << " total basis (all atoms) for descriptor= " << std::endl;
    
    //initialize the density matrix (dm)
    delete[] this->dm_double;
    this->dm_double = new double[GlobalV::NLOCAL * GlobalV::NLOCAL];
    ModuleBase::GlobalFunc::ZEROS(this->dm_double, GlobalV::NLOCAL * GlobalV::NLOCAL);
    
    //init S_mu_alpha**
    this->S_mu_alpha = new double* [this->inlmax];    //inl
    for (int inl = 0;inl < this->inlmax;inl++)
    {
        this->S_mu_alpha[inl] = new double[GlobalV::NLOCAL * (2 * this->lmaxd + 1)];     //GlobalV::NLOCAL*nm
        ModuleBase::GlobalFunc::ZEROS(S_mu_alpha[inl], GlobalV::NLOCAL * (2 * this->lmaxd+ 1));
    }

    //init pdm**
    const int pdm_size = (this->lmaxd * 2 + 1) * (this->lmaxd * 2 + 1);
    this->pdm = new double* [this->inlmax];
    for (int inl = 0;inl < this->inlmax;inl++)
    {
        this->pdm[inl] = new double[pdm_size];
        ModuleBase::GlobalFunc::ZEROS(this->pdm[inl], pdm_size);
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
    this->alpha_index = new ModuleBase::IntArray[GlobalC::ucell.ntype];
    delete[] this->inl_index;
    this->inl_index = new ModuleBase::IntArray[GlobalC::ucell.ntype];
    delete[] this->inl_l;
    this->inl_l = new int[this->inlmax];
    ModuleBase::GlobalFunc::ZEROS(this->inl_l, this->inlmax);

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

        GlobalV::ofs_running << " Type " << it + 1
                    << " number_of_atoms " << GlobalC::ucell.atoms[it].na << std::endl;

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
    GlobalV::ofs_running << " descriptors_per_atom " << this->des_per_atom << std::endl;
    GlobalV::ofs_running << " total_descriptors " << this->n_descriptor << std::endl;
	return;
}

//this subroutine calculates the inner product between projectors and atomic basis
//<alpha|chi> as well as its derivative d/dX <alpha|chi>
//the former is recorded in array S_mu_alpha; the latter in arrays DS_mu_alpha_x,y,z
void LCAO_Descriptor::build_S_descriptor(const bool& calc_deri)
{
    ModuleBase::TITLE("LCAO_Descriptor", "build_S_descriptor");
    //array to store data

    double olm[3] = {0.0, 0.0, 0.0};
    ModuleBase::Vector3<double> tau1, tau2, dtau;
    ModuleBase::Vector3<double> dtau1, dtau2, tau0;
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

                        int iw1_local=GlobalC::ParaD.trace_loc_orb[iw1_all];
                        if(iw1_local<0)
                        {
                            ++iw1_all;
                        }
                        else
                        {
                            for (int L2 = 0; L2 <= GlobalC::ORB.Alpha[0].getLmax(); ++L2)
                            {
                                for (int N2 = 0; N2 < GlobalC::ORB.Alpha[0].getNchi(L2); ++N2)
                                {
                                    for (int m2 = 0; m2 < 2 * L2 + 1; ++m2)
                                    {
                                        olm[0] = olm[1] = olm[2] = 0.0;
                                        if (!calc_deri)
                                        {
                                            GlobalC::UOT.snap_psipsi(GlobalC::ORB, olm, 0, 'D', tau1,
                                                    T1, L1, m1, N1, GlobalC::GridD.getAdjacentTau(ad),
                                                    T2, L2, m2, N2, GlobalV::NSPIN);
                                            this->set_S_mu_alpha(iw1_all, inl_index[T2](I2,L2,N2), m2, olm[0]);
                                        }
                                        else
                                        {
                                            GlobalC::UOT.snap_psipsi(GlobalC::ORB, olm, 1, 'D', tau1,
                                                T1, L1, m1, N1, GlobalC::GridD.getAdjacentTau(ad),
                                                T2, L2, m2, N2, GlobalV::NSPIN);
                                                this->set_DS_mu_alpha(iw1_all, inl_index[T2](I2,L2,N2), m2, olm[0], olm[1], olm[2]);
                                        }

                                    } //m2
                                }     //N2
                            }         //nw2(L2)
                            ++iw1_all;
                        }
                    } // nw1
                }     // distance
            }         // ad
        } // I1
    }     // T1

#ifdef __MPI
    if(calc_deri)
    {
        GlobalC::ParaD.allsum_deepks(this->inlmax,GlobalV::NLOCAL*(2*this->lmaxd+1),this->DS_mu_alpha_x);
        GlobalC::ParaD.allsum_deepks(this->inlmax,GlobalV::NLOCAL*(2*this->lmaxd+1),this->DS_mu_alpha_y);
        GlobalC::ParaD.allsum_deepks(this->inlmax,GlobalV::NLOCAL*(2*this->lmaxd+1),this->DS_mu_alpha_z);
    }
    else
    {
        GlobalC::ParaD.allsum_deepks(this->inlmax,GlobalV::NLOCAL*(2*this->lmaxd+1),this->S_mu_alpha);
    }
#endif

/*    
    for(int inl=0;inl<this->inlmax;inl++)
    {
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"inl:",inl);
        for(int j=0;j<GlobalV::NLOCAL*(2*this->lmaxd+1);j++)
        {
            GlobalV::ofs_running << "j,s_mu_alpha: " << j << " " << this->S_mu_alpha[inl][j] << std::endl;
        }
    }
*/    
    return;
}

void LCAO_Descriptor::set_S_mu_alpha(
	const int &iw1_all, 
	const int &inl, 
	const int &im, 
	const double &v)
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
    this->S_mu_alpha[inl][index] += v;
    return;
}


//this subroutine performs the calculation of projected density matrices
//pdm_m,m'=\sum_{mu,nu} rho_{mu,nu} <chi_mu|alpha_m><alpha_m'|chi_nu>
void LCAO_Descriptor::cal_projected_DM(const ModuleBase::matrix &dm)
{
    ModuleBase::TITLE("LCAO_Descriptor", "cal_projected_DM");
    ModuleBase::timer::tick("LCAO_Descriptor","cal_projected_DM"); 
    const int pdm_size = (this->lmaxd * 2 + 1) * (this->lmaxd * 2 + 1);

    if(GlobalV::NPROC>1)
    {
#ifdef __MPI
        //This is for first SCF iteration, when density matrix is not available yet
        if(dm.nr == 0 && dm.nc ==0)
        {
            ModuleBase::timer::tick("LCAO_Descriptor","cal_projected_DM"); 
            return;
        }
        //step 1: get S_alpha_mu and S_nu_beta
        double **ss = this->S_mu_alpha;

        //step 2 : multiply: cal A=ST*DM*S
        //A(im1,im2) = sum iw1 sum iw2 S(iw1,im1) * dm(iw1,iw2) * S(iw2,im2)
        // = sum iw1 S(iw1,im1) * X(iw1,im2)
        // where X(iw1,im2) = sum iw2 dm(iw1,iw2) * S(iw2,im2)
        for (int inl = 0;inl < inlmax;inl++)
        {
            ModuleBase::GlobalFunc::ZEROS(this->pdm[inl], pdm_size);
            int nm = 2 * inl_l[inl] + 1;
            const int tmp_pdm_size = GlobalV::NLOCAL * nm;
            double* tmp_pdm = new double[tmp_pdm_size]; // saves X(iw1,im2)

            //for each pair index1=(iw1,im2)
            for(int iw1 = 0; iw1 < GlobalV::NLOCAL; iw1++)
            {
                int iw1_local = GlobalC::ParaO.trace_loc_col[iw1];
                if(iw1_local < 0) continue;
                const int ir1 = iw1;

                ModuleBase::GlobalFunc::ZEROS(tmp_pdm, tmp_pdm_size);

                for(int im2=0;im2<nm;im2++)
                {
                    const int ic2 = im2;

                    int index1;
                    if (GlobalV::KS_SOLVER == "genelpa" || GlobalV::KS_SOLVER == "scalapack_gvx") // save the matrix as column major format
                    {
                        index1 = ic2 * GlobalV::NLOCAL + ir1;
                    }
                    else
                    {
                        index1 = ir1 * nm  + ic2; //row: lcao orbitals; col: descriptor basis                        
                    }

                    //calculates X(iw1,im2) = sum iw2 dm(iw1,iw2) * S(iw2,im2)
                    for(int iw2=0; iw2 < GlobalV::NLOCAL; iw2++)
                    {
                        int iw2_local = GlobalC::ParaO.trace_loc_row[iw2];
                        if(iw2_local < 0) continue;
                        const int ir2 = iw2;

                        int index2;
                        if (GlobalV::KS_SOLVER == "genelpa" || GlobalV::KS_SOLVER == "scalapack_gvx") // save the matrix as column major format
                        {
                            index2 = ic2 * GlobalV::NLOCAL + ir2;
                        }
                        else
                        {
                            index2 = ir2 * nm  + ic2; //row: lcao orbitals; col: descriptor basis                        
                        }
                        double element = ss[inl][index2]* dm(iw1_local,iw2_local);
                        tmp_pdm[index1] += element;
                    }

                    //for each im1 : accumulates S(iw1,im1) * X(iw1,im2)
                    for(int im1=0;im1<nm;im1++)
                    {
                        const int ic1 = im1;
                        int index3;
                        if (GlobalV::KS_SOLVER == "genelpa" || GlobalV::KS_SOLVER == "scalapack_gvx") // save the matrix as column major format
                        {
                            index3 = ic1 * GlobalV::NLOCAL + ir1;
                        }
                        else
                        {
                            index3 = ir1 * nm  + ic1; //row: lcao orbitals; col: descriptor basis                        
                        }
                        double element = tmp_pdm[index1] * ss[inl][index3];
                        int ind = im1 + im2 * nm;
                        this->pdm[inl][ind] += element;
                    }
                }
            }
            delete[] tmp_pdm;
        }

        //step 3 : gather from all ranks
        GlobalC::ParaD.allsum_deepks(this->inlmax,pdm_size,this->pdm);
#endif
    }
    else //serial; or mpi with nproc=1
    {
// for checking purpose
/*
        std::ifstream ifs("dm");
        ModuleBase::matrix dm1(dm.nr,dm.nc);

        if(dm.nc>0)
        {
            for (int mu=0;mu<GlobalV::NLOCAL;mu++)
            {
                for (int nu=0;nu<GlobalV::NLOCAL;nu++)
                {
                    double c;
                    ifs >> c;
                    dm1(mu,nu)=c;
                }
            }
        }
*/      
//
        //step 1: get dm: the coefficient of wfc, not charge density
        //now,  dm is an input arg of this func, but needed converting to double*
        this->getdm_double(dm);

        //step 2: get S_alpha_mu and S_nu_beta
        double **ss = this->S_mu_alpha;

        //step 3 : multiply: cal ST*DM*S
        
        //init tmp_pdm*
        const int tmp_pdm_size = GlobalV::NLOCAL * (lmaxd*2+1);
        double* tmp_pdm = new double[tmp_pdm_size];
        ModuleBase::GlobalFunc::ZEROS(tmp_pdm, tmp_pdm_size);
        for (int inl = 0;inl < inlmax;inl++)
        {   
            int nm = 2 * inl_l[inl] + 1;   //1,3,5,...
            const char t = 'T';  //transpose
            const char nt = 'N'; //non transpose
            const double alpha = 1;
            const double beta = 0;
            double *a = this->dm_double;
            double *b = ss[inl];
            double *c = tmp_pdm;
            dgemm_(&nt, &nt, &GlobalV::NLOCAL, &nm, &GlobalV::NLOCAL, &alpha, a, &GlobalV::NLOCAL, b, &GlobalV::NLOCAL, &beta, c, &GlobalV::NLOCAL); //DM*S
            a = ss[inl];
            b = c;
            c = this->pdm[inl];
            dgemm_(&t, &nt, &nm, &nm, &GlobalV::NLOCAL, &alpha, a, &GlobalV::NLOCAL, b, &GlobalV::NLOCAL, &beta, c, &nm); //ST*DM*S
        }
        delete[] tmp_pdm;
    }

// for checking purpose
/*
    for(int inl=0;inl<inlmax;inl++) 
    {
        int dim = 2 * inl_l[inl] + 1;
        if(dim>1) continue; //print s orbitals for checking
        for (int m = 0; m < dim; m++)
        {
            for (int m2 = 0; m2 < dim; m2++)
            {
                int index = m * dim + m2;
                GlobalV::ofs_running << std::setprecision(10) << pdm[inl][index] << " ";
            }
        }
        GlobalV::ofs_running << std::endl;
    }
*/
//
    ModuleBase::timer::tick("LCAO_Descriptor","cal_projected_DM"); 
    return;
}

void LCAO_Descriptor::cal_projected_DM_k(const std::vector<ModuleBase::ComplexMatrix>& dm)
{

    std::complex<double> **pdm_complex;
    pdm_complex = new std::complex<double>* [this->inlmax];
    const int pdm_size = (this->lmaxd * 2 + 1) * (this->lmaxd * 2 + 1);
    for(int inl=0;inl<this->inlmax;inl++)
    {
        pdm_complex[inl] = new std::complex<double> [pdm_size];
        ModuleBase::GlobalFunc::ZEROS(pdm_complex[inl],pdm_size);
    }
    const double Rcut_Alpha = GlobalC::ORB.Alpha[0].getRcut();

    for (int T0 = 0; T0 < GlobalC::ucell.ntype; T0++)
    {
		Atom* atom0 = &GlobalC::ucell.atoms[T0]; 
        for (int I0 =0; I0< atom0->na; I0++)
        {
            const int iat = GlobalC::ucell.itia2iat(T0,I0);
            const ModuleBase::Vector3<double> tau0 = atom0->tau[I0];
            GlobalC::GridD.Find_atom(GlobalC::ucell, atom0->tau[I0] ,T0, I0);

            for (int ad1=0; ad1<GlobalC::GridD.getAdjacentNum()+1 ; ++ad1)
            {
                const int T1 = GlobalC::GridD.getType(ad1);
                const int I1 = GlobalC::GridD.getNatom(ad1);
                const int start1 = GlobalC::ucell.itiaiw2iwt(T1, I1, 0);
                const ModuleBase::Vector3<double> tau1 = GlobalC::GridD.getAdjacentTau(ad1);
				const Atom* atom1 = &GlobalC::ucell.atoms[T1];
				const int nw1_tot = atom1->nw*GlobalV::NPOL;
				const double Rcut_AO1 = GlobalC::ORB.Phi[T1].getRcut();

                ModuleBase::Vector3<double> dR1(GlobalC::GridD.getBox(ad1).x, GlobalC::GridD.getBox(ad1).y, GlobalC::GridD.getBox(ad1).z); 

				for (int ad2=0; ad2 < GlobalC::GridD.getAdjacentNum()+1 ; ad2++)
				{
					const int T2 = GlobalC::GridD.getType(ad2);
					const int I2 = GlobalC::GridD.getNatom(ad2);
					const int start2 = GlobalC::ucell.itiaiw2iwt(T2, I2, 0);
					const ModuleBase::Vector3<double> tau2 = GlobalC::GridD.getAdjacentTau(ad2);
					const Atom* atom2 = &GlobalC::ucell.atoms[T2];
					const int nw2_tot = atom2->nw*GlobalV::NPOL;
                    ModuleBase::Vector3<double> dR2(GlobalC::GridD.getBox(ad2).x, GlobalC::GridD.getBox(ad2).y, GlobalC::GridD.getBox(ad2).z);
					
					const double Rcut_AO2 = GlobalC::ORB.Phi[T2].getRcut();
                	const double dist1 = (tau1-tau0).norm() * GlobalC::ucell.lat0;
                	const double dist2 = (tau2-tau0).norm() * GlobalC::ucell.lat0;

					if (dist1 > Rcut_Alpha + Rcut_AO1
							|| dist2 > Rcut_Alpha + Rcut_AO2)
					{
						continue;
					}

					for (int iw1=0; iw1<nw1_tot; ++iw1)
					{
						const int iw1_all = start1 + iw1;
						const int iw1_local = GlobalC::ParaO.trace_loc_row[iw1_all];
						if(iw1_local < 0)continue;
						const int iw1_0 = iw1/GlobalV::NPOL;
						for (int iw2=0; iw2<nw2_tot; ++iw2)
						{
							const int iw2_all = start2 + iw2;
							const int iw2_local = GlobalC::ParaO.trace_loc_col[iw2_all];
							if(iw2_local < 0)continue;
							const int iw2_0 = iw2/GlobalV::NPOL;

                            std::vector<double> nlm1 = this->nlm_save[iat][ad1][iw1_all][0];
                            std::vector<double> nlm2 = this->nlm_save[iat][ad2][iw2_all][0];

                            assert(nlm1.size()==nlm2.size());
                            for(int ik=0;ik<GlobalC::kv.nks;ik++)
                            {
                                const double arg = ( GlobalC::kv.kvec_d[ik] * (dR2-dR1) ) * ModuleBase::TWO_PI;
                                const std::complex<double> kphase = std::complex <double> ( cos(arg),  sin(arg) );
                                int ib=0;
                                for (int L0 = 0; L0 <= GlobalC::ORB.Alpha[0].getLmax();++L0)
                                {
                                    for (int N0 = 0;N0 < GlobalC::ORB.Alpha[0].getNchi(L0);++N0)
                                    {
                                        const int inl = this->inl_index[T0](I0, L0, N0);
                                        const int nm = 2*L0+1;
                                        for (int m1 = 0;m1 < 2 * L0 + 1;++m1)
                                        {
                                            for (int m2 = 0; m2 < 2 * L0 + 1; ++m2)
                                            {
                                                pdm_complex[inl][m1*nm+m2] += nlm1[ib+m1]*nlm2[ib+m2]*dm[ik](iw2_local,iw1_local);
                                            }
                                        }
                                        ib+=nm;
                                    }
                                }
                                assert(ib==nlm1.size());
                            }//ik
						}//iw2
					}//iw1
				}//ad2
			}//ad1
        }//I0
    }//T0

    for(int inl=0;inl<inlmax;inl++)
    {
        for(int ind=0;ind<pdm_size;ind++)
        {
            if(pdm_complex[inl][ind].imag()>1.0e-8)
            {
                ModuleBase::WARNING_QUIT("pdm_k","pdm_complex not real!");
            }
            this->pdm[inl][ind] = pdm_complex[inl][ind].real();
        }
        delete[] pdm_complex[inl];
    }
    delete[] pdm_complex;
#ifdef __MPI
        GlobalC::ParaD.allsum_deepks(this->inlmax,pdm_size,this->pdm);
#endif

    return;
    
}

//the eigenvalues of pdm are descriptors
void LCAO_Descriptor::cal_descriptor(void)
{
    ModuleBase::TITLE("LCAO_Descriptor", "cal_descriptor");
    delete[] d;

    d = new double[this->n_descriptor];
    const int lmax = GlobalC::ORB.get_lmax_d();
    int id = 0;
    for (int it = 0; it < GlobalC::ucell.ntype; it++)
    {
        for (int ia = 0; ia < GlobalC::ucell.atoms[it].na; ia++)
        {
            for (int l = 0; l <= lmax; l++)
            {
                int nmax = GlobalC::ORB.Alpha[0].getNchi(l);
                for (int n = 0; n < nmax; n++)
                {
                    const int dim = 2 * l + 1;
                    const int inl = inl_index[it](ia, l, n);
                    // descriptor for atom (it, ia)
                    ModuleBase::ComplexMatrix des(dim, dim);
                    for (int m = 0; m < dim; m++)
                    {
                        for (int m2 = 0; m2 < dim; m2++)
                        {
                            int index = m * dim + m2;
                            complex<double> tmp(this->pdm[inl][index], 0);
                            des(m, m2) += tmp;
                        }
                    }
                    /*
                    stringstream ss;
                    ss << winput::spillage_outdir << "/"<< "PDM.dat";
                    if (GlobalV::MY_RANK == 0)
                    {
                        ofstream ofs;
                        ofs.open(ss.str().c_str());
                        this->print_projected_DM(ofs, des, it, ia, l, n);
                        ofs.close();
                    }
                    */
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
    this->print_descriptor();
    return;
}


//for checking purpose
void LCAO_Descriptor::print_projected_DM(
	ofstream& ofs, 
	ModuleBase::ComplexMatrix& des, 
	const int& it, // index for atom type
	const int& ia, // index for atoms
	const int& l, // index for angular momentum quantum number L
	const int& n) // index for principal quantum number n
{
    ofs << "L=" << l << "   N=" << n << std::endl;
    for (int i = 0; i < 2 * l + 1; i++)
    {
        for (int j = 0; j < 2 * l + 1; j++)
        {
            ofs << des(i, j).real() << " ";
        }
        ofs << std::endl;
    }
    return;
}


void LCAO_Descriptor::print_descriptor(void)
{
    ModuleBase::TITLE("LCAO_Descriptor", "print_descriptor");
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
            ofs << GlobalC::ucell.atoms[it].label << " atom_index " << ia + 1 << " n_descriptor " << this->des_per_atom << std::endl;
            int id0 = this->alpha_index[it](ia, 0, 0, 0);
            for (int id = id0; id < id0 + this->des_per_atom; ++id)
            {
                if ((id - id0) > 0 && (id - id0) % 8 == 0)
                    ofs << std::endl;
                ofs << std::setprecision(10) << d[id] << " ";
            }
            ofs << std::endl << std::endl;
        }
    }
    GlobalV::ofs_running << " Descriptors have been printed to " << ss.str() << std::endl;

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


void LCAO_Descriptor::getdm_double(const ModuleBase::matrix &dm)
{
    for (int i = 0; i < dm.nr; i++)
    {
        for (int j = 0; j < dm.nc; j++)
        {
            this->dm_double[i * GlobalV::NLOCAL + j] = dm(i, j); //only consider default GlobalV::NSPIN = 1
        }
    }
	return;
}

//this subroutine calculates the gradient of projected density matrices
//gdmx_m,m = d/dX sum_{mu,nu} rho_{mu,nu} <chi_mu|alpha_m><alpha_m'|chi_nu>
void LCAO_Descriptor::cal_gdmx(const ModuleBase::matrix &dm)
{
    ModuleBase::TITLE("LCAO_Descriptor", "cal_gdmx");
    //get DS_alpha_mu and S_nu_beta

// for checking purpose
/*    
    std::ifstream ifs("dm");
    ModuleBase::matrix dm1(dm.nr,dm.nc);

    for (int mu=0;mu<GlobalV::NLOCAL;mu++)
    {
        for (int nu=0;nu<GlobalV::NLOCAL;nu++)
        {
            double c;
            ifs >> c;
            dm1(mu,nu)=c;
        }
    }
*/
//
    double** ss = this->S_mu_alpha;
    double** dsx = this->DS_mu_alpha_x;
    double** dsy = this->DS_mu_alpha_y;
    double** dsz = this->DS_mu_alpha_z;
    int nlmax = inlmax/GlobalC::ucell.nat;

    for (int inl = 0;inl < inlmax;inl++)
    {
        //dE/dD will be multiplied in cal_f_delta, here only calculate dD/dx_I
        int nm = 2 * inl_l[inl] + 1;   //1,3,5,...
        const int ibt = inl/nlmax; //the atom on which alpha is located

        for (int mu =0; mu < GlobalV::NLOCAL;++mu) 
        {
            const int iat = GlobalC::ucell.iwt2iat[mu];//the atom on which chi_mu is located
            int iw1 = GlobalC::ParaO.trace_loc_col[mu];
            if(iw1 < 0) continue;
            
            for (int nu= 0;nu < GlobalV::NLOCAL; ++nu)
            {
                int iw2 = GlobalC::ParaO.trace_loc_row[nu];
                if(iw2 < 0) continue;

                for (int m1 = 0;m1 < nm;++m1)
                {
                    for (int m2 = 0;m2 < nm;++m2)
                    {
                        for (int is = 0;is < GlobalV::NSPIN;++is)
                        {
                            //  save the matrix as column major format
                            if (GlobalV::KS_SOLVER == "genelpa" || GlobalV::KS_SOLVER == "scalapack_gvx")
                            {
                                //mu->iw1,nu->iw2
                                //(<d/dX chi_mu|alpha_m>)<chi_nu|alpha_m'>
                                gdmx[iat][inl][m1*nm + m2] += 
                                dsx[inl][m1*GlobalV::NLOCAL + mu] * dm(iw1, iw2) * ss[inl][m2*GlobalV::NLOCAL + nu];
                                gdmy[iat][inl][m1*nm + m2] += 
                                dsy[inl][m1*GlobalV::NLOCAL + mu] * dm(iw1, iw2) * ss[inl][m2*GlobalV::NLOCAL + nu];
                                gdmz[iat][inl][m1*nm + m2] += 
                                dsz[inl][m1*GlobalV::NLOCAL + mu] * dm(iw1, iw2) * ss[inl][m2*GlobalV::NLOCAL + nu];

                                //(<d/dX chi_nu|alpha_m'>)<chi_mu|alpha_m>
                                gdmx[iat][inl][m2*nm + m1] += 
                                dsx[inl][m1*GlobalV::NLOCAL + mu] * dm(iw1, iw2) * ss[inl][m2*GlobalV::NLOCAL + nu];
                                gdmy[iat][inl][m2*nm + m1] += 
                                dsy[inl][m1*GlobalV::NLOCAL + mu] * dm(iw1, iw2) * ss[inl][m2*GlobalV::NLOCAL + nu];
                                gdmz[iat][inl][m2*nm + m1] += 
                                dsz[inl][m1*GlobalV::NLOCAL + mu] * dm(iw1, iw2) * ss[inl][m2*GlobalV::NLOCAL + nu];

                                //(<chi_mu|d/dX alpha_m>)<chi_nu|alpha_m'> = -(<d/dX chi_mu|alpha_m>)<chi_nu|alpha_m'>
                                gdmx[ibt][inl][m1*nm + m2] -= 
                                dsx[inl][m1*GlobalV::NLOCAL + mu] * dm(iw1, iw2) * ss[inl][m2*GlobalV::NLOCAL + nu];
                                gdmy[ibt][inl][m1*nm + m2] -= 
                                dsy[inl][m1*GlobalV::NLOCAL + mu] * dm(iw1, iw2) * ss[inl][m2*GlobalV::NLOCAL + nu];
                                gdmz[ibt][inl][m1*nm + m2] -= 
                                dsz[inl][m1*GlobalV::NLOCAL + mu] * dm(iw1, iw2) * ss[inl][m2*GlobalV::NLOCAL + nu];

                                //(<chi_nu|d/dX alpha_m'>)<chi_mu|alpha_m> = -(<d/dX chi_nu|alpha_m'>)<chi_mu|alpha_m>
                                gdmx[ibt][inl][m2*nm + m1] -= 
                                dsx[inl][m1*GlobalV::NLOCAL + mu] * dm(iw1, iw2) * ss[inl][m2*GlobalV::NLOCAL + nu];
                                gdmy[ibt][inl][m2*nm + m1] -= 
                                dsy[inl][m1*GlobalV::NLOCAL + mu] * dm(iw1, iw2) * ss[inl][m2*GlobalV::NLOCAL + nu];
                                gdmz[ibt][inl][m2*nm + m1] -= 
                                dsz[inl][m1*GlobalV::NLOCAL + mu] * dm(iw1, iw2) * ss[inl][m2*GlobalV::NLOCAL + nu];                                
                            }
                            else
                            {
                                gdmx[iat][inl][m1*nm + m2] += 
                                dsx[inl][mu*nm + m1] * dm(iw1, iw2) * ss[inl][nu*nm + m2];
                                gdmy[iat][inl][m1*nm + m2] += 
                                dsy[inl][mu*nm + m1] * dm(iw1, iw2) * ss[inl][nu*nm + m2];
                                gdmz[iat][inl][m1*nm + m2] += 
                                dsz[inl][mu*nm + m1] * dm(iw1, iw2) * ss[inl][nu*nm + m2];
                                
                                gdmx[iat][inl][m2*nm + m1] += 
                                dsx[inl][mu*nm + m1] * dm(iw1, iw2) * ss[inl][nu*nm + m2];
                                gdmy[iat][inl][m2*nm + m1] += 
                                dsy[inl][mu*nm + m1] * dm(iw1, iw2) * ss[inl][nu*nm + m2];
                                gdmz[iat][inl][m2*nm + m1] += 
                                dsz[inl][mu*nm + m1] * dm(iw1, iw2) * ss[inl][nu*nm + m2];

                                gdmx[ibt][inl][m1*nm + m2] -= 
                                dsx[inl][mu*nm + m1] * dm(iw1, iw2) * ss[inl][nu*nm + m2];
                                gdmy[ibt][inl][m1*nm + m2] -= 
                                dsy[inl][mu*nm + m1] * dm(iw1, iw2) * ss[inl][nu*nm + m2];
                                gdmz[ibt][inl][m1*nm + m2] -= 
                                dsz[inl][mu*nm + m1] * dm(iw1, iw2) * ss[inl][nu*nm + m2];
                                
                                gdmx[ibt][inl][m2*nm + m1] -= 
                                dsx[inl][mu*nm + m1] * dm(iw1, iw2) * ss[inl][nu*nm + m2];
                                gdmy[ibt][inl][m2*nm + m1] -= 
                                dsy[inl][mu*nm + m1] * dm(iw1, iw2) * ss[inl][nu*nm + m2];
                                gdmz[ibt][inl][m2*nm + m1] -= 
                                dsz[inl][mu*nm + m1] * dm(iw1, iw2) * ss[inl][nu*nm + m2];
                            }
                        }
                    }//end m2
                } //end m1
            }//end nu
        }//end mu     
    }//end inl
#ifdef __MPI
    const int gdm_size = (this->lmaxd * 2 + 1) * (this->lmaxd * 2 + 1);
    for(int iat=0;iat<GlobalC::ucell.nat;iat++)
    {
        GlobalC::ParaD.allsum_deepks(this->inlmax,gdm_size,this->gdmx[iat]);
        GlobalC::ParaD.allsum_deepks(this->inlmax,gdm_size,this->gdmy[iat]);
        GlobalC::ParaD.allsum_deepks(this->inlmax,gdm_size,this->gdmz[iat]);
    }
#endif   
// for checking purpose
/*
    GlobalV::ofs_running << "gdmx" << std::endl;
    for(int iat=0;iat<GlobalC::ucell.nat;iat++)
    {
        GlobalV::ofs_running << iat << std::endl;
        for(int inl = 0;inl < inlmax;inl++)
        {
            int nm = 2 * inl_l[inl] + 1;
            if(nm>1) continue; //print s orbitals for checking
            for (int m1 = 0;m1 < nm;++m1)
            {
                for (int m2 = 0;m2 < nm;++m2)
                {
                    int index = m1 * nm + m2;
                    GlobalV::ofs_running << std::setprecision(10) << gdmx[iat][inl][index] << " ";
                }
            }
            GlobalV::ofs_running << std::endl;
        }
    }
    GlobalV::ofs_running << "gdmy" << std::endl;
    for(int iat=0;iat<GlobalC::ucell.nat;iat++)
    {
        GlobalV::ofs_running << iat << std::endl;
        for(int inl = 0;inl < inlmax;inl++)
        {
            int nm = 2 * inl_l[inl] + 1;
            if(nm>1) continue; //print s orbitals for checking
            for (int m1 = 0;m1 < nm;++m1)
            {
                for (int m2 = 0;m2 < nm;++m2)
                {
                    int index = m1 * nm + m2;
                    GlobalV::ofs_running << std::setprecision(10) << gdmy[iat][inl][index] << " ";
                }
            }
            GlobalV::ofs_running << std::endl;
        }
    }
    GlobalV::ofs_running << "gdmz" << std::endl;
    for(int iat=0;iat<GlobalC::ucell.nat;iat++)
    {
        GlobalV::ofs_running << iat << std::endl;
        for(int inl = 0;inl < inlmax;inl++)
        {
            int nm = 2 * inl_l[inl] + 1;
            if(nm>1) continue; //print s orbitals for checking
            for (int m1 = 0;m1 < nm;++m1)
            {
                for (int m2 = 0;m2 < nm;++m2)
                {
                    int index = m1 * nm + m2;
                    GlobalV::ofs_running << std::setprecision(10) << gdmz[iat][inl][index] << " ";
                }
            }
            GlobalV::ofs_running << std::endl;
        }
    }
*/
//
    return;
}


void LCAO_Descriptor::init_gdmx(void)
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
            ModuleBase::GlobalFunc::ZEROS(gdmx[iat][inl], (2 * lmaxd + 1) * (2 * lmaxd + 1));
            ModuleBase::GlobalFunc::ZEROS(gdmy[iat][inl], (2 * lmaxd + 1) * (2 * lmaxd + 1));
            ModuleBase::GlobalFunc::ZEROS(gdmz[iat][inl], (2 * lmaxd + 1) * (2 * lmaxd + 1));
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

#endif
