#include "../src_parallel/parallel_common.h"
#include "src_parallel/parallel_orbitals.h"
#include "src_pdiag/pdiag_double.h"
#include "../src_pw/global.h"
#include "../src_io/wf_local.h"
#include "../module_base/lapack_connector.h"
#include "../module_base/memory.h"


void ORB_control::set_parameters(void)
{
    ModuleBase::TITLE("ORB_control","set_parameters");

    Parallel_Orbitals* pv = &this->ParaV;
    // set loc_size
	if(GlobalV::GAMMA_ONLY_LOCAL)//xiaohui add 2014-12-21
	{
		pv->loc_size=GlobalV::NBANDS/GlobalV::DSIZE;

		// mohan add 2012-03-29
		if(pv->loc_size==0)
		{
			GlobalV::ofs_warning << " loc_size=0" << " in proc " << GlobalV::MY_RANK+1 << std::endl;
			ModuleBase::WARNING_QUIT("ORB_control::set_parameters","NLOCAL < GlobalV::DSIZE");
		}

		if (GlobalV::DRANK<GlobalV::NBANDS%GlobalV::DSIZE) pv->loc_size+=1;
		if(pv->testpb)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"local size",pv->loc_size);

		// set loc_sizes
		delete[] pv->loc_sizes;
		pv->loc_sizes = new int[GlobalV::DSIZE];
		ModuleBase::GlobalFunc::ZEROS(pv->loc_sizes, GlobalV::DSIZE);

		pv->lastband_in_proc = 0;
		pv->lastband_number = 0;
		int count_bands = 0;
		for (int i=0; i<GlobalV::DSIZE; i++)
		{
			if (i<GlobalV::NBANDS%GlobalV::DSIZE)
			{
				// mohan modify 2010-07-05
				pv->loc_sizes[i]=GlobalV::NBANDS/GlobalV::DSIZE+1;
			}
			else
			{
				pv->loc_sizes[i]=GlobalV::NBANDS/GlobalV::DSIZE;
			}
			count_bands += pv->loc_sizes[i];
			if (count_bands >= GlobalV::NBANDS)
			{
				pv->lastband_in_proc = i;
				pv->lastband_number = GlobalV::NBANDS - (count_bands - pv->loc_sizes[i]);
				break;
			}
		}
	}
	else
	{
		pv->loc_size=GlobalV::NLOCAL/GlobalV::DSIZE;

		// mohan add 2012-03-29
		if(pv->loc_size==0)
		{
			GlobalV::ofs_warning << " loc_size=0" << " in proc " << GlobalV::MY_RANK+1 << std::endl;
			ModuleBase::WARNING_QUIT("ORB_control::set_parameters","NLOCAL < GlobalV::DSIZE");
		}

		if (GlobalV::DRANK<GlobalV::NLOCAL%GlobalV::DSIZE) 
		{
            pv->loc_size += 1;
        }
		if(pv->testpb) ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"local size",pv->loc_size);

		// set loc_sizes
		delete[] pv->loc_sizes;
		pv->loc_sizes = new int[GlobalV::DSIZE];
		ModuleBase::GlobalFunc::ZEROS(pv->loc_sizes, GlobalV::DSIZE);

		pv->lastband_in_proc = 0;
		pv->lastband_number = 0;
		int count_bands = 0;
		for (int i=0; i<GlobalV::DSIZE; i++)
		{
			if (i<GlobalV::NLOCAL%GlobalV::DSIZE)
			{
				// mohan modify 2010-07-05
				pv->loc_sizes[i]=GlobalV::NLOCAL/GlobalV::DSIZE+1;
			}
			else
			{
				pv->loc_sizes[i]=GlobalV::NLOCAL/GlobalV::DSIZE;
			}
			count_bands += pv->loc_sizes[i];
			if (count_bands >= GlobalV::NBANDS)
			{
				pv->lastband_in_proc = i;
				pv->lastband_number = GlobalV::NBANDS - (count_bands - pv->loc_sizes[i]);
				break;
			}
		}
	}//xiaohui add 2014-12-21

    if (GlobalV::KS_SOLVER=="hpseps") //LiuXh add 2021-09-06, clear memory, Z_LOC only used in hpseps solver
    {
	    pv->Z_LOC = new double*[GlobalV::NSPIN];
	    for(int is=0; is<GlobalV::NSPIN; is++)
	    {
		    pv->Z_LOC[is] = new double[pv->loc_size * GlobalV::NLOCAL];
		    ModuleBase::GlobalFunc::ZEROS(pv->Z_LOC[is], pv->loc_size * GlobalV::NLOCAL);
	    }
	    pv->alloc_Z_LOC = true;//xiaohui add 2014-12-22
    }

    if(pv->testpb)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"lastband_in_proc", pv->lastband_in_proc);
    if(pv->testpb)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"lastband_number", pv->lastband_number);

    return;
}


#ifdef __MPI
// creat the 'comm_2D' stratege.
void ORB_control::mpi_creat_cart(MPI_Comm *comm_2D, int prow, int pcol)
{
    ModuleBase::TITLE("ORB_control","mpi_creat_cart");
    // the matrix is divided as ( dim[0] * dim[1] )
    int dim[2];
    int period[2]={1,1};
    int reorder=0;
    dim[0]=prow;
    dim[1]=pcol;

    if(this->ParaV.testpb)GlobalV::ofs_running << " dim = " << dim[0] << " * " << dim[1] << std::endl;

    MPI_Cart_create(DIAG_WORLD,2,dim,period,reorder,comm_2D);
    return;
}
#endif

#ifdef __MPI
void ORB_control::mat_2d(MPI_Comm vu,
                         const int &M_A,
                         const int &N_A,
                         const int &nb,
                         LocalMatrix &LM)
{
    ModuleBase::TITLE("ORB_control", "mat_2d");
    
    Parallel_Orbitals* pv = &this->ParaV;
    
    int dim[2];
    int period[2];
    int coord[2];
    int i,j,k,end_id;
    int block;

    // (0) every processor get it's id on the 2D comm
    // : ( coord[0], coord[1] )
    MPI_Cart_get(vu,2,dim,period,coord);

    // (1.1) how many blocks at least
    // eg. M_A = 6400, nb = 64;
    // so block = 10;
    block=M_A/nb;

    // (1.2) If data remain, add 1.
    if (block*nb<M_A)
    {
        block++;
    }

    if(pv->testpb)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"Total Row Blocks Number",block);

	// mohan add 2010-09-12
	if(dim[0]>block)
	{
		GlobalV::ofs_warning << " cpu 2D distribution : " << dim[0] << "*" << dim[1] << std::endl;
		GlobalV::ofs_warning << " but, the number of row blocks is " << block << std::endl;
		ModuleBase::WARNING_QUIT("ORB_control::mat_2d","some processor has no row blocks, try a smaller 'nb2d' parameter.");
	}

    // (2.1) row_b : how many blocks for this processor. (at least)
    LM.row_b=block/dim[0];

    // (2.2) row_b : how many blocks in this processor.
    // if there are blocks remain, some processors add 1.
    if (coord[0]<block%dim[0])
    {
        LM.row_b++;
    }

    if(pv->testpb)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"Local Row Block Number",LM.row_b);

    // (3) end_id indicates the last block belong to
    // which processor.
    if (block%dim[0]==0)
    {
        end_id=dim[0]-1;
    }
    else
    {
        end_id=block%dim[0]-1;
    }

    if(pv->testpb)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"Ending Row Block in processor",end_id);

    // (4) row_num : how many rows in this processors :
    // the one owns the last block is different.
    if (coord[0]==end_id)
    {
        LM.row_num=(LM.row_b-1)*nb+(M_A-(block-1)*nb);
    }
    else
    {
        LM.row_num=LM.row_b*nb;
    }

    if(pv->testpb)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"Local rows (including nb)",LM.row_num);

    // (5) row_set, it's a global index :
    // save explicitly : every row in this processor
    // belongs to which row in the global matrix.
    delete[] LM.row_set;
    LM.row_set= new int[LM.row_num];
    j=0;
    for (i=0; i<LM.row_b; i++)
    {
        for (k=0; k<nb&&(coord[0]*nb+i*nb*dim[0]+k<M_A); k++,j++)
        {
            LM.row_set[j]=coord[0]*nb+i*nb*dim[0]+k;
           // GlobalV::ofs_running << " j=" << j << " row_set=" << LM.row_set[j] << std::endl;
        }
    }

    // the same procedures for columns.
    block=N_A/nb;
    if (block*nb<N_A)
    {
        block++;
    }

    if(pv->testpb)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"Total Col Blocks Number",block);

	if(dim[1]>block)
	{
		GlobalV::ofs_warning << " cpu 2D distribution : " << dim[0] << "*" << dim[1] << std::endl;
		GlobalV::ofs_warning << " but, the number of column blocks is " << block << std::endl;
		ModuleBase::WARNING_QUIT("ORB_control::mat_2d","some processor has no column blocks.");
	}

    LM.col_b=block/dim[1];
    if (coord[1]<block%dim[1])
    {
        LM.col_b++;
    }

    if(pv->testpb)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"Local Row Block Number",LM.col_b);

    if (block%dim[1]==0)
    {
        end_id=dim[1]-1;
    }
    else
    {
        end_id=block%dim[1]-1;
    }

    if(pv->testpb)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"Ending Row Block in processor",end_id);

    if (coord[1]==end_id)
    {
        LM.col_num=(LM.col_b-1)*nb+(N_A-(block-1)*nb);
    }
    else
    {
        LM.col_num=LM.col_b*nb;
    }

    if(pv->testpb)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"Local columns (including nb)",LM.row_num);

    delete[] LM.col_set;
    LM.col_set = new int[LM.col_num];

    j=0;
    for (i=0; i<LM.col_b; i++)
    {
        for (k=0; k<nb&&(coord[1]*nb+i*nb*dim[1]+k<N_A); k++,j++)
        {
            LM.col_set[j]=coord[1]*nb+i*nb*dim[1]+k;
        }
    }
    LM.col_pos=0;
    LM.row_pos=0;
    return;
}
#endif


#ifdef __MPI
// A : contains total matrix element in processor.
void ORB_control::data_distribution(
    MPI_Comm comm_2D,
    const std::string &file,
    const int &n,
    const int &nb,
    double *A,
    const LocalMatrix &LM)
{
    ModuleBase::TITLE("ORB_control", "data_distribution");
    Parallel_Orbitals* pv = &this->ParaV;
    MPI_Comm comm_row;
    MPI_Comm comm_col;
    MPI_Status status;

    int dim[2];
    int period[2];
    int coord[2];
    MPI_Cart_get(comm_2D,2,dim,period,coord);

    if(pv->testpb) GlobalV::ofs_running << "\n dim = " << dim[0] << " * " << dim[1] << std::endl;
    if(pv->testpb) GlobalV::ofs_running << " coord = ( " << coord[0] << " , " << coord[1] << ")." << std::endl;
    if(pv->testpb) GlobalV::ofs_running << " n = " << n << std::endl;

    mpi_sub_col(comm_2D,&comm_col);
    mpi_sub_row(comm_2D,&comm_row);

    // total number of processors
    const int myid = coord[0]*dim[1]+coord[1];

    // the matrix is n * n
    double* ele_val = new double[n];
    double* val = new double[n];
    int* sends = new int[dim[1]];
    int* fpt = new int[dim[1]];
    int* snd = new int[dim[1]];
    int* temp = new int[dim[1]];

    ModuleBase::GlobalFunc::ZEROS(ele_val, n);
    ModuleBase::GlobalFunc::ZEROS(val, n);
    ModuleBase::GlobalFunc::ZEROS(sends, dim[1]);
    ModuleBase::GlobalFunc::ZEROS(fpt, dim[1]);
    ModuleBase::GlobalFunc::ZEROS(snd, dim[1]);
    ModuleBase::GlobalFunc::ZEROS(temp, dim[1]);

    // the columes of matrix is divided by 'dim[1]' 'rows of processors'.
    // collect all information of each 'rows of processors'
    // collection data is saved in 'sends'
    snd[coord[1]] = LM.col_num;
    MPI_Allgather(&snd[coord[1]],1,MPI_INT,sends,1,MPI_INT,comm_row);

    // fpt : start column index after applied 'mat_2d' reorder algorithms
    // to matrix.
    fpt[0] = 0;
    for (int i=1; i<dim[1]; i++)
    {
        fpt[i]=fpt[i-1]+sends[i-1];
//      GlobalV::ofs_running << " col_pro = " << i << " start_col = " << fpt[i] << std::endl;
    }

//    GlobalV::ofs_running << "\n myid = " << myid << std::endl;

    int cur_i = 0;

    int iacol;
    int iarow;
    int ai;
    int aj;
    int tag = 0;

    bool find = true;
    if (myid==0)
    {
        FILE *fp;
        fp=fopen(file.c_str(),"rb");
        if (fp==NULL)
        {
            std::cout << " Can't find file : " << file << std::endl;
            find = false;
        }
        else
        {
            GlobalV::ofs_running << " Open file : " << file << std::endl;
            int dim = 0;
            fread(&dim,sizeof(int),1,fp);
            if (dim!=n)
            {
                find = false;
            }
            GlobalV::ofs_running << " Read in dimension = " << dim << std::endl;
        }
        int nrow = 0;
        while (nrow<n && !feof(fp))
        {
            ModuleBase::GlobalFunc::ZEROS(ele_val, n);
            ModuleBase::GlobalFunc::ZEROS(val, n);

            // read om one row elements.
//            GlobalV::ofs_running << "\n nrow = " << nrow << std::endl;

            for (int i=nrow; i<n; i++)
            {
                //if ((i-nrow)%8==0)GlobalV::ofs_running << std::endl;
                fread(&ele_val[i],sizeof(double),1,fp);
                //			GlobalV::ofs_running << " " << ele_val[i];
            }

            // start position of col_pro.
            for (int i=0; i<dim[1]; i++)
            {
                temp[i] = fpt[i];
            }

            for (int k=0; k<n; k++)
            {
                // calculate iarow and iacol.
                // belong to which col_pro.
                indxg2p(comm_2D,nb,nrow,k,&iarow,&iacol);
                val[temp[iacol]]=ele_val[k];
                temp[iacol]++;
            }

            indxg2l(nrow,0,nb,dim[0],dim[1],&ai,&aj);
            indxg2p(comm_2D,nb,nrow,0,&iarow,&iacol);

            const int incx = 1;
            if (iarow==0&&iacol==0)
            {
                BlasConnector::copy(LM.col_num,val,incx,&A[ai*LM.col_num],incx);
                for (int i=1; i<dim[1]; i++)
                {
//					GlobalV::ofs_running << " send to processor " << iarow*dim[1]+i << std::endl;
                    MPI_Send(&val[fpt[i]],sends[i],MPI_DOUBLE,iarow*dim[1]+i,tag,DIAG_WORLD);
                }
            }
            else
            {
                for (int i=0; i<dim[1]; i++)
                {
//					GlobalV::ofs_running << " else, send to processor " << iarow*dim[1]+i << std::endl;
                    MPI_Send(&val[fpt[i]],sends[i],MPI_DOUBLE,iarow*dim[1]+i,tag,DIAG_WORLD);
                }
            }
            nrow++;
        }// end read in nrow

        fclose(fp);
    }
    else
    {
        for (int j=0; j<LM.row_num; j++)
        {
//			GlobalV::ofs_running << " receive row = " <<  j << std::endl;
            MPI_Recv(&A[j*LM.col_num],LM.col_num,MPI_DOUBLE,0,tag,DIAG_WORLD,&status);
        }
    }

    /*
    for (int i=0; i<LM.row_num; i++)
    {
        GlobalV::ofs_running << "\n\n Row = " << i << std::endl;
        for (int j=0; j<LM.col_num; j++)
        {
            if (j%8==0) GlobalV::ofs_running << std::endl;
            GlobalV::ofs_running << " " << A[j*LM.col_num+i];
        }
    }
    */

    MPI_Comm_free(&comm_row);
    MPI_Comm_free(&comm_col);

    delete[] ele_val;
    delete[] val;
    delete[] sends;
    delete[] fpt;
    delete[] snd;
    delete[] temp;

#ifdef __MPI
    Parallel_Common::bcast_bool(find);
#endif

    //ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"Find the H/S file",find);

    if (!find)
    {
        ModuleBase::WARNING_QUIT("ORB_control::data_distribution","Can't find the H/S file");
    }

    return;
}
#endif

#ifdef __MPI
#include "../src_pw/occupy.h"
void Pdiag_Double::gath_eig_complex(MPI_Comm comm,int n,std::complex<double> **cc,std::complex<double> *Z, const int &ik)
{
    ModuleBase::TITLE("Pdiag_Double","gath_eig_complex");
    time_t time_start = time(NULL);
    //GlobalV::ofs_running << " Start gath_eig_complex Time : " << ctime(&time_start);
    const Parallel_Orbitals* pv = this->ParaV;

    int i,j,k;
    int nprocs,myid;
    MPI_Status status;
    MPI_Comm_size(comm,&nprocs);
    MPI_Comm_rank(comm,&myid);

    std::complex<double> **ctot;

	// mohan add 2010-07-03
	// the occupied bands are useless
	// for calculating charge density.
	if(GlobalV::DRANK> pv->lastband_in_proc)
	{
		delete[] Z;
	}

	// first we need to collect all
	// the occupied bands.
	// GlobalV::NBANDS * GlobalV::NLOCAL	
	if(GlobalV::DRANK==0)
	{
		ctot = new std::complex<double>*[GlobalV::NBANDS];
    	for (int i=0; i<GlobalV::NBANDS; i++)
    	{
        	ctot[i] = new std::complex<double>[GlobalV::NLOCAL];
        	ModuleBase::GlobalFunc::ZEROS(ctot[i], GlobalV::NLOCAL);
    	}
    	ModuleBase::Memory::record("Pdiag_Double","ctot",GlobalV::NBANDS*GlobalV::NLOCAL,"cdouble");
	}

	k=0;
    if (myid==0)
    {
        // mohan add nbnd0 2010-07-02
        int nbnd0 = -1;
        if (GlobalV::NBANDS < pv->loc_sizes[0])
        {
			// means all bands in this processor
			// is needed ( is occupied)
            nbnd0 = GlobalV::NBANDS;
        }
        else
        {
			// means this processor only save
			// part of GlobalV::NBANDS.
            nbnd0 = pv->loc_sizes[0];
        }
        if(pv->testpb)GlobalV::ofs_running << " nbnd in processor 0 is " << nbnd0 << std::endl;

        for (i=0; i<nbnd0; i++)
        {
            for (j=0; j<GlobalV::NLOCAL; j++)
            {
				// change the order in processor 0.
				// the contribution from processor 0.
                ctot[k][j]=Z[j*pv->loc_sizes[0]+i];
            }
            k++;
        }
		// Z is useless in processor 0 now.
		delete[] Z;
    }
    MPI_Barrier(comm);

	for (i=1; i<= pv->lastband_in_proc; i++)
    {
        // mohan fix bug 2010-07-02
        // rows indicates the data structure of Z.
        // mpi_times indicates the data distribution
        // time, each time send a band.
        int rows = pv->loc_sizes[i];
        int mpi_times;
        if (i==pv->lastband_in_proc)
        {
            mpi_times = pv->lastband_number;
        }
        else
        {
            mpi_times = pv->loc_sizes[i];
        }
        if(pv->testpb)GlobalV::ofs_running << " nbnd in processor " << i << " is " << mpi_times << std::endl;
        if (myid==i)
        {
            for (j=0; j<mpi_times; j++)
            {
                int tag = j;
                std::complex<double> *send = new std::complex<double>[n];
                int count = 0;

                for (int m=0; m<rows*n; m++)
                {
                    if (m%rows==j)
                    {
                        send[count] = Z[m];
                        ++count;
                    }
                }

				// send the data to processor 0.
                MPI_Send(send,n,mpicomplex,0,tag,comm);

                delete[] send;
            }
			// third part to delete Z;
			delete[] Z;
        }
        else if (myid==0)
        {
            int col=0;
            for (j=0; j<mpi_times; j++)
            {
                std::complex<double> *ctmp = new std::complex<double>[GlobalV::NLOCAL];
                ModuleBase::GlobalFunc::ZEROS(ctmp, GlobalV::NLOCAL);
                int tag = j;
                
				// Processor 0 receive the data from other processors.
				MPI_Recv(ctmp,n,mpicomplex,i,tag,comm,&status);

                for (int m=0; m<GlobalV::NLOCAL; m++)
                {
                    ctot[k][m]=ctmp[m];
//					GlobalV::ofs_running << " receive Z=" << ctmp[m] << std::endl;
                }
                k++;

                delete[] ctmp;
            }
        }
        //MPI_Barrier(comm);
    }
    if(pv->testpb)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"Final k",k);

	// output the wave function if required.
	// this is a bad position to output wave functions.
	// but it works!
	std::stringstream ss;
	ss << GlobalV::global_out_dir << "LOWF_K_" << ik+1 << ".dat";
    if(this->out_lowf)
	{
//		std::cout << " write the wave functions" << std::endl;
		WF_Local::write_lowf_complex( ss.str(), ctot, ik );//mohan add 2010-09-09        
	}

	// mohan add 2010-09-10
	// distribution of local wave functions 
	// to each processor.
	WF_Local::distri_lowf_complex( ctot, cc);
	
	// clean staff.
	if(GlobalV::DRANK==0)
	{
    	for (int i=0; i<GlobalV::NBANDS; i++)
    	{
        	delete[] ctot[i];
    	}
    	delete[] ctot;
	}

    time_t time_end = time(NULL);
//    GlobalV::ofs_running << " End   gath_eig_complex Time : " << ctime(&time_end);

	ModuleBase::GlobalFunc::OUT_TIME("gather eigenvalues",time_start,time_end);
    return;
}
#endif

#ifdef __MPI
void Pdiag_Double::gath_full_eig(MPI_Comm comm,int n,double **c,double *Z)
{
    ModuleBase::TITLE("Pdiag_Double","gath_full_eig");

    time_t time_start = time(NULL);
    //GlobalV::ofs_running << " Start gath_full_eig Time : " << ctime(&time_start);

    int i,j,k,incx=1;
    int *loc_sizes,loc_size,nprocs,myid;
    MPI_Status status;
    MPI_Comm_size(comm,&nprocs);
    MPI_Comm_rank(comm,&myid);
    loc_sizes =(int*)malloc(sizeof(int)*nprocs);
    loc_size=n/nprocs;
    for (i=0; i<nprocs; i++)
    {
        if (i<n%nprocs)
        {
            loc_sizes[i]=loc_size+1;
        }
        else
        {
            loc_sizes[i]=loc_size;
        }
    }
    if (myid==0)
    {
        k=0;
        for (i=0; i<loc_sizes[0]; i++)
        {
            for (j=0; j<n; j++)
            {
                // i : column index;
                // j : row index.
//              c[k*n+j]=Z[i*n+j];
                c[k][j]=Z[j*loc_sizes[0]+i];
                //			GlobalV::ofs_running << " Z=" << Z[i*n+j] << std::endl;
            }
            k++;
        }
    }
    MPI_Barrier(comm);

    for (i=1; i<nprocs; i++)
    {
        const int rows = loc_sizes[i];
        const int mpi_times = rows;
        if (myid==i)
        {
            for (j=0; j<mpi_times; j++)
            {
                int tag = j;

                double *send = new double[n];

                int count = 0;
                for (int m=0; m<rows*n; m++)
                {
                    if (m%rows==j)
                    {
                        send[count] = Z[m];
                        ++count;
                    }
                }

                MPI_Send(send,n,MPI_DOUBLE,0,tag,comm);

                delete[] send;
            }
        }
        else if (myid==0)
        {
            int col=0;
            for (j=0; j<mpi_times; j++)
            {
                double *ctmp = new double[n];
                ModuleBase::GlobalFunc::ZEROS(ctmp, n);
                int tag = j;
                MPI_Recv(ctmp,n,MPI_DOUBLE,i,tag,comm,&status);

                for (int m=0; m<n; m++)
                {
                    c[k][m]=ctmp[m];
//					GlobalV::ofs_running << " receive Z=" << ctmp[m] << std::endl;
                }
                k++;

                delete[] ctmp;
            }
        }
        MPI_Barrier(comm);
    }

    for (int i=0; i<GlobalV::NLOCAL; i++)
    {
        Parallel_Common::bcast_double(c[i],GlobalV::NLOCAL);
    }

    time_t time_end = time(NULL);
//    GlobalV::ofs_running << " End   gath_full_eig Time : " << ctime(&time_end);

    return;
}

void Pdiag_Double::gath_full_eig_complex(MPI_Comm comm,int n,std::complex<double> **c,std::complex<double> *Z)
{
    ModuleBase::TITLE("Pdiag_Double","gath_full_eig_complex");

    time_t time_start = time(NULL);
    //GlobalV::ofs_running << " Start gath_full_eig_complex Time : " << ctime(&time_start);

    int i,j,k,incx=1;
    int *loc_sizes,loc_size,nprocs,myid;
    MPI_Status status;
    MPI_Comm_size(comm,&nprocs);
    MPI_Comm_rank(comm,&myid);
    loc_sizes =(int*)malloc(sizeof(int)*nprocs);
    loc_size=n/nprocs;
    for (i=0; i<nprocs; i++)
    {
        if (i<n%nprocs)
        {
            loc_sizes[i]=loc_size+1;
        }
        else
        {
            loc_sizes[i]=loc_size;
        }
    }
    if (myid==0)
    {
        k=0;
        for (i=0; i<loc_sizes[0]; i++)
        {
            for (j=0; j<n; j++)
            {
                // i : column index;
                // j : row index.
//              c[k*n+j]=Z[i*n+j];
                c[k][j]=Z[j*loc_sizes[0]+i];
                //			GlobalV::ofs_running << " Z=" << Z[i*n+j] << std::endl;
            }
            k++;
        }
    }
    MPI_Barrier(comm);

	for (i=1; i<nprocs; i++)
	{
		const int rows = loc_sizes[i];
		const int mpi_times = rows;
		if (myid==i)
		{
			for (j=0; j<mpi_times; j++)
			{
				int tag = j;

				std::complex<double> *send = new std::complex<double>[n];

				int count = 0;
				for (int m=0; m<rows*n; m++)
				{
					if (m%rows==j)
					{
						send[count] = Z[m];
						++count;
					}
				}

				MPI_Send(send,n,mpicomplex,0,tag,comm);

				delete[] send;
			}
		}
		else if (myid==0)
		{
			int col=0;
			for (j=0; j<mpi_times; j++)
			{
				std::complex<double> *ctmp = new std::complex<double>[n];
				ModuleBase::GlobalFunc::ZEROS(ctmp, n);
				int tag = j;
				MPI_Recv(ctmp,n,mpicomplex,i,tag,comm,&status);

				for (int m=0; m<n; m++)
				{
					c[k][m]=ctmp[m];
					//					GlobalV::ofs_running << " receive Z=" << ctmp[m] << std::endl;
				}
				k++;

				delete[] ctmp;
			}
		}
		MPI_Barrier(comm);
	}

	for (int i=0; i<GlobalV::NLOCAL; i++)
	{
		Parallel_Common::bcast_complex_double(c[i],GlobalV::NLOCAL);
	}

    time_t time_end = time(NULL);
    //GlobalV::ofs_running << " End   gath_full_eig Time : " << ctime(&time_end);
	ModuleBase::GlobalFunc::OUT_TIME("gather full eigenvalues",time_start,time_end);

    return;
}
#endif

//LiuXh add 2021-09-06, clear memory, totwfc and WFC_GAMMA_aug not used now
#ifdef __MPI
#include "../src_pw/occupy.h"
void Pdiag_Double::gath_eig(MPI_Comm comm,int n,double *Z)
{
    ModuleBase::TITLE("Pdiag_Double","gath_eig");
    time_t time_start = time(NULL);
//  GlobalV::ofs_running << " Start gath_eig Time : " << ctime(&time_start);

    const Parallel_Orbitals* pv = this->ParaV;
    
    int i, j, k;
    int nprocs,myid;
    MPI_Status status;
    MPI_Comm_size(comm,&nprocs);
    MPI_Comm_rank(comm,&myid);

    double **ctot;

	// mohan add 2010-07-03
	// the occupied bands are useless
	// for calculating charge density.
	if(GlobalV::DRANK > pv->lastband_in_proc)
	{
		delete[] Z;
	}

	// first we need to collect all
	// the occupied bands.
	// GlobalV::NBANDS * GlobalV::NLOCAL	
	if(myid==0)
	{
		ctot = new double*[GlobalV::NBANDS];
    	for (int i=0; i<GlobalV::NBANDS; i++)
    	{
        	ctot[i] = new double[GlobalV::NLOCAL];
        	ModuleBase::GlobalFunc::ZEROS(ctot[i], GlobalV::NLOCAL);
    	}
    	ModuleBase::Memory::record("Pdiag_Double","ctot",GlobalV::NBANDS*GlobalV::NLOCAL,"double");
	}

    k=0;
    if (myid==0)
    {
        // mohan add nbnd0 2010-07-02
        int nbnd0 = -1;
        if (GlobalV::NBANDS < pv->loc_sizes[0])
        {
			// means all bands in this processor
			// is needed ( is occupied)
            nbnd0 = GlobalV::NBANDS;
        }
        else
        {
			// means this processor only save
			// part of GlobalV::NBANDS.
            nbnd0 = pv->loc_sizes[0];
        }
        if(pv->testpb) GlobalV::ofs_running << " nbnd in processor 0 is " << nbnd0 << std::endl;

//printf("from 0 to %d\n",nbnd0-1);
        for (i=0; i<nbnd0; i++)
        {
            for (j=0; j<GlobalV::NLOCAL; j++)
            {
				// change the order in processor 0.
				// the contribution from processor 0.
                ctot[k][j]=Z[j*pv->loc_sizes[0]+i];
            }
            k++;
        }
		// Z is useless in processor 0 now.
		delete[] Z;
    }
    MPI_Barrier(comm);


    for (i=1; i<= pv->lastband_in_proc; i++)
    {
        // mohan fix bug 2010-07-02
        // rows indicates the data structure of Z.
        // mpi_times indicates the data distribution
        // time, each time send a band.
        int rows = pv->loc_sizes[i];
        int mpi_times;
        if (i==pv->lastband_in_proc)
        {
            mpi_times = pv->lastband_number;
        }
        else
        {
            mpi_times = pv->loc_sizes[i];
        }
        if(pv->testpb)GlobalV::ofs_running << " nbnd in processor " << i << " is " << mpi_times << std::endl;
        if (myid==i)
        {
            for (j=0; j<mpi_times; j++)
            {
                int tag = j;
                double *send = new double[n];
                int count = 0;

                for (int m=0; m<rows*n; m++)
                {
                    if (m%rows==j)
                    {
                        send[count] = Z[m];
                        ++count;
                    }
                }

				// send the data to processor 0.
                MPI_Send(send,n,MPI_DOUBLE,0,tag,comm);

                delete[] send;
            }
			// third part to delete Z;
			delete[] Z;
        }
        else if (myid==0)
        {
            int col=0;
            for (j=0; j<mpi_times; j++)
            {
                double *ctmp = new double[GlobalV::NLOCAL];
                ModuleBase::GlobalFunc::ZEROS(ctmp, GlobalV::NLOCAL);
                int tag = j;
                
				// Processor 0 receive the data from other processors.
				MPI_Recv(ctmp,n,MPI_DOUBLE,i,tag,comm,&status);

                for (int m=0; m<GlobalV::NLOCAL; m++)
                {
                    ctot[k][m]=ctmp[m];
//					GlobalV::ofs_running << " receive Z=" << ctmp[m] << std::endl;
                }
                k++;

                delete[] ctmp;
            }
        }
        //MPI_Barrier(comm);
    }
    if(pv->testpb)ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"Final k",k);
/*
if(myid==0){
	double *vect=new double[GlobalV::NLOCAL*GlobalV::NBANDS];
	double *form=new double[GlobalV::NBANDS*GlobalV::NBANDS];
	int x,y;
	for(x=0;x<GlobalV::NBANDS;x++){
		for(y=0;y<GlobalV::NLOCAL;y++){
			vect[x*GlobalV::NLOCAL+y]=ctot[x][y];
		}
	}
		char chT='T';
		char chN='N';
		int	ne = GlobalV::NBANDS;
		int m1, n1, k1;
		double ONE=1.0,ZERO=0.0;
		m1 = ne;
		n1 = ne;
		k1 = GlobalV::NLOCAL;
		dgemm_(&chT, &chN, &m1, &n1, &k1, &ONE, vect, &k1, vect, &k1, &ZERO, form, &m1);
		double di=0.0,oth=0.0;
		for(x=0;x<ne;x++){
			for(y=0;y<ne;y++){
				if(x==y){
					di+=fabs(form[x*ne+y]);
				}else{
					oth+=fabs(form[x*ne+y]);
				}
			}
		}
		di-=ne;
		printf("\n\n\ndi=%.16lf\n\n\nother=%.16lf\n\n\n",di,oth);
		//assert(0>1);
}
*/
MPI_Barrier(comm);

	// mohan add 2010-09-10
	// output the wave function if required.
	// this is a bad position to output wave functions.
	// but it works!
    if(this->out_lowf)
	{
		// read is in ../src_algorithms/wf_local.cpp
		std::stringstream ss;
		ss << GlobalV::global_out_dir << "LOWF_GAMMA_S" << GlobalV::CURRENT_SPIN+1 << ".dat";
		// mohan add 2012-04-03, because we need the occupations for the
		// first iteration. 
		Occupy::calculate_weights();
		WF_Local::write_lowf( ss.str(), ctot );//mohan add 2010-09-09        
	}

	// mohan add 2010-09-10
	// distribution of local wave functions 
	// to each processor.
	// only used for GlobalV::GAMMA_ONLY_LOCAL
	//WF_Local::distri_lowf( ctot, wfc);


	// clean staff.
	if(myid==0)
	{
    	for (int i=0; i<GlobalV::NBANDS; i++)
    	{
        	delete[] ctot[i];
    	}
    	delete[] ctot;
	}

    time_t time_end = time(NULL);
    //GlobalV::ofs_running << " End   gath_eig Time : " << ctime(&time_end);
	ModuleBase::GlobalFunc::OUT_TIME("gather eigenvalues",time_start,time_end);
    return;
}
#endif
