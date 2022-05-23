//=========================================================
//REFACTOR : Peize Lin, 2021.06.28
//=========================================================
#include "gint_gamma.h"
#include "gint_tools.h"
#include "grid_technique.h"
#include "../module_orbital/ORB_read.h"
#include "../src_pw/global.h"
#include "../module_base/blas_connector.h"
#include "../module_base/memory.h"
#include "../module_base/timer.h"

#include "global_fp.h" // mohan add 2021-01-30

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef __MKL
#include <mkl_service.h>
#endif

extern "C"
{
    void Cblacs_gridinfo(int icontxt, int* nprow, int *npcol, int *myprow, int *mypcol);
    void Cblacs_pinfo(int *myid, int *nprocs);
    void Cblacs_pcoord(int icontxt, int pnum, int *prow, int *pcol);
}

void Gint_Gamma::gint_kernel_vlocal(
	const int na_grid,
	const int grid_index,
	const double delta_r,
	double* vldr3,
	const int LD_pool,
    const int lgd_now,
	double* pvpR_grid_in)
{

    int * block_iw, * block_index, * block_size;
    Gint_Tools::get_block_info(na_grid, grid_index, block_iw, block_index, block_size);
    //------------------------------------------------------
    // whether the atom-grid distance is larger than cutoff
    //------------------------------------------------------
    bool **cal_flag = Gint_Tools::get_cal_flag(na_grid, grid_index);

    //------------------------------------------------------------------
    // compute atomic basis phi(r) with both radial and angular parts
    //------------------------------------------------------------------
    Gint_Tools::Array_Pool<double> psir_ylm(GlobalC::pw.bxyz, LD_pool);
    Gint_Tools::cal_psir_ylm(
        na_grid, grid_index, delta_r,
        block_index, block_size, 
        cal_flag,
        psir_ylm.ptr_2D);

    const Gint_Tools::Array_Pool<double> psir_vlbr3 = Gint_Tools::get_psir_vlbr3(
        na_grid, LD_pool, block_index, cal_flag, vldr3, psir_ylm.ptr_2D);

    this->cal_meshball_vlocal(
        na_grid, LD_pool, block_iw, block_size, block_index, cal_flag,
        vldr3, psir_ylm.ptr_2D, psir_vlbr3.ptr_2D, lgd_now, pvpR_grid_in);
        
    delete[] block_iw;
    delete[] block_index;
    delete[] block_size;
    for(int ib=0; ib<GlobalC::pw.bxyz; ++ib)
    {
        delete[] cal_flag[ib];
    }
    delete[] cal_flag;
}

void Gint_Gamma::cal_meshball_vlocal(
	const int na_grid,  					    // how many atoms on this (i,j,k) grid
	const int LD_pool,
	const int*const block_iw,				    // block_iw[na_grid],	index of wave functions for each block
	const int*const block_size, 			    // block_size[na_grid],	number of columns of a band
	const int*const block_index,		    	// block_index[na_grid+1], count total number of atomis orbitals
	const bool*const*const cal_flag,	    	// cal_flag[GlobalC::pw.bxyz][na_grid],	whether the atom-grid distance is larger than cutoff
	const double*const vldr3,			    	// vldr3[GlobalC::pw.bxyz]
	const double*const*const psir_ylm,		    // psir_ylm[GlobalC::pw.bxyz][LD_pool]
	const double*const*const psir_vlbr3,	    // psir_vlbr3[GlobalC::pw.bxyz][LD_pool]
	const int lgd_now,
	double*const GridVlocal) const	    // GridVlocal[lgd_now][lgd_now]
{
	const char transa='N', transb='T';
	const double alpha=1, beta=1;

	for(int ia1=0; ia1<na_grid; ++ia1)
	{
		const int iw1_lo=block_iw[ia1];
		const int m=block_size[ia1];
		for(int ia2=0; ia2<na_grid; ++ia2)
		{
			const int iw2_lo=block_iw[ia2];
			if(iw1_lo<=iw2_lo)
			{
                int first_ib=0;
                for(int ib=0; ib<GlobalC::pw.bxyz; ++ib)
                {
                    if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
                    {
                        first_ib=ib;
                        break;
                    }
                }
                int last_ib=0;
                for(int ib=GlobalC::pw.bxyz-1; ib>=0; --ib)
                {
                    if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
                    {
                        last_ib=ib+1;
                        break;
                    }
                }
                const int ib_length = last_ib-first_ib;
                if(ib_length<=0) continue;

                int cal_pair_num=0;
                for(int ib=first_ib; ib<last_ib; ++ib)
                {
                    cal_pair_num += cal_flag[ib][ia1] && cal_flag[ib][ia2];
                }

                const int n=block_size[ia2];
                if(cal_pair_num>ib_length/4)
                {
                    dgemm_(&transa, &transb, &n, &m, &ib_length, &alpha,
                        &psir_vlbr3[first_ib][block_index[ia2]], &LD_pool,
                        &psir_ylm[first_ib][block_index[ia1]], &LD_pool,
                        &beta, &GridVlocal[iw1_lo*lgd_now+iw2_lo], &lgd_now);   
                }
                else
                {
                    for(int ib=first_ib; ib<last_ib; ++ib)
                    {
                        if(cal_flag[ib][ia1] && cal_flag[ib][ia2])
                        {
                            int k=1;                            
                            dgemm_(&transa, &transb, &n, &m, &k, &alpha,
                                &psir_vlbr3[ib][block_index[ia2]], &LD_pool,
                                &psir_ylm[ib][block_index[ia1]], &LD_pool,
                                &beta, &GridVlocal[iw1_lo*lgd_now+iw2_lo], &lgd_now);                          
                        }
                    }
                }
			}
		}
	}
}

#ifdef __MPI
//------------------------------------------------------------------
// mohan add notes: 2021-03-11
// this subroutine is used to transform data from grid integrals
// to 2D-block distribution
// s stands for 'sender' and r stands for 'receiver'
//------------------------------------------------------------------
inline int setBufferParameter(
	MPI_Comm comm_2D,
	int blacs_ctxt,
	int nblk,
	int& s_index_siz,
	int*& s_local_index,
	int*& s_siz_pro,
	int*& s_dis_pro,
	int& s_siz,
	double*& s_buffer,
	int& r_index_siz,
	int*& r_global_index,
	int*& r_siz_pro,
	int*& r_dis_pro,
	int& r_siz,
	double*& r_buffer)
{
	//-----------------------------------------
    // setup blacs parameters
	//-----------------------------------------
    int nprows, npcols, nprocs;
    int myprow, mypcol, myproc;

    Cblacs_gridinfo(blacs_ctxt, &nprows, &npcols, &myprow, &mypcol);

	//-----------------------------------------
	// set index of current proor: myproc
	// set number of total proors: nprocs
	//-----------------------------------------
    Cblacs_pinfo(&myproc, &nprocs);

    // initialize data arrays
    delete[] s_siz_pro;
    delete[] s_dis_pro;
    delete[] r_siz_pro;
    delete[] r_dis_pro;

    s_siz_pro=new int[nprocs];
    s_dis_pro=new int[nprocs];
    r_siz_pro=new int[nprocs];
    r_dis_pro=new int[nprocs];

	//---------------------------------------------------------------------
    // build the local index to be sent to other pro (s_local_index),
    // the global index to be received from other pro (r_global_index),
    // the send/receive siz/dis for data exchange by MPI_Alltoall
	//---------------------------------------------------------------------
    s_index_siz=GlobalC::GridT.lgd*GlobalC::GridT.lgd*2;

    delete[] s_local_index;
    s_local_index=new int[s_index_siz];

    int *s_global_index=new int[s_index_siz];

    int pos=0;
    s_siz_pro[0]=0;
    for(int iproc=0; iproc<nprocs; ++iproc)
    {
        s_dis_pro[iproc]=pos;

        int iprow=0;
		int ipcol=0;
        Cblacs_pcoord(blacs_ctxt, iproc, &iprow, &ipcol);

        // find out the global index and local index of elements
		// in each pro based on 2D block cyclic distribution
        for(int irow=0, grow=0; grow<GlobalV::NLOCAL; ++irow)
        {
            grow=Local_Orbital_wfc::globalIndex(irow, nblk, nprows, iprow);
            if (grow >= GlobalV::NLOCAL)
                continue;
            int lrow = GlobalC::GridT.trace_lo[grow];
            if (lrow < 0)
                continue;

            for(int icol=0, gcol=0; gcol<GlobalV::NLOCAL; ++icol)
            {
                gcol=Local_Orbital_wfc::globalIndex(icol,nblk, npcols, ipcol);
                if (gcol >= GlobalV::NLOCAL)
                    continue;
                int lcol = GlobalC::GridT.trace_lo[gcol];
                if (lcol < 0)
                    continue;
                s_global_index[pos]=grow;
                s_global_index[pos+1]=gcol;
                s_local_index[pos]=lrow;
                s_local_index[pos+1]=lcol;
                pos+=2;
            }
        }
        s_siz_pro[iproc]=pos-s_dis_pro[iproc];
    }

    MPI_Alltoall(s_siz_pro, 1, MPI_INT,
                 r_siz_pro, 1, MPI_INT, comm_2D);

    r_index_siz=r_siz_pro[0];
    r_dis_pro[0]=0;
    for(int i=1; i<nprocs; ++i)
    {
        r_index_siz+=r_siz_pro[i];
        r_dis_pro[i]=r_dis_pro[i-1]+r_siz_pro[i-1];
    }

	delete[] r_global_index;
	r_global_index=new int[r_index_siz];

    // send the global index in sendBuffer to recvBuffer
    MPI_Alltoallv(s_global_index, s_siz_pro, s_dis_pro, MPI_INT,
                  r_global_index, r_siz_pro, r_dis_pro, MPI_INT, comm_2D);

    delete [] s_global_index;

    // the s_siz_pro, s_dis_pro, r_siz_pro,
    // and r_dis_pro will be used in transfer s_buffer, which
    // is half siz of s_global_index
    // we have to rebuild the siz and dis for each pro
    for (int iproc=0; iproc < nprocs; ++iproc)
    {
        s_siz_pro[iproc]=s_siz_pro[iproc]/2;
        s_dis_pro[iproc]=s_dis_pro[iproc]/2;
        r_siz_pro[iproc]=r_siz_pro[iproc]/2;
        r_dis_pro[iproc]=r_dis_pro[iproc]/2;
    }

    s_siz=s_index_siz/2;
	delete[] s_buffer;
	s_buffer=new double[s_siz];

    r_siz=r_index_siz/2;
	delete[] r_buffer;
	r_buffer=new double[r_siz];

    return 0;
}
#endif

void Gint_Gamma::vl_grid_to_2D(const int lgd_now, LCAO_Matrix &lm)
{
    // setup send buffer and receive buffer size
    // OUT(GlobalV::ofs_running, "Start transforming vlocal from grid distribute to 2D block");
    if(GlobalC::CHR.get_new_e_iteration())
    {
        ModuleBase::timer::tick("Gint_Gamma","distri_vl_index");
        #ifdef __MPI
        setBufferParameter(lm.ParaV->comm_2D, lm.ParaV->blacs_ctxt, lm.ParaV->nb,
                           this->sender_index_size, this->sender_local_index,
                           this->sender_size_process, this->sender_displacement_process,
                           this->sender_size, this->sender_buffer,
                           this->receiver_index_size, this->receiver_global_index,
                           this->receiver_size_process, this->receiver_displacement_process,
                           this->receiver_size, this->receiver_buffer);
        #endif
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "vlocal exchange index is built");
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "buffer size(M):", (this->sender_size+this->receiver_size)*sizeof(double)/1024/1024);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "buffer index size(M):", (this->sender_index_size+this->receiver_index_size)*sizeof(int)/1024/1024);
        ModuleBase::timer::tick("Gint_Gamma","distri_vl_index");
    }

    ModuleBase::timer::tick("Gint_Gamma","distri_vl_value");

    // put data to send buffer
    for(int i=0; i<this->sender_index_size; i+=2)
    {
        const int irow=this->sender_local_index[i];
        const int icol=this->sender_local_index[i+1];
        if(irow<=icol)
		{
            this->sender_buffer[i/2]=pvpR_grid[irow*lgd_now+icol];
		}
        else
		{
            this->sender_buffer[i/2]=pvpR_grid[icol*lgd_now+irow];
		}
    }
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "vlocal data are put in sender_buffer, size(M):", this->sender_size*8/1024/1024);

    // use mpi_alltoall to get local data
    #ifdef __MPI
    MPI_Alltoallv(this->sender_buffer, this->sender_size_process, this->sender_displacement_process, MPI_DOUBLE,
                  this->receiver_buffer, this->receiver_size_process,
					this->receiver_displacement_process, MPI_DOUBLE, lm.ParaV->comm_2D);
    #endif
    ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running, "vlocal data are exchanged, received size(M):", this->receiver_size*8/1024/1024);

    // put local data to H matrix
    for(int i=0; i<this->receiver_index_size; i+=2)
    {
        const int g_row=this->receiver_global_index[i];
        const int g_col=this->receiver_global_index[i+1];
        // if(g_col<0 || g_col>=GlobalV::NLOCAL||g_row<0 || g_row>=GlobalV::NLOCAL)
        // {
        //     OUT(GlobalV::ofs_running, "index error, i:", i);
        //     OUT(GlobalV::ofs_running, "index：", this->receiver_global_index[i]);
        //     OUT(GlobalV::ofs_running, "g_col:", g_col);
        //     OUT(GlobalV::ofs_running, "g_col:", g_col);
        // }
        lm.set_HSgamma(g_row,g_col,this->receiver_buffer[i/2],'L', lm.Hloc.data());
    }

    ModuleBase::timer::tick("Gint_Gamma","distri_vl_value");
    ModuleBase::timer::tick("Gint_Gamma","distri_vl");
}
