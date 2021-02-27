#include "parallel_kpoints.h"
#include "parallel_global.h"

Parallel_Kpoints::Parallel_Kpoints()
{
    nks_pool = new int[1];
    startk_pool = new int[1];
    whichpool = new int[1];
}

Parallel_Kpoints::~Parallel_Kpoints()
{
    delete[] nks_pool;
    delete[] startk_pool;
    delete[] whichpool;
}

void Parallel_Kpoints::init_pools(void)
{
#ifdef __MPI
//----------------------------------------------------------
// CALL Function : divide_pools
//----------------------------------------------------------
    this->divide_pools();

// for test
// turn on when you want to check the index of pools.
/*
    if (MY_RANK==0)
    {
        cout << "\n     " << setw(8) << "MY_RANK"
             << setw(8) << "MY_POOL"
             << setw(13) << "RANK_IN_POOL"
             << setw(6) << "NPROC"
             << setw(6) << "NPOOL"
             << setw(14) << "NPROC_IN_POOL" << endl;
    }
    for (int i=0; i<NPROC; i++)
    {
        if (MY_RANK == i)
        {
            cout << " I'm " << setw(8) << MY_RANK
                 << setw(8) << MY_POOL
                 << setw(13) << RANK_IN_POOL
                 << setw(6) << NPROC
                 << setw(6) << NPOOL
                 << setw(14) << NPROC_IN_POOL << endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    if (MY_RANK != 0 )
    {
        cout.rdbuf(NULL);
    }
*/

    return;
#endif
}

#ifdef __MPI
void Parallel_Kpoints::divide_pools(void)
{
    //cout<<"\n ==> mpi_split()"<<endl;
    int i=0;
    int j=0;
    if (NPROC<NPOOL)
    {
        cout<<"\n NPROC=" << NPROC << " NPOOL=" << NPOOL;
        cout<<"Error : Too many pools !"<<endl;
        exit(0);
    }
    //if(kv.nkstot<NPOOL) cout<<"Error !"<<endl;

    // (1) per process in each pool
    NPROC_IN_POOL = NPROC/NPOOL;
    if (MY_RANK < (NPROC%NPOOL)*(NPROC_IN_POOL+1))
    {
        NPROC_IN_POOL++;
    }

    // (2) To know how many process in pool j.
    nproc_pool = new int[NPOOL];
    ZEROS(nproc_pool, NPOOL);
    for (i=0; i<NPROC; i++)
    {
        j = i%NPOOL;
        nproc_pool[j]++;
    }

    // (3) To know start proc index in each pool.
    startpro_pool = new int[NPOOL];
    ZEROS(startpro_pool, NPOOL);
    for (i=1; i<NPOOL; i++)
    {
        startpro_pool[i]=startpro_pool[i-1]+nproc_pool[i-1];
    }

    // use 'MY_RANK' to know 'MY_POOL'.
    for (i=0; i<NPOOL; i++)
    {
        if (MY_RANK >= startpro_pool[i])
        {
            MY_POOL=i;
        }
    }

    int key = 1;
    RANK_IN_POOL = MY_RANK-startpro_pool[MY_POOL];

    //========================================================
    // MPI_Comm_Split: Creates new communicators based on
    // colors(2nd parameter) and keys(3rd parameter)
    // Note: The color must be non-negative or MPI_UNDEFINED.
    //========================================================
    MPI_Comm_split(MPI_COMM_WORLD,MY_POOL,key,&POOL_WORLD);

    return;
}
#endif


// the kpoints here are reduced after symmetry applied.
void Parallel_Kpoints::kinfo(int &nkstot)
{
#ifdef __MPI
    Parallel_Common::bcast_int(nkstot);
    this->get_nks_pool(nkstot);
    this->get_startk_pool(nkstot);
    this->get_whichpool(nkstot);
#endif
    return;
}

#ifdef __MPI
void Parallel_Kpoints::get_whichpool(const int &nkstot)
{
    delete[] whichpool;
    this->whichpool = new int[nkstot];
    ZEROS(whichpool, nkstot);

	//cout << " calculate : whichpool" << endl;
	//cout << " nkstot is " << nkstot << endl;

    for (int i=0; i<NPOOL; i++)
    {
        for (int ik=0; ik< this->nks_pool[i]; ik++)
        {
            const int k_now = ik + startk_pool[i];
            this->whichpool[k_now] = i;
            //ofs_running << "\n whichpool[" << k_now <<"] = " << whichpool[k_now];
        }
    }

    return;
}

void Parallel_Kpoints::get_nks_pool(const int &nkstot)
{
    delete[] nks_pool;
    this->nks_pool = new int[NPOOL];
    ZEROS(nks_pool, NPOOL);

    const int nks_ave = nkstot/NPOOL;
    const int remain = nkstot%NPOOL;

    //ofs_running << "\n nkstot = " << nkstot;
    //ofs_running << "\n NPOOL = " << NPOOL;
    //ofs_running << "\n nks_ave = " << nks_ave;

    for (int i=0; i<NPOOL; i++)
    {
        this->nks_pool[i] = nks_ave;
        if (i<remain)
        {
            nks_pool[i]++;
        }
        //ofs_running << "\n nks_pool[i] = " << nks_pool[i];
    }
    return;
}

void Parallel_Kpoints::get_startk_pool(const int &nkstot)
{
    delete[] startk_pool;
    startk_pool = new int[NPOOL];
    //const int remain = nkstot%NPOOL;

    startk_pool[0] = 0;
    for (int i=1; i<NPOOL; i++)
    {
        startk_pool[i] = startk_pool[i-1] + nks_pool[i-1];
        //ofs_running << "\n startk_pool[i] = " << startk_pool[i];
    }
    return;
}
#endif


void Parallel_Kpoints::pool_collection(double &value, const double *wk, const int &ik)
{
#ifdef __MPI

    if (RANK_IN_POOL==0)
    {
        const int ik_now = ik - this->startk_pool[MY_POOL];
        //ofs_running << "\n\n ik=" << ik << " ik_now=" << ik_now;
        if (this->whichpool[ik] == MY_POOL)
        {
            if (MY_POOL > 0)
            {
                value = wk[ik_now];
                MPI_Send(&value, 1, MPI_DOUBLE, 0, ik, MPI_COMM_WORLD);
                //ofs_running << "\n send wk[" << ik << "]=" << value;
            }
            else
            {
                value = wk[ik_now];
                //ofs_running << "\n wk[" << ik << "]=" << value;
            }
        }
        else
        {
            if (MY_RANK==0)
            {
                MPI_Status ierror;
                const int iproc = this->startpro_pool[ this->whichpool[ik] ];
                MPI_Recv(&value, 1, MPI_DOUBLE, iproc, ik, MPI_COMM_WORLD,&ierror);
                //ofs_running << "\n receive wk[" << ik << "]=" << value << " from proc=" << iproc;
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
#else
    value = wk[ik];
#endif
    return;
}


void Parallel_Kpoints::pool_collection(double *valuea, double *valueb, const realArray &a, const realArray &b, const int &ik)
{
    const int dim2 = a.getBound2();
    const int dim3 = a.getBound3();
    const int dim4 = a.getBound4();
    assert( a.getBound2() == b.getBound2() );
    assert( a.getBound3() == b.getBound3() );
    assert( a.getBound4() == b.getBound4() );
    const int dim = dim2 * dim3 * dim4;
#ifdef __MPI
    const int ik_now = ik - this->startk_pool[MY_POOL];
    const int begin = ik_now * dim2 * dim3 * dim4;
    double* pa = &a.ptr[begin];
    double* pb = &b.ptr[begin];

    const int pool = this->whichpool[ik];

	ofs_running << "\n ik=" << ik;

    if (RANK_IN_POOL==0)
    {
        if (MY_POOL==0)
        {
            if (pool==0)
            {
                for (int i=0; i<dim ; i++)
                {
                    valuea[i] = *pa;
                    valueb[i] = *pb;
                    ++pa;
                    ++pb;
                }
            }
            else
            {
				ofs_running << " receive data.";
                MPI_Status ierror;
                MPI_Recv(valuea, dim, MPI_DOUBLE, this->startpro_pool[pool], ik*2+0, MPI_COMM_WORLD,&ierror);
                MPI_Recv(valueb, dim, MPI_DOUBLE, this->startpro_pool[pool], ik*2+1, MPI_COMM_WORLD,&ierror);
            }
        }
        else
        {
            if (MY_POOL == pool)
            {
				ofs_running << " send data.";
                MPI_Send(pa, dim, MPI_DOUBLE, 0, ik*2+0, MPI_COMM_WORLD);
                MPI_Send(pb, dim, MPI_DOUBLE, 0, ik*2+1, MPI_COMM_WORLD);
            }
        }
    }
	else
	{
		ofs_running << "\n do nothing.";
	}
    MPI_Barrier(MPI_COMM_WORLD);

    /*


    	if(this->whichpool[ik] == MY_POOL)
    	{
    		if(MY_POOL > 0 && RANK_IN_POOL == 0)
    		{
    			// data transfer ends.
    			MPI_Send(pa, dim, MPI_DOUBLE, 0, ik*2,   MPI_COMM_WORLD);
    			MPI_Send(pb, dim, MPI_DOUBLE, 0, ik*2+1, MPI_COMM_WORLD);
    		}
    		else if(MY_POOL == 0 && MY_RANK == 0)
    		{
    //			cout << "\n ik = " << ik << endl;
    			// data transfer begin.
    			for(int i=0; i<dim; i++)
    			{
    				valuea[i] = *pa;
    				valueb[i] = *pb;
    				++pa;
    				++pb;
    			}
    			// data transfer ends.
    		}
    	}
    	else
    	{
    		if(MY_RANK==0)
    		{
    			MPI_Status* ierror;
    			const int iproc = this->startpro_pool[ this->whichpool[ik] ];
    			MPI_Recv(valuea, dim, MPI_DOUBLE, iproc, ik*2,   MPI_COMM_WORLD,ierror);
    			MPI_Recv(valueb, dim, MPI_DOUBLE, iproc, ik*2+1, MPI_COMM_WORLD,ierror);
    		}
    	}
    	*/
#else
    // data transfer ends.
    const int begin = ik * dim2 * dim3 * dim4;
    double* pa = &a.ptr[begin];
    double* pb = &b.ptr[begin];
    for (int i=0; i<dim; i++)
    {
        valuea[i] = *pa;
        valueb[i] = *pb;
        ++pa;
        ++pb;
    }
    // data transfer ends.
#endif
    return;
}
