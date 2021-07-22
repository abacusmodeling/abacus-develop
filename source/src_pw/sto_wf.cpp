#include "sto_wf.h"
#include "global.h"

Stochastic_WF::Stochastic_WF()
{
    chiortho  = new ComplexMatrix[1];
	chi0  = new ComplexMatrix[1];
	emax_sto = 1;
	emin_sto = -1;
	stotype = "pw";
}

Stochastic_WF::~Stochastic_WF()
{ 
    delete[] chiortho;
    delete[] chi0;
}

void Stochastic_WF::init(void)
{
    //wait for init

    int ndim=0;
    if(stotype == "pw")
    {
        ndim = kv.ngk[0]; // only GAMMA point temporarily
    }
    else
    {
        ndim = pw.nrxx;
    }
    int totnpw=0;
#ifdef __MPI
    int * npwip = new int [NPROC_IN_POOL];
    int *rec = new int [NPROC_IN_POOL];
    int *displ = new int [NPROC_IN_POOL];
    for(int ip = 0 ; ip < NPROC_IN_POOL ;++ip)
    {
        rec[ip] = 1;
        displ[ip] = ip;
    }
    MPI_Allgatherv(&kv.ngk[0], 1, MPI_INT, npwip, rec, displ, MPI_INT, POOL_WORLD);
    for(int ip = 0; ip < NPROC_IN_POOL; ++ip)
    {
        totnpw += npwip[ip];
    }
#else
    totnpw = ndim;
#endif
    
    
    //distribute nchi for each process
    bool allbase = false;
    if(nchi == 0)
    {
        nchi = totnpw-NBANDS;
        cout<<"Using all normal bases: "<<totnpw<<endl;
        allbase = true;
    }
    nchip = int(nchi/NPOOL);
    if(NPOOL - MY_POOL - 1 < nchi%NPOOL) ++nchip;

    complex<double> ui(0,1);

    //We temporarily init one group of orbitals for all k points.
    delete[] chi0;
    chi0 = new ComplexMatrix[1]; 
    
    
    //srand((unsigned)time(NULL)+MY_RANK*10000);
    srand((unsigned)MY_RANK*10000);
    //srand((unsigned)0);
    
    if(allbase)
    {
        chi0[0].create(nchip,ndim,true);
        int re = NPOOL - nchi % NPOOL;
        int ip = 0, ig0 = 0;
        int ig;
        for(int i = 0 ; i < nchip ; ++i)
        {
#ifdef __MPI
            if(MY_POOL < re)
			{
                ig = MY_POOL * nchip + i - ig0;
			}
            else
			{
                ig = MY_POOL * nchip - re + i - ig0;
			}
            while(ig >= npwip[ip])
            {
                ig -= npwip[ip];
                ig0 += npwip[ip];
                ip++;
            }
            if(RANK_IN_POOL == ip)
			{
                chi0[0](i , ig) = 1;
			}
#else
                chi0[0](i , ig) = 1;
#endif
        }
    }
    else
    {
        //init with random number
        chi0[0].create(nchip,ndim,false);
        for(int i=0; i<chi0[0].size; ++i)
        {
            chi0[0].c[i]=exp(2*PI*rand()/double(RAND_MAX)*ui) / sqrt(double(nchi));
        }
    }
    

    //

    

    delete[] chiortho;
#ifdef __MPI
    delete[] npwip;
    delete[] rec;
    delete[] displ;
#endif
    int nkk = 1; // We temporarily use gamma k point.
    chiortho = new ComplexMatrix[1];
    if(NBANDS > 0)
    {
        chiortho[0].create(nchip,ndim,false);
    }
    
    return;
}


