#include "sto_wf.h"
#include "time.h"

//---------Temporary------------------------------------
#include "global.h"
#include "../module_base/global_function.h"
#include "../module_base/complexmatrix.h"
//------------------------------------------------------

Stochastic_WF::Stochastic_WF()
{
    chiortho  = NULL;
	chi0  = NULL;
    shchi = NULL;
}

Stochastic_WF::~Stochastic_WF()
{ 
    if(chi0 != NULL)        delete[] chi0;
    if(shchi != NULL)       delete[] shchi;
    if(chiortho != NULL)    delete[] chiortho;
    
}

void Stochastic_WF::init(const int nks)
{
    chi0 = new ModuleBase::ComplexMatrix[GlobalC::kv.nks];
    shchi = new ModuleBase::ComplexMatrix[GlobalC::kv.nks];
    chiortho = new ModuleBase::ComplexMatrix[GlobalC::kv.nks];
}

void Init_Sto_Orbitals(Stochastic_WF& sto_wf, const int seed_in)
{
    if(seed_in == 0 || seed_in == -1)
        srand((unsigned)time(NULL)+GlobalV::MY_RANK*10000);
    else
    {
        srand((unsigned)abs(seed_in)+GlobalV::MY_RANK*10000);
    }

    bool firstrankmore = false; 
    int igroup;
    //I am not sure which is better.
    //former processor calculate more bands
    if(firstrankmore)       igroup = GlobalV::MY_STOGROUP;
    //latter processor calculate more bands
    else                    igroup = GlobalV::NSTOGROUP - GlobalV::MY_STOGROUP - 1;
    
    const int nchi = sto_wf.nchi;
    const int ndim = GlobalC::wf.npwx;
    const int ngroup = GlobalV::NSTOGROUP;
    const int nks = GlobalC::kv.nks;

    int tmpnchip = int(nchi / ngroup);
    if(igroup < nchi%ngroup) ++tmpnchip;
    for(int ik = 0 ; ik < nks ; ++ik)
    {
        sto_wf.nchip[ik] = tmpnchip; 
        sto_wf.chi0[ik].create(tmpnchip,ndim,false);
        if(seed_in >= 0)
            for(int i = 0 ; i < ndim ; ++i)
            {
                const double phi = 2 * ModuleBase::PI * rand()/double(RAND_MAX);
                sto_wf.chi0[ik].c[i] = complex<double>(cos(phi), sin(phi)) / sqrt(double(nchi));
            }
        else
            for(int i = 0; i < ndim ; ++i)
            {
                if(rand()/double(RAND_MAX) < 0.5)
                    sto_wf.chi0[ik].c[i]=-1.0 / sqrt(double(nchi));
                else
                    sto_wf.chi0[ik].c[i]=1.0 / sqrt(double(nchi));
            }
    }
}

#ifdef __MPI
void Init_Com_Orbitals(Stochastic_WF& sto_wf, K_Vectors& kv)
{
    const bool firstrankmore = false; 
    int igroup;
    //former processor calculate more bands
    if(firstrankmore)       igroup = GlobalV::MY_STOGROUP;
    //latter processor calculate more bands
    else                    igroup = GlobalV::NSTOGROUP - GlobalV::MY_STOGROUP - 1;

    const int nks = kv.nks;
    const int npool = GlobalV::KPAR;
    const int ngroup = GlobalV::NSTOGROUP;
    const int n_in_pool = GlobalV::NPROC_IN_POOL;
    const int i_in_group = GlobalV::RANK_IN_STOGROUP;
    const int i_in_pool = GlobalV::RANK_IN_POOL;
    const int ndim = GlobalC::wf.npwx;

    int* totnpw = new int[nks];
    for(int ik = 0; ik < nks ; ++ik)
    {
        int* npwip = new int [n_in_pool];
        int* rec = new int [n_in_pool];
        int* displ = new int [n_in_pool];
        const int npw = kv.ngk[ik];

        for(int i_in_p = 0 ; i_in_p < n_in_pool ;++i_in_p)
        {
            rec[i_in_p] = 1;
            displ[i_in_p] = i_in_p;
        }
        MPI_Allgatherv(&npw, 1, MPI_INT, npwip, rec, displ, MPI_INT, POOL_WORLD);
        for(int ip = 0; ip < npool; ++ip)
        {
            totnpw[ik] += npwip[ip];
        }

        int tmpnchip = int(totnpw[ik]/ngroup);
        if(igroup < totnpw[ik] % ngroup)     ++tmpnchip;
        sto_wf.nchip[ik] = tmpnchip;
        sto_wf.chi0[ik].create(tmpnchip,ndim,true);

        const int re = totnpw[ik] % ngroup;
        int ip = 0, ig0 = 0;
        //give value to orbitals in one parallel group one by one.
        for(int ichi = 0 ; ichi < tmpnchip ; ++ichi)
        {
            int ig;
            if(igroup < re)
			{
                //It has more nchip.
                ig = igroup * tmpnchip + ichi - ig0;
			}
            else
			{
                //It has less nchip and should add re.
                ig = igroup * tmpnchip + re + ichi - ig0;
			}
            //get which ip stores this ig.
            while(ig >= npwip[ip])
            {
                ig -= npwip[ip];
                ig0 += npwip[ip];
                ++ip;
            }
            if(i_in_pool == ip)
			{
                sto_wf.chi0[ik](ichi , ig) = 1;
			}
        }
        
        delete[] npwip;
        delete[] rec;
        delete[] displ;
    }
    delete[] totnpw;

    // int iktot = Pkpoints.startk_pool[MY_POOL]+1;
    // if(MY_STOGROUP==0)
    // {
    //     for(int ip = 0 ; ip < NPOOL ; ++ip)
    //     {
    //         if(MY_POOL == ip)
    //         {   
    //             if(ip!=0 && RANK_IN_POOL==0) cout.clear();
    //             for(int ik = 0; ik < kv.nks ; ++ik)
    //             {
    //                 cout<<setw(7)<<"Kpoint "<<iktot<<" uses a compelete basis including "<<totnpw[ik]<<" plane waves"<<endl;
    //                 iktot++;
    //             }
    //             if(ip!=0 && RANK_IN_POOL==0) cout.setstate(ios::failbit);
    //         }
    //         MPI_Barrier(STO_WORLD);
    //     }
    // }
    
}
#else
void Init_Com_Orbitals(Stochastic_WF& sto_wf, K_Vectors& kv)
{
    const int ndim = GlobalC::wf.npwx;
    for(int ik = 0 ; ik < kv.nks ; ++ik)
    {
        
        chi0[ik].create(nchip[ik],ndim,true);
        for(int ichi = 0 ; ichi < kv.ngk[ik] ; ++ichi)
        {
            sto_wf.chi0[ik](ichi, ichi) = 1;
        }
    }
}
#endif

