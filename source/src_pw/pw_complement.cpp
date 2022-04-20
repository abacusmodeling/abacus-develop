#include "pw_complement.h"
#include "../module_base/mymath.h"
#include "../module_base/memory.h"
#include "../module_base/timer.h"

int PW_complement::get_total_pw_number(
    const double& ggcut_start,
    const double& ggcut_end,
    const int& nx,
    const int& ny,
    const int& nz,
    const ModuleBase::Matrix3& GGT)
{
    if (GlobalV::test_pw) ModuleBase::TITLE("PW_complement","get_total_pw_number");
    if (ggcut_end<=0.0)
    {
        ModuleBase::WARNING_QUIT("PW_complement::get_total_pw_number","ggcut <= 0.0");
    }
    // First, figure out the number of G-vectors within the cutoff G2max.
    int ibox[3]={0,0,0};

    // mohan modify 2008-3-25 nx -> ncx ,ny-> ncy, nz->ncz
    ibox[0] = int(nx / 2) + 1;
    ibox[1] = int(ny / 2) + 1;
    ibox[2] = int(nz / 2) + 1;

    // first create the auxiliary arrays for the 1D G vectors
    ModuleBase::Vector3<double> f;
    int ngm = 0 ;
    for (int i = -ibox[0]; i <= ibox[0]; i++)
    {
        for (int j = -ibox[1]; j <= ibox[1]; j++)
        {
            for (int k = -ibox[2]; k <= ibox[2]; k++)
            {
                f.x = i;
                f.y = j;
                f.z = k;
                double g2 = f * (GGT * f);  //G2= |f|^2 in the unit of (2Pi/lat0)^2
                
				if (g2 < ggcut_end && g2 >= ggcut_start)
                {
                    ngm++;
                }
            }
        }
    }
    return ngm;
}

void PW_complement::get_total_pw(
    double* gg,
    ModuleBase::Vector3<double> *ig,
    const double& ggcut_start,
    const double& ggcut_end,
    const int& nx,
    const int& ny,
    const int& nz,
    const ModuleBase::Matrix3& GGT, // GGT = G*GT.
    int& ngm// number of total plane waves.
)
{
    if (GlobalV::test_pw) ModuleBase::TITLE("PW_complement","get_total_pw");
    ModuleBase::timer::tick("PW_complement","get_total_pw");
    if (ggcut_end<=0.0)
    {
        ModuleBase::WARNING_QUIT("PW_complement::get_total_pw","ggcut <= 0.0");
    }
    // First, figure out the number of G-vectors within the cutoff G2max.
    int ibox[3]={0,0,0};

    // mohan modify 2008-3-25 nx -> ncx ,ny-> ncy, nz->ncz
    ibox[0] = int(nx / 2) + 1;
    ibox[1] = int(ny / 2) + 1;
    ibox[2] = int(nz / 2) + 1;

    // first create the auxiliary arrays for the 1D G vectors
    ModuleBase::Vector3<double> f;
    int ng = 0;
    for (int i = -ibox[0]; i <= ibox[0]; i++)
    {
        for (int j = -ibox[1]; j <= ibox[1]; j++)
        {
            for (int k = -ibox[2]; k <= ibox[2]; k++)
            {
                f.x = i;
                f.y = j;
                f.z = k;
                double g2 = f * (GGT * f);  //G2= |f|^2 in the unit of (2Pi/lat0)^2
                if (g2 < ggcut_end && g2 >= ggcut_start)
                {
                    /*** g std::vector indices f=(i,j,k) ***/
                    ig[ng] = f ;
                    gg[ng] = g2;
                    //std::cout<<std::setw(12)<<f.x<<std::setw(12)<<f.y<<std::setw(12)<<f.z<<std::setw(12)<<g2<<std::endl;
                    ++ng;
                }
            }
        }
    }

    //std::cout << "\n ng = " << ng;
    ModuleBase::timer::tick("PW_complement","get_total_pw");
    return;
}

void PW_complement::get_FFT_dimension(
    const ModuleBase::Matrix3 &latvec,
    const double &ggcut,
    int &nx_tmp,
    int &ny_tmp,
    int &nz_tmp,
	const int &bx,
	const int &by,
	const int &bz)
{
    if (GlobalV::test_pw) ModuleBase::TITLE("PW_complement","get_FFT_dimension");
    // read in the FFT dimension (Nx, Ny, Nz) from INPUT,
    // if Nx*Ny*Nz >0,use the input
    // FFT grid, otherwise generate the FFT grid in the code.

    int i = 0;
    ModuleBase::Vector3<double> lat;
    int ibox[3]={0,0,0};

    // ibox[i] are the minimal FFT dimensions,
    lat.x = latvec.e11;
    lat.y = latvec.e12;
    lat.z = latvec.e13;
    ibox[0] = 2 * int(sqrt(ggcut) * sqrt(lat * lat)) + 1;

    lat.x = latvec.e21;
    lat.y = latvec.e22;
    lat.z = latvec.e23;
    ibox[1] = 2 * int(sqrt(ggcut) * sqrt(lat * lat)) + 1;

    lat.x = latvec.e31;
    lat.y = latvec.e32;
    lat.z = latvec.e33;
    ibox[2] = 2 * int(sqrt(ggcut) * sqrt(lat * lat)) + 1;
    //lat*lat=lat.x*lat.x+lat.y*lat.y+lat.z+lat.z

    // Find the minimal FFT box size the factors into the primes (2,3,5,7).
    for (i = 0; i < 3; i++)
    {
    	int b = 0;
        int n2 = 0;
        int n3 = 0;
        int n5 = 0;
        //int n7 = 0;
        bool done_factoring = false;
	
		// mohan add 2011-04-22
		int s;
		if(i==0) s=bx;
		else if(i==1) s=by;
		else if(i==2) s=bz;
		int ns = 0;

        // increase ibox[i] by 1 until it is totally factorizable by (2,3,5,7) 
        do
        {
			b = ibox[i];

			// mohan add 2011-04-22            
			if( ibox[i] % s != 0) 
			{
				b = -1; // meaning less
			}
			else
			{
				//n2 = n3 = n5 = n7 = 0;
				n2 = n3 = n5 = ns = 0;
				done_factoring = false;

				while (!done_factoring)
				{
					if (b % 2 == 0) 
					{
						n2++;
						b /= 2;
						continue;
					}
					if (b % 3 == 0) 
					{
						n3++;
						b /= 3;
						continue;
					}
					if (b % 5 == 0) 
					{
						n5++;
						b /= 5;
						continue;
					}
					//if (b%7==0) { n7++; b /= 7; continue; }
					done_factoring = true;
				}
			}//
            ibox[i] += 1;
        }
        while (b != 1);
        ibox[i] -= 1;
        //  b==1 means fftbox[i] is (2,3,5,7) factorizable 
    }

    // Nx, Ny, Nz are the FFT dimensions

    nx_tmp = ibox[0];
    ny_tmp = ibox[1];
    nz_tmp = ibox[2];

    int nxyz_tmp = 0;
    nxyz_tmp = nx_tmp * ny_tmp * nz_tmp;

    if (GlobalV::test_pw > 1)
    {
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"nx",nx_tmp);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"ny",ny_tmp);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"nz",nz_tmp);
        ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"nxyz",nxyz_tmp);
    }
    return;
}


//==========================================================
// MEMBER FUNCTION :
// NAME : PW_Basis::setup_GVectors
// Second part of the initialization : find out all G
// vectors that |G|^2<=G2max and std::map it into a one
// dimentional array G1d in the increase order of |G|^2.
// Next generate the indices between the 1D array and
// the 3D G-grid and the FFT grid.
// generate : gg_global, g_global, ig_global
//==========================================================
void PW_complement::setup_GVectors(
    const ModuleBase::Matrix3& G,
    const int &ngmc_g_in,
    double* gg,
    ModuleBase::Vector3<double>* ig,
    ModuleBase::Vector3<double>* g)
{
    if (GlobalV::test_pw) ModuleBase::TITLE("PW_complement","setup_GVectors");
    ModuleBase::timer::tick("PW_complement","setup_GVectors");

    int *ind = new int[ngmc_g_in];// auxiliary array for the 1d G std::vector index
    ModuleBase::GlobalFunc::ZEROS( ind, ngmc_g_in );
    ind[0] = 0;//ind[0]=0, meaning ind[] is not initialized.

    //----------------------------------------------------------
    // CALL GLOBAL FUNCTION :
    // sort the norm of the G vectors of G1d[] in the ascending
    // order
    //----------------------------------------------------------
    ModuleBase::heapsort(ngmc_g_in, gg, ind);

    /***************************************************************/
    // adjust G1d[NG] according to the new order.
    // establish the link between the 1d array and 3d grid,
    // i.e. give the position "i" in 1d array,
    // G1d[i] is the std::vector in 3d grid, G1d2[i] is its norm.
    /***************************************************************/

    ModuleBase::Vector3<double> *igsort = new ModuleBase::Vector3<double>[ngmc_g_in];
    for (int i=0;i<ngmc_g_in;i++)
    {
        igsort[i] = ig[ind[i]];
	//	std::cout << i << " " << ind[i] << " " << ig[ind[i]].x << " " << ig[ind[i]].y << " " << ig[ind[i]].z << std::endl;
    }

	/* remain to be done by someone. mohan note 2011-07-23
	ModuleBase::Vector3<double> change;
	for(int i=0; i<ngmc_g_in; ++i)
	{
		for(int j=i; j<ngmc_g_in; ++j)
		{
			if( igsort[j].norm() == igsort[i].norm() )
			{
				if( igsort[j].x < igsort[i].x )
				{
					change = igsort[j];
					igsort[j] = igsort[i];
					igsort[i] = change;
				}
				
				{

				}
			}
			*/



    for (int i=0; i<ngmc_g_in; i++)
    {
        ig[i] = igsort[i];
    }

    delete[] igsort;
    delete[] ind;

    //----------------------------------------------------------
    // EXPLAIN : store the G vectors in 1d array (Cartian
    // coordinate)
    //----------------------------------------------------------
    for (int i=0; i<ngmc_g_in; i++)
    {
        g[i] = ig[i] * G;
        // g std::vector (in 2*pi/lat0), E_k=g*g*(2*pi/lat0)^2
	//	std::cout << i << " " << ig[i].x << " " << ig[i].y << " " << ig[i].z << std::endl;
    }


    ModuleBase::timer::tick("PW_complement","setup_GVectors");
    return;
}


// #ifndef __MPI
void PW_complement::get_ngmw(const int &ngmc, const double& ggwfc2, const double* gg_global, int &ngmw)
{
    int ng = 0;
    for (int ig=0; ig<ngmc; ig++)
    {
        //mohan modify 2008-3-25
        if (gg_global[ng] <= ggwfc2)
        {
            ng++;
        }
    }
    ngmw = ng;
	ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"ngmw",ngmw);
    return;
}

void PW_complement::get_ig2fftw(const int &ngmw, const int &nx, const int &ny, const int &nz,
                                const ModuleBase::Vector3<double> *gvec_direct, int *ig2fftw)
{
    //=====================================================================
    // establish the mapping between 1d array and the FFT box.for wfc
    //=====================================================================
    // firt , get ngmw = number of plane waves for wave functions.
    if (GlobalV::test_pw) ModuleBase::TITLE("PW_complement","get_ig2fftw");

    for (int ig=0; ig<ngmw; ig++)
    {
        int n1 = int(gvec_direct[ig].x);
        int n2 = int(gvec_direct[ig].y);
        int n3 = int(gvec_direct[ig].z);
        //ig2fftw[i]= n1 + n2*nx + n3*nx*ny;
        if (n1 < 0) n1 = n1 + nx;
        if (n2 < 0) n2 = n2 + ny;
        if (n3 < 0) n3 = n3 + nz;
        ig2fftw[ig] = n1 * ny * nz + n2 * nz + n3;
    }
    return;
}


void PW_complement::get_ig2fftc(const int &ngmc, const int &ncx, const int &ncy, const int &ncz,
                                const ModuleBase::Vector3<double> *ig, int* ig1, int *ig2, int *ig3, int *ig2fftc)
{
    //=================================================================
    // set ig2fftc & ig2fftw with the correct fft correspondence
    //=================================================================
    if (GlobalV::test_pw) ModuleBase::TITLE("PW_Basis","get_ig2fftc");
    // see ggen.f90, set ig2fftc ith the correct fft correspondence Page 4/6

    ModuleBase::GlobalFunc::ZEROS(ig2fftc, ngmc);
    ModuleBase::GlobalFunc::ZEROS(ig1, ngmc);
    ModuleBase::GlobalFunc::ZEROS(ig2, ngmc);
    ModuleBase::GlobalFunc::ZEROS(ig3, ngmc);

    ModuleBase::Memory::record("PW_complement","ig2fftc",ngmc,"int");
    ModuleBase::Memory::record("PW_complement","ig1",ngmc,"int");
    ModuleBase::Memory::record("PW_complement","ig2",ngmc,"int");
    ModuleBase::Memory::record("PW_complement","ig3",ngmc,"int");

    for (int i = 0; i < ngmc;i++)
    {
        int n1 = int(ig[i].x);//ig --> f --> (i,j,k)
        int n2 = int(ig[i].y);
        int n3 = int(ig[i].z);
        ig1[i] = n1 + ncx;
        ig2[i] = n2 + ncy;
        ig3[i] = n3 + ncz;
        if (n1 < 0)n1 = n1 + ncx;
        if (n2 < 0)n2 = n2 + ncy;
        if (n3 < 0)n3 = n3 + ncz;
        ig2fftc[i]=n1 * ncy * ncz + n2 * ncz + n3;
    }
    return;
}
// #endif

//LiuXh add a new function here,
//20180515
void PW_complement::get_total_pw_after_vc(
        double* gg0,
        double* gg,
        ModuleBase::Vector3<double> *ig,
        const double& ggcut_start,
        const double& ggcut_end,
        const int& nx,
        const int& ny,
        const int& nz,
        const ModuleBase::Matrix3& GGT, // GGT = G*GT.
        const ModuleBase::Matrix3& GGT0,
        int& ngm// number of total plane waves.
)
{
    if (GlobalV::test_pw) ModuleBase::TITLE("PW_complement","get_total_pw");
    ModuleBase::timer::tick("PW_complement","get_total_pw");
    if (ggcut_end<=0.0)
    {
        ModuleBase::WARNING_QUIT("PW_complement::get_total_pw","ggcut <= 0.0");
    }

    int ibox[3]={0,0,0};

    ibox[0] = int(nx / 2) + 1;
    ibox[1] = int(ny / 2) + 1;
    ibox[2] = int(nz / 2) + 1;

    ModuleBase::Vector3<double> f;
    int ng = 0;
    //int ng2 = 0;
    for (int i = -ibox[0]; i <= ibox[0]; i++)
    {
        for (int j = -ibox[1]; j <= ibox[1]; j++)
        {
            for (int k = -ibox[2]; k <= ibox[2]; k++)
            {
                f.x = i;
                f.y = j;
                f.z = k;
                double g2 = f * (GGT0 * f);
                //double g22 = f * (GGT * f);
                if (g2 < ggcut_end && g2 >= ggcut_start)
                {
                    ig[ng] = f ;
                    gg0[ng] = g2;
                    g2 = f * (GGT * f);
                    gg[ng] = g2;
                    ++ng;
                }
            }
        }
    }
    ModuleBase::timer::tick("PW_complement","get_total_pw");
    return;
}
