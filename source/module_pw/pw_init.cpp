#include "./pw_basis.h"

namespace ModulePW
{
//
//Init the grids for FFT
//Input: lattice vectors of the cell, Energy cut off for G^2/2
//Output: nx, ny, nz, nxyz, latvec, G, GT, GGT
//
void PW_Basis:: initgrids(
        ModuleBase::Matrix3 latvec_in, // Unitcell lattice vectors
        double gridecut
)
{
    //init latice
    this->latvec = latvec_in;
    this->GT = latvec.Inverse();
	this->G  = GT.Transpose();
	this->GGT = G * GT;
    

    //------------------------------------------------------------
    //-------------------------init grids-------------------------
    //------------------------------------------------------------
     ModuleBase::Vector3<double> lat;
    int *ibox = new int[3];// ibox[i] are the minimal FFT dimensions,
    
    lat.x = latvec.e11;
    lat.y = latvec.e12;
    lat.z = latvec.e13;
    ibox[0] = 2 * int(sqrt(gridecut) * sqrt(lat * lat)) + 1;

    lat.x = latvec.e21;
    lat.y = latvec.e22;
    lat.z = latvec.e23;
    ibox[1] = 2 * int(sqrt(gridecut) * sqrt(lat * lat)) + 1;

    lat.x = latvec.e31;
    lat.y = latvec.e32;
    lat.z = latvec.e33;
    ibox[2] = 2 * int(sqrt(gridecut) * sqrt(lat * lat)) + 1;
    //lat*lat=lat.x*lat.x+lat.y*lat.y+lat.z+lat.z

    // Find the minimal FFT box size the factors into the primes (2,3,5,7).
    for (int i = 0; i < 3; i++)
    {
    	int b = 0;
        int n2 = 0;
        int n3 = 0;
        int n5 = 0;
        //int n7 = 0;
        bool done_factoring = false;
	
		int ns = 0;
        // increase ibox[i] by 1 until it is totally factorizable by (2,3,5,7) 
        do
        {
            ibox[i] += 1;
			b = ibox[i];          

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
        }
        while (b != 1);
        //  b==1 means fftbox[i] is (2,3,5,7) factorizable 
    }
    this->nx = ibox[0];
    this->ny = ibox[1];
    this->nz = ibox[2];
    this->nxyz = this->nx * this->ny * this->nz;
    this->nxy = this->nx * this->ny;
    delete[] ibox;

    
    return;
}

//
//Init the grids for FFT
//Input: lattice vectors of the cell, nx, ny, nz
//Output: nx, ny, nz, nxyz, latvec, G, GT, GGT
//
void PW_Basis:: initgrids(
    ModuleBase::Matrix3 latvec_in, // Unitcell lattice vectors
    int nx_in, int ny_in, int nz_in
)
{
    this->latvec = latvec_in;
    this->GT = latvec.Inverse();
	this->G  = GT.Transpose();
	this->GGT = G * GT;
    this->nx = nx_in;
    this->ny = ny_in;
    this->nz = nz_in;
    this->nxyz = this->nx * this->ny * this->nz;
    return;
}


//Init some parameters
void PW_Basis:: initparameters(
    bool gamma_only_in,
    double ggecut_in,
    int poolnproc_in,
    int poolrank_in,
    int distribution_type_in
)
{
    this->gamma_only = gamma_only_in;
    this->ggecut = ggecut_in;
    this->poolnproc = poolnproc_in;
    this->poolrank = poolrank_in;
    this->distribution_type = distribution_type_in;
}

}