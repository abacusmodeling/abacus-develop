#ifndef PW_BASIS_BIG_H
#define PW_BASIS_BIG_H
#include "../module_base/constants.h"
#include "../module_base/global_function.h"

// temporary class, because previous ABACUS consider big grid for fft grids 
// which are used for grid integration in LCAO.
// In fact, it is unnecessary. It will be moved after grid integration is refactored.
namespace ModulePW
{

class PW_Basis_Big: public PW_Basis
{
public:
    
    // combine [bx,by,bz] FFT grids into a big one
	// typical values are bx=2, by=2, bz=2
	// nbx=nx/bx, nby=ny/by, nbz=nz/bz, 
    PW_Basis_Big(){
        bx = 1;
        by = 1;
        bz = 1;
    }
    ~PW_Basis_Big(){};
    void setbxyz(const int bx_in, const int by_in, const int bz_in)
    {
        bx = bx_in;
        by = by_in;
        bz = bz_in;
    }
    int bx,by,bz;
    int nbx, nby, nbz;
    virtual void initgrids(const double lat0_in,const ModuleBase::Matrix3 latvec_in,
        const double gridecut, const int poolnproc_in, const int poolrank_in){
        //init lattice
    this->lat0 = lat0_in;
    this->latvec = latvec_in;
    this->GT = latvec.Inverse();
	this->G  = GT.Transpose();
	this->GGT = G * GT;
    this->poolnproc = poolnproc_in;
    this->poolrank = poolrank_in;
    

    //------------------------------------------------------------
    //-------------------------init grids-------------------------
    //------------------------------------------------------------
    double tpiba2 = ModuleBase::TWO_PI * ModuleBase::TWO_PI / this->lat0 / this->lat0;
    const double gridecut_lat = gridecut / tpiba2;
    ModuleBase::Vector3<double> lat;
    int *ibox = new int[3];
    
    lat.x = latvec.e11;
    lat.y = latvec.e12;
    lat.z = latvec.e13;
    ibox[0] = 2 * int(sqrt(gridecut_lat) * sqrt(lat * lat)) + 1;

    lat.x = latvec.e21;
    lat.y = latvec.e22;
    lat.z = latvec.e23;
    ibox[1] = 2 * int(sqrt(gridecut_lat) * sqrt(lat * lat)) + 1;

    lat.x = latvec.e31;
    lat.y = latvec.e32;
    lat.z = latvec.e33;
    ibox[2] = 2 * int(sqrt(gridecut_lat) * sqrt(lat * lat)) + 1;

    // Find the minimal FFT box size the factors into the primes (2,3,5,7).
    for (int i = 0; i < 3; i++)
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
			n2 = n3 = n5 = 0;
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
            ibox[i] += 1;
        }
        while (b != 1);
        ibox[i] -= 1;
        //  b==1 means fftbox[i] is (2,3,5,7) factorizable 
    }
    this->nx = ibox[0];
    this->ny = ibox[1];
    this->nz = ibox[2];
    this->nxy =this->nx * this->ny;
    this->nxyz = this->nxy * this->nz;
    this->nbx = this->nx / bx;
    this->nby = this->ny / by;
    this->nbz = this->nz / bz;

    delete[] ibox;    
    return;

    }
    virtual void distribute_r()
    {
       if(this->numz!=nullptr) delete[] this->numz; this->numz = new int[this->poolnproc];
       if(this->startz!=nullptr) delete[] this->startz; this->startz = new int[this->poolnproc];
        ModuleBase::GlobalFunc::ZEROS(this->numz, this->poolnproc);
        ModuleBase::GlobalFunc::ZEROS(this->startz, this->poolnproc);

        int npbz = this->nbz / this->poolnproc;
        int modbz = this->nbz % this->poolnproc;
        this->startz[0] = 0;
        for(int ip = 0 ; ip < this->poolnproc ; ++ip)
        {
            this->numz[ip] = npbz*this->bz;
            if(ip < modbz)   this->numz[ip]+=this->bz;
            if(ip < this->poolnproc - 1)   this->startz[ip+1] = this->startz[ip] + numz[ip];
            if(ip == this->poolrank) this->nplane = numz[ip];
        }
        this->nrxx = this->numz[this->poolrank] * this->nxy;
        return;
    }

};
}
#endif