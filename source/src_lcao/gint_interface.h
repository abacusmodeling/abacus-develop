#ifndef GINT_INTERFACE
#define GINT_INTERFACE

#include "gint_tools.h"

class Gint_Interface
{
    public:

    // the unified interface to grid integration
	void cal_gint(Gint_inout *inout);

    // preparing FFT grid
    void prep_grid(
        const int &nbx_in,
        const int &nby_in,
        const int &nbz_in,
        const int &nbz_start_in,
        const int& ncxyz_in);

    protected:
    // variables related to FFT grid
 	int nbx;
	int nby;
	int nbz;
	int ncxyz;
	int nbz_start;

    // dimension: [GlobalC::LNNR.nnrg] 
    // save the < phi_0i | V | phi_Rj > in sparse H matrix.
    bool pvpR_alloc_flag;
    double** pvpR_reduced;
};

#endif