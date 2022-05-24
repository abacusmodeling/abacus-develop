#include "gint_gamma.h"
#include "../src_pw/global.h"
#include "../module_base/ylm.h"
#include "../module_neighbor/sltk_atom_arrange.h"
#include "../module_base/timer.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef __MKL
#include <mkl_service.h>
#endif

Gint_Gamma::Gint_Gamma()
{
   
    sender_index_size = 1;
	sender_local_index = new int[1];
    sender_size_process = new int[1];
    sender_displacement_process = new int[1];
    sender_size=1;
    sender_buffer=new double[1];

    receiver_index_size=1;
    receiver_global_index = new int[1];
    receiver_size_process = new int[1];
    receiver_displacement_process = new int[1];
    receiver_size=1;
    receiver_buffer=new double[1];
}

Gint_Gamma::~Gint_Gamma()
{

    delete[] sender_local_index;
    delete[] sender_size_process;
    delete[] sender_displacement_process;
    delete[] sender_buffer;

    delete[] receiver_global_index;
    delete[] receiver_size_process;
    delete[] receiver_displacement_process;
    delete[] receiver_buffer;
}

void Gint_Gamma::cal_vlocal(Gint_inout *inout)
{
	this->cal_gint(inout);
	const int max_size = GlobalC::GridT.max_atom;
	const int lgd = GlobalC::GridT.lgd;

	if(inout->job==Gint_Tools::job_type::vlocal)
	{
		vl_grid_to_2D(lgd,inout->lm[0]);
		if (max_size >0 && lgd > 0) delete[] pvpR_grid;
	}
}