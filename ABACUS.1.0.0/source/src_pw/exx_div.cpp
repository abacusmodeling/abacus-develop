#include "exx_div.h"
#include "global.h"

Exx_Divergence::Exx_Divergence()
{
	this->factor = new double[1];	
}

Exx_Divergence::~Exx_Divergence()
{
	delete[] factor;	
}

void Exx_Divergence::init(void)
{
	// default : use gamma_extrapolation

	this->gamma_extra = true;

	if(this->gamma_extra)
	{
		this->grid_fac = 8.0 / 7.0;
	}
	else
	{
		this->grid_fac = 1.0;
	}	
	return;
}

// set the factor = 1/|G + k - q|^2
// and then deal with the divergence term.
void Exx_Divergence::get_factor(const int &ik, const int &iq)
{
	delete[] factor;
	this->factor = new double[wf.npw];

	static double e2 = 2.0; // in Rydberg unit.
	const double multi = e2 * FOUR_PI / ucell.tpiba2;

	if(pw.gstart==1)
	{
		factor[0] = 141.966625802957; // tmp
		//factor[0] = 0.0;
	}

	for(int ig=pw.gstart; ig<wf.npw; ig++)
	{
		this->factor[ig] = multi / ( pw.gcar[ig] + kv.kvec_c[ik] - kv.kvec_c[iq] ).norm2() * grid_fac;
	//	cout << " ig=" << ig << " factor=" << factor[ig] << endl;
	}
	//QUIT();
	

	return;
}
