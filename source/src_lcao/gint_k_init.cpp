#include "gint_k_init.h"
#include "../src_pw/global.h"

Gint_k_init::Gint_k_init()
{ 
	nbx = nby = nbz =0;
	ncxyz = 0;
	nbz_start = 0;
	nnn = new int[1];
}

Gint_k_init::~Gint_k_init()
{
	delete[] nnn;
}

void Gint_k_init::init(
	const int &nbx_in,
	const int &nby_in,
	const int &nbz_in,
	const int &nbz_start_in,
	const int &ncxyz_in
	)
{
	TITLE(GlobalV::ofs_running,"Gint_k_init","init");

	this->nbx = nbx_in;
	this->nby = nby_in;
	this->nbz = nbz_in;
	this->ncxyz = ncxyz_in;
	this->nbz_start = nbz_start_in;
	assert(nbx>0);
	assert(nby>0);
	assert(nbz>0);
	assert(ncxyz>0);

	// claculate the maximal orbital numbers 
	// for each type of atom.
	delete[] this->nnn;
	this->nnn = new int[ucell.ntype];
	for(int it=0; it<ucell.ntype; it++)
	{
		this->nnn[it] = (ucell.atoms[it].nwl+1) * (ucell.atoms[it].nwl+1);
	}

	assert( ucell.omega > 0.0);

	return;
}

