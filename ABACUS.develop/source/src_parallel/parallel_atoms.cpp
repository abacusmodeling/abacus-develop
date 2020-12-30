#include "parallel_atoms.h"
#include "../src_pw/global.h"

Parallel_Atoms::Parallel_Atoms()
{
	keep_this_atom = new bool[1];
	nat = 0;
}

Parallel_Atoms::~Parallel_Atoms()
{
	delete[] keep_this_atom;
}


void Parallel_Atoms::cut_atoms(void)
{
	TITLE("Parallel_Atoms","cut_atoms");
	OUT(ofs_running,"TotalAtomNumber",ucell.nat);
	OUT(ofs_running,"AtomDistribution",ATOM_DISTRIBUTION);

	delete[] keep_this_atom;
	keep_this_atom = new bool[ucell.nat];
	for(int i=0; i<ucell.nat; i++)
	{
		keep_this_atom[i] = false;
	}
	
//	if(ATOM_DISTRIBUTION==0) return;

	// the simplest dividsion of atoms.
	// the alternative method is using chenkai's Hilbert space method.
	int iat = 0;
	//add by chenkai
	int* nap = new int[NPROC];
	for(int in = 0; in < NPROC; in++)
	{
		nap[in] = 0;
	}

	for(int it=0; it<ucell.ntype; it++)
	{
		for(int ia=0; ia<ucell.atoms[it].na; ia++)
		{
			const int proc = iat % NPROC;
			nap[ proc ]++;
			ofs_running << " Atom" << iat << " is in processor: " << proc+1 << endl; 
			++iat;
		}
	}
	iat = 0;
	for(int ip = 0; ip < NPROC; ip++)
	{
		for(int j = 0; j <nap[ip]; j++)
		{
			if(ip == MY_RANK)
			{
				keep_this_atom[iat] = true;
			}
			iat++;
		}
	}
	assert(iat == ucell.nat);
		
	this->nat = 0;
	this->nlocal = 0;
	ofs_running << " Atom index in this processor." << endl;
	for(int iat=0; iat<ucell.nat; iat++)
	{
		if(keep_this_atom[iat]) 
		{
			ofs_running << " iat=" << iat+1 << endl;
			++nat;

			int it = ucell.iat2it[iat];
			nlocal += ucell.atoms[it].nw;
		}
	}

	OUT(ofs_running,"Parallel_Atoms::nat",nat);
	OUT(ofs_running,"Parallel_Atoms::nlocal",nlocal);
	
	return;
}

// the function is called in Parallel_Orbitals, as one part of ParaO.set_trace.
void Parallel_Atoms::set_trace(int *trace_loc_row, int *trace_loc_col, int &nrow, int &ncol)
{		
	int iat=0;
	int i=0; // the real trace
	int nc=NLOCAL; 
	int nr=0;

	for(int iw=0; iw<NLOCAL; iw++)
	{
		trace_loc_col[iw] = iw;
	}

	for(int it=0; it<ucell.ntype; it++)
	{
		int nw = ucell.atoms[it].nw;
		if(NSPIN==4) nw *= 2;//added by zhengdy-soc
		for(int ia=0; ia<ucell.atoms[it].na; ia++)
		{
			if( this->keep_this_atom[iat] )
			{
				const int start = ucell.itiaiw2iwt(it, ia, 0);
				for(int iw=0; iw<nw; iw++)
				{
					const int iw_all = iw+start;
					trace_loc_row[iw_all] = i;
					++i;
				}
				nr += nw;
			}
			else
			{
				// do nothing.
			}
			++iat;
		}
	}	
	nrow = nr;
	ncol = nc;


/*
	ofs_running << setw(8) << "i" << setw(8) << "j" << setw(8) << "trace_i" << setw(8) << "trace_j" << endl;
	for(int i=0; i<NLOCAL; i++)
	{
		for(int j=0; j<NLOCAL; j++)
		{
			if(trace_loc_row[i]>=0 && trace_loc_col[j]>=0)
			{
				ofs_running << setw(8) << i+1 << setw(8) << j+1 
				<< setw(8) << trace_loc_row[i] << setw(8) << trace_loc_col[j] << endl;
			}
		}
	}
	*/

	return;
}
