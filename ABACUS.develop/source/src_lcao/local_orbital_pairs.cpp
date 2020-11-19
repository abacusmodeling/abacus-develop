#include "local_orbital_pairs.h"
#include "lcao_orbitals.h"

Local_Orbital_Pairs::Local_Orbital_Pairs()
{
	// xinguo: lmax of orbital + 1
	lmax_aux = 3;	
	naux = 0;
}
Local_Orbital_Pairs::~Local_Orbital_Pairs(){}

void Local_Orbital_Pairs::init()
{
	TITLE("Local_Orbital_Pairs","init");

	// calculate the number of auxiliary basis.
	this->naux = 0;
	for(int it=0; it<ORB.get_ntype(); it++)
	{
		for(int L0=0; L0<lmax_aux+1; L0++)
		{
			for(int L1=0; L1<ORB.Phi[it].getLmax()+1; L1++)
			{
				for(int L2=0; L2<ORB.Phi[it].getLmax()+1; L2++)
				{
					if( abs(L1-L2)<=L0 && L1+L2>=L0)
					{
						for(int N1=0; N1<ORB.Phi[it].getNchi(L1); N1++)
						{
							for(int N2=0; N2<ORB.Phi[it].getNchi(L2); N2++)
							{
								cout << "\n L0=" << L0 << " L1=" << L1 << " L2=" << L2 
								<< " N1=" << N1 << " N2=" << N2;
								++naux;
							}
						}
					}
				}
			}
		}
	}	
	cout << "\n Number of auxiliary basis = " << naux << endl;

	this->aux_basis = new double*[naux];
	for(int i=0; i<naux; i++)
	{
		aux_basis[i] = new double[mesh];
	}
	
	for(int it=0; it<ORB.get_ntype(); it++)
	{
		for(int L0=0; L0<lmax_aux+1; L0++)
		{
			for(int L1=0; L1<ORB.Phi[it].getLmax()+1; L1++)
			{
				for(int L2=0; L2<ORB.Phi[it].getLmax()+1; L2++)
				{
					if( abs(L1-L2)<=L0 && L1+L2>=L0)
					{
						for(int N1=0; N1<ORB.Phi[it].getNchi(L1); N1++)
						{
							for(int N2=0; N2<ORB.Phi[it].getNchi(L2); N2++)
							{

							}
						}
					}
				}
			}
		}
	}
			
	delete[] this->aux_basis;
}


