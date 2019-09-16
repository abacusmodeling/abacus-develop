#include "magnetism.h"
#include "global.h"

Magnetism::Magnetism()
{
    this->tot_magnetization = 0.0;
    this->abs_magnetization = 0.0;
    this->start_magnetization = new double[10];
}

Magnetism::~Magnetism()
{
    delete[] this->start_magnetization;
}

void Magnetism::compute_magnetization()
{
    if (NSPIN==2)
    {
        this->tot_magnetization = 0.00;
        this->abs_magnetization = 0.00;

		//chr.check_ne(chr.rho[0]);
		//chr.check_ne(chr.rho[1]);
        for (int ir=0; ir<pw.nrxx; ir++)
        {
            double diff = chr.rho[0][ir] - chr.rho[1][ir];
            this->tot_magnetization += diff;
            this->abs_magnetization += abs(diff);
        }

        Parallel_Reduce::reduce_double_pool( this->tot_magnetization );
        Parallel_Reduce::reduce_double_pool( this->abs_magnetization );
        this->tot_magnetization *= ucell.omega / pw.ncxyz;
        this->abs_magnetization *= ucell.omega / pw.ncxyz;

		OUT(ofs_running,"total magnetism (Bohr mag/cell)",this->tot_magnetization);
		OUT(ofs_running,"absolute magnetism (Bohr mag/cell)",this->abs_magnetization);
		
		if(TWO_EFERMI)
		{
			OUT(ofs_running,"nelup",get_nelup());
			OUT(ofs_running,"neldw",get_neldw());
		}
		else
		{
			OUT(ofs_running,"nelec",ucell.nelec);
		}

//        cout << "\n tot_mag = " << setprecision(6) << this->tot_magnetization << " Bohr mag/cell" << endl;
  //      cout << " abs_mag = " << setprecision(6) << this->abs_magnetization << " Bohr mag/cell" << endl;
    }
	// noncolliear :
	else if(NSPIN==4)
	{
		for(int i=0;i<3;i++)this->tot_magnetization_nc[i] = 0.00;
		this->abs_magnetization = 0.00;
		for (int ir=0; ir<pw.nrxx; ir++)
		{
			double diff = sqrt(pow(chr.rho[1][ir], 2) + pow(chr.rho[2][ir], 2) +pow(chr.rho[3][ir], 2));
 
			for(int i=0;i<3;i++)this->tot_magnetization_nc[i] += chr.rho[i+1][ir];
			this->abs_magnetization += abs(diff);
		}
		Parallel_Reduce::reduce_double_pool( this->tot_magnetization_nc, 3 );
		Parallel_Reduce::reduce_double_pool( this->abs_magnetization );

		for(int i=0;i<3;i++)this->tot_magnetization_nc[i] *= ucell.omega / pw.ncxyz;
		this->abs_magnetization *= ucell.omega / pw.ncxyz;
		ofs_running<<"total magnetism (Bohr mag/cell)"<<'\t'<<this->tot_magnetization_nc[0]<<'\t'<<this->tot_magnetization_nc[1]<<'\t'<<this->tot_magnetization_nc[2]<<'\n';
		OUT(ofs_running,"absolute magnetism (Bohr mag/cell)",this->abs_magnetization);
	}

    return;
}


double Magnetism::get_nelup(void)
{
	double nelup = 0.0;
	if(TWO_EFERMI)
	{
//===============================================================
//  this type of electrons are used as "fixed" magnetization.
//===============================================================
		nelup = 0.5 * ucell.nelec + 0.5 * tot_magnetization;
	}
	else
	{
		nelup = 0.5 * ucell.nelec;
	}
    return nelup;

    // for constrained magnetism calculation : not used now
    // nelup = ( nelec + mcons(3,1) ) * 0.5D0

//	double nelup = 0.0;
//	for(int i=0; i<pw.ntype; i++)
//	{
//		nelup += ucell.nelec * (1.0+start_magnetization[i])/2.0/pw.ntype;
//	}
//	return nelup;
}


double Magnetism::get_neldw(void)
{
	double neldw = 0.0;
	if(TWO_EFERMI)
	{
//===============================================================
//  this type of electrons are used as "fixed" magnetization.
//===============================================================
		neldw = 0.5 * ucell.nelec - 0.5 * tot_magnetization;
	}
	else
	{
		neldw = 0.5 * ucell.nelec;
	}
    return neldw ;

//	double neldw = 0.0;
//	for(int i=0; i<pw.ntype; i++)
//	{
//		neldw += pw.nelec * (1.0-start_magnetization[i])/2.0/pw.ntype;
//	}
//	return neldw;
}
