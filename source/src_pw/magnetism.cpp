#include "magnetism.h"
#ifndef __CELL
#include "global.h"
#endif

Magnetism::Magnetism()
{
    this->tot_magnetization = 0.0;
    this->abs_magnetization = 0.0;
    this->start_magnetization = new double[10];

	m_loc_ = new Vector3<double> [1];
	angle1_ = new double[1];
	angle2_ = new double[1];
}

Magnetism::~Magnetism()
{
    delete[] this->start_magnetization;
	delete[] m_loc_;
	delete[] angle1_;
	delete[] angle2_;
}

#ifndef __CELL
void Magnetism::compute_magnetization()
{
    if (GlobalV::NSPIN==2)
    {
        this->tot_magnetization = 0.00;
        this->abs_magnetization = 0.00;

		//CHR.check_ne(CHR.rho[0]);
		//CHR.check_ne(CHR.rho[1]);
        for (int ir=0; ir<GlobalC::pw.nrxx; ir++)
        {
            double diff = CHR.rho[0][ir] - CHR.rho[1][ir];
            this->tot_magnetization += diff;
            this->abs_magnetization += abs(diff);
        }

        Parallel_Reduce::reduce_double_pool( this->tot_magnetization );
        Parallel_Reduce::reduce_double_pool( this->abs_magnetization );
        this->tot_magnetization *= GlobalC::ucell.omega / GlobalC::pw.ncxyz;
        this->abs_magnetization *= GlobalC::ucell.omega / GlobalC::pw.ncxyz;

		OUT(GlobalV::ofs_running,"total magnetism (Bohr mag/cell)",this->tot_magnetization);
		OUT(GlobalV::ofs_running,"absolute magnetism (Bohr mag/cell)",this->abs_magnetization);
		
		if(GlobalV::TWO_EFERMI)
		{
			OUT(GlobalV::ofs_running,"nelup",get_nelup());
			OUT(GlobalV::ofs_running,"neldw",get_neldw());
		}
		else
		{
			OUT(GlobalV::ofs_running,"nelec",CHR.nelec);
		}

//        cout << "\n tot_mag = " << setprecision(6) << this->tot_magnetization << " Bohr mag/cell" << endl;
  //      cout << " abs_mag = " << setprecision(6) << this->abs_magnetization << " Bohr mag/cell" << endl;
    }
	// noncolliear :
	else if(GlobalV::NSPIN==4)
	{
		for(int i=0;i<3;i++)this->tot_magnetization_nc[i] = 0.00;
		this->abs_magnetization = 0.00;
		for (int ir=0; ir<GlobalC::pw.nrxx; ir++)
		{
			double diff = sqrt(pow(CHR.rho[1][ir], 2) + pow(CHR.rho[2][ir], 2) +pow(CHR.rho[3][ir], 2));
 
			for(int i=0;i<3;i++)this->tot_magnetization_nc[i] += CHR.rho[i+1][ir];
			this->abs_magnetization += abs(diff);
		}
		Parallel_Reduce::reduce_double_pool( this->tot_magnetization_nc, 3 );
		Parallel_Reduce::reduce_double_pool( this->abs_magnetization );

		for(int i=0;i<3;i++)this->tot_magnetization_nc[i] *= GlobalC::ucell.omega / GlobalC::pw.ncxyz;
		this->abs_magnetization *= GlobalC::ucell.omega / GlobalC::pw.ncxyz;
		GlobalV::ofs_running<<"total magnetism (Bohr mag/cell)"<<'\t'<<this->tot_magnetization_nc[0]<<'\t'<<this->tot_magnetization_nc[1]<<'\t'<<this->tot_magnetization_nc[2]<<'\n';
		OUT(GlobalV::ofs_running,"absolute magnetism (Bohr mag/cell)",this->abs_magnetization);
	}

    return;
}


double Magnetism::get_nelup(void)
{
	double nelup = 0.0;
	if(GlobalV::TWO_EFERMI)
	{
//===============================================================
//  this type of electrons are used as "fixed" magnetization.
//===============================================================
		nelup = 0.5 * CHR.nelec + 0.5 * tot_magnetization;
	}
	else
	{
		nelup = 0.5 * CHR.nelec;
	}
    return nelup;

    // for constrained magnetism calculation : not used now
    // nelup = ( nelec + mcons(3,1) ) * 0.5D0

//	double nelup = 0.0;
//	for(int i=0; i<GlobalC::pw.ntype; i++)
//	{
//		nelup += CHR.nelec * (1.0+start_magnetization[i])/2.0/GlobalC::pw.ntype;
//	}
//	return nelup;
}


double Magnetism::get_neldw(void)
{
	double neldw = 0.0;
	if(GlobalV::TWO_EFERMI)
	{
//===============================================================
//  this type of electrons are used as "fixed" magnetization.
//===============================================================
		neldw = 0.5 * CHR.nelec - 0.5 * tot_magnetization;
	}
	else
	{
		neldw = 0.5 * CHR.nelec;
	}
    return neldw ;

//	double neldw = 0.0;
//	for(int i=0; i<GlobalC::pw.ntype; i++)
//	{
//		neldw += GlobalC::pw.nelec * (1.0-start_magnetization[i])/2.0/GlobalC::pw.ntype;
//	}
//	return neldw;
}

void Magnetism::cal_ux(const int ntype)
{
	double amag, uxmod;
	int starting_it;
	bool is_paraller;
	//do not sign feature in teh general case
	lsign_ = false;
	ZEROS(ux_, 3);

	starting_it = 0;
	for(int it = 0;it<ntype;it++)
	{
		amag = pow(m_loc_[it].x,2) + pow(m_loc_[it].y,2) + pow(m_loc_[it].z,2);
		if(amag > 1e-6)
		{
			ux_[0] = m_loc_[it].x;
			ux_[1] = m_loc_[it].y;
			ux_[2] = m_loc_[it].z;
			starting_it = it;
			lsign_ = true;
			break;
		}
	}
	//initial magnetizations should be parallel
	for(int it = starting_it+1; it<ntype;it++)
	{
		lsign_ = lsign_ && judge_parallel(ux_, m_loc_[it]);
	}
	if(lsign_)
	{
		uxmod =  pow(ux_[0],2) + pow(ux_[1],2) +pow(ux_[2],2);
		if(uxmod<1e-6) 
		{
			WARNING_QUIT("cal_ux","wrong uxmod");
		}
		for(int i = 0;i<3;i++)
		{
			ux_[i] *= 1/sqrt(uxmod);
		}
		//       cout<<"    Fixed quantization axis for GGA: "
		//<<setw(10)<<ux[0]<<"  "<<setw(10)<<ux[1]<<"  "<<setw(10)<<ux[2]<<endl;
	}
	return;
}

bool Magnetism::judge_parallel(double a[3], Vector3<double> b)
{
   bool jp=false;
   double cross;
   cross = pow((a[1]*b.z-a[2]*b.y),2) +  pow((a[2]*b.x-a[0]*b.z),2) + pow((a[0]*b.y-a[1]*b.x),2);
   jp = (fabs(cross)<1e-6);
   return jp;
}
#endif