#include "wavefunc_in_pw.h"
#include <cstring>		// Peize Lin fix bug about strcmp 2016-08-02

void Wavefunc_in_pw::make_table_q(std::vector<string> &fn, realArray &table_local)
{
	TITLE("Wavefunc_in_pw","make_table_q");

	if( fn.size() != ucell.ntype )
	{
		WARNING_QUIT("Wavefunc_in_pw::make_table_q","maybe NUMERICAL_ORBITAL is not read in, please check.");
	}

	for(int it=0; it<ucell.ntype; it++)
	{
		ifstream in(fn[it].c_str());
		if(!in) 
		{
			ofs_warning << " File name : " << fn[it] << endl;
			WARNING_QUIT("Wavefunc_in_pw::make_table_q","Can not find file.");
		}
		else 
		{
			stringstream ss;
			ss << "Orbital of species " << ucell.atoms[it].label;
			OUT(ofs_running,ss.str(),fn[it]);
		}
		in.close();	
	}

	table_local.zero_out();
	for(int it=0; it<ucell.ntype; it++)
	{	
		int ic=0;
		for(int L=0; L<ucell.atoms[it].nwl+1; L++)
		{
			for(int N=0; N<ucell.atoms[it].l_nchi[L]; N++)
			{
				ofs_running << " L=" << L << " N=" << N;
				ifstream in(fn[it].c_str());
				if (!in) 
				{
					ofs_warning << " File name : " << fn[it] << endl;
					WARNING_QUIT("Wavefunc_in_pw::make_table_q","Can not find file.");
				}
				int meshr;
				double dr; // only used in uniform grid
                                char word[80];     // pengfei Li add 15-1-31
                                while (in.good())
                                {
                                    in >> word;
                                    if (std::strcmp(word , "END") == 0)		// Peize Lin fix bug about strcmp 2016-08-02
                                    {
                                        break;
                                    }
                                }

				CHECK_NAME(in, "Mesh");
				in >> meshr;
				int meshr_read = meshr;
				if(meshr%2==0)
				{
					++meshr;
				}
				ofs_running << " meshr=" << meshr; 

				CHECK_NAME(in, "dr");
				in >> dr;
				ofs_running << " dr=" << dr;
				
				double* radial = new double[meshr];
				double *psi = new double[meshr];
				double* psir = new double[meshr];
				double* rab = new double[meshr];

				ZEROS(radial, meshr);
				ZEROS(psi, meshr);
				ZEROS(psir, meshr);
				ZEROS(rab, meshr);
				for(int ir=0; ir<meshr; ir++)
				{
					rab[ir] = dr;
					// plus one because we can't read in r = 0 term now.
					radial[ir] = ir*dr;  //mohan modify 2010-04-19
				}
				ofs_running << " Rmax(Angstrom)=" << radial[meshr-1] << endl;

				string name1, name2, name3;
				int tmp_it, tmp_l ,tmp_n;
				bool find = false;

				while( !find )
				{
					if(in.eof())
					{
						ofs_warning << "\n Can't find l="
						<< L << " n=" << N << " orbital." << endl;
						WARNING_QUIT("Control_Overlap","Read_PAO");
					}
					in >> name1 >> name2 >> name3;
					assert( name1 == "Type" );
					in >> tmp_it >> tmp_l >> tmp_n;
					if( L == tmp_l && N == tmp_n )
					{
						// meshr_read is different from meshr if meshr is even number.
						for(int ir=0; ir<meshr_read; ir++)
						{
							in >> psi[ir];
							//psi[ir] = 1.0; //hahaha
							psir[ir] = psi[ir] * radial[ir];
						}
						find = true;
					}
					else
					{
						 double no_use;
						 for(int ir=0; ir<meshr_read; ir++)
						 {
							 in >> no_use;
						 }
					}		
				}
				double* table = new double[NQX];
				Wavefunc_in_pw::integral(meshr, psir, radial, rab, L, table);
				for(int iq=0; iq<NQX; iq++)
				{
					double energy_q = pow(iq * DQ,2);
					table_local(it,ic,iq) = table[iq];//* Wavefunc_in_pw::smearing(energy_q,150,0.666666);
				}
				delete[] table;
				delete[] radial;
				delete[] rab;
				delete[] psi;
				delete[] psir;
				++ic;
			}// N 
		}// L
	}// T


	if(MY_RANK==0)
	{
		for(int it=0; it<ucell.ntype; it++)
		{
			stringstream ss;
			ss << global_out_dir << ucell.atoms[it].label << "/LOCAL_G.dat";
			ofstream ofs(ss.str().c_str());
			for(int iq=0; iq<NQX; iq++)
			{
				int ic=0;
				double energy_q = pow((double)iq*DQ,2);
				ofs << energy_q; // unit (Ry)
				for(int L=0; L<ucell.atoms[it].nwl+1; L++)
				{
					for(int N=0; N<ucell.atoms[it].l_nchi[L]; N++)
					{
						ofs << " " << table_local(it,ic,iq);
						++ic;
					}
				}
				ofs << endl;
			}
			ofs.close();
		}
	}
		
	return;
}

double Wavefunc_in_pw::smearing(const double &energy_x,
                               const double &ecut,
                               const double &beta)
{
    double w = 0.0;
    const double beta_e = beta * ecut ;

    if (beta >= 1.0 || beta<0 )
    {
        WARNING_QUIT("wavefunc_in_pw::smearing", "beta must between 0 ~ 1 ");
    }

    if (energy_x < beta_e)
    {
        w = 1.0;
    }
    else if (energy_x >= beta_e && energy_x <= ecut)
    {
        const double arg = PI * (ecut - energy_x) * 0.5 / (1-beta) / ecut ;
        // const double sin_arg = sin(arg);  // gong 2009. 7. 12 , correct
        // w = sin_arg*sin_argi ;
        w = 0.5 * (1 - cos(2.0 * arg));
    }
    else if (energy_x > ecut)
    {
        w = 0.0 ;
    }
    return w ;

}




void Wavefunc_in_pw::integral(const int meshr, const double *psir, const double *r,
const double *rab, const int &l, double* table)
{
	const double pref = FOUR_PI / sqrt(ucell.omega);
	
	double *inner_part = new double[meshr];
	for(int ir=0; ir<meshr; ir++)
	{
		inner_part[ir] = psir[ir] * psir[ir];
	}
	
	double unit = 0.0;
	Mathzone::Simpson_Integral(meshr, inner_part, rab, unit);
	delete[] inner_part;
	OUT(ofs_running,"normalize unit",unit);

	double *aux = new double[meshr];
	double *vchi = new double[meshr];
	for (int iq=0; iq<NQX; iq++)
	{
		const double q = DQ * iq;
		Mathzone::Spherical_Bessel(meshr, r, q, l, aux);
		for (int ir = 0;ir < meshr;ir++)
		{
			vchi[ir] = psir[ir] * aux[ir] * r[ir];
		}
		
		double vqint = 0.0;
		Mathzone::Simpson_Integral(meshr, vchi, rab, vqint);

		table[iq] =  vqint * pref;
	}
	delete[] aux;
	delete[] vchi;
	return;
}


void Wavefunc_in_pw::produce_local_basis_in_pw(const int &ik,ComplexMatrix &psi, const realArray &table_local)
{
	TITLE("Wavefunc_in_pw","produce_local_basis_in_pw");
	assert(ik>=0);
	const int npw = kv.ngk[ik];
	const int total_lm = ( ucell.lmax + 1) * ( ucell.lmax + 1);
	matrix ylm(total_lm, npw);

	Vector3<double> *gk = new Vector3<double>[npw];
	for(int ig=0;ig<npw;ig++)
	{
		gk[ig] = wf.get_1qvec_cartesian(ik, ig);
	}
	
	Mathzone::Ylm_Real(total_lm, npw, gk, ylm);

	int index = 0;
	double *flq = new double[npw];
	int iwall=0;
	for (int it = 0;it < ucell.ntype;it++)
	{
		for (int ia = 0;ia < ucell.atoms[it].na;ia++)
		{
			complex<double> *sk = wf.get_sk(ik, it, ia);
			int ic=0;
			for(int L = 0; L < ucell.atoms[it].nwl+1; L++)
			{
				complex<double> lphase = pow(NEG_IMAG_UNIT, L); //mohan 2010-04-19
				for(int N=0; N < ucell.atoms[it].l_nchi[L]; N++)
				{
//					ofs_running << " it=" << it << " ia=" << ia << " L=" << L << " N=" << N << endl;
					
					for(int ig=0; ig<npw; ig++)
					{
						flq[ig] = Mathzone::Polynomial_Interpolation(table_local,
						it, ic, NQX, DQ, gk[ig].norm() * ucell.tpiba );
					}

					for(int m=0; m<2*L+1; m++)
					{
						const int lm = L*L+m;
						for(int ig=0; ig<npw; ig++)
						{
							psi(iwall, ig) =
							lphase * sk[ig] * ylm(lm, ig) * flq[ig];
						}	
						++iwall;
					}
					++ic;
				}
			}
			delete[] sk;
		}
	}
	assert(iwall == NLOCAL);
	delete[] flq;
	delete[] gk;
}


void Wavefunc_in_pw::produce_local_basis_q_in_pw(const int &ik, ComplexMatrix &psi, const realArray &table_local, Vector3<double> q)   // pengfei 2016-11-23
{
	TITLE("Wavefunc_in_pw","produce_local_basis_in_pw");
	assert(ik>=0);
	const int npw = kv.ngk[ik];
	const int total_lm = ( ucell.lmax + 1) * ( ucell.lmax + 1);
	matrix ylm(total_lm, npw);

	Vector3<double> *gkq = new Vector3<double>[npw];
        
	for(int ig=0;ig<npw;ig++)
	{
		gkq[ig] = wf.get_1qvec_cartesian(ik, ig) + q;
	}

	Mathzone::Ylm_Real(total_lm, npw, gkq, ylm);

	int index = 0;
	double *flq = new double[npw];
	int iwall=0;
	for (int it = 0;it < ucell.ntype;it++)
	{
		for (int ia = 0;ia < ucell.atoms[it].na;ia++)
		{
			complex<double> *skq = wf.get_skq(ik, it, ia, q);
			int ic=0;
			for(int L = 0; L < ucell.atoms[it].nwl+1; L++)
			{
				complex<double> lphase = pow(NEG_IMAG_UNIT, L); //mohan 2010-04-19
				for(int N=0; N < ucell.atoms[it].l_nchi[L]; N++)
				{
					for(int ig=0; ig<npw; ig++)
					{
                                                if(gkq[ig].norm() * ucell.tpiba > ((NQX-4) * DQ) )
                                                {
                                                   flq[ig] = 0.0;
                                                }
                                                else
                                                {
						   flq[ig] = Mathzone::Polynomial_Interpolation(table_local, it, ic, NQX, DQ, gkq[ig].norm() * ucell.tpiba );
                                                }
					}

					for(int m=0; m<2*L+1; m++)
					{
						const int lm = L*L+m;
						for(int ig=0; ig<npw; ig++)
						{
							psi(iwall, ig) =
							lphase * skq[ig] * ylm(lm, ig) * flq[ig];
						}	
						++iwall;
					}
					++ic;
				}
			}
			delete[] skq;
		}
	}

	assert(iwall == NLOCAL);
	delete[] flq;
	delete[] gkq;
}
