#include "hs_matrix.h"
#include "../src_pw/global.h"
#include "../src_lcao/local_orbital_elec.h"

void HS_Matrix::saving_HS(const double *Hloc, const double* Sloc, bool bit, const int &out_hs)
{	
	if(out_hs==1)
	{
		save_HS(Hloc, Sloc, bit);
	}
	else if(out_hs==2)
	{
		save_HS(Hloc, Sloc, bit);
	}
	else if(out_hs==3)
	{
		//please call individually
	}
	else if(out_hs==0)
	{
		// do nothing.
	}
	else
	{
		WARNING("Diago_LCAO_Matrix","unrecorganized out_hs value.");
	}
	return;
}


void HS_Matrix::save_HS_ccf(const int &iter, const int &Hnnz, const int *colptr_H, const int *rowind_H, 
	const double *nzval_H, const double *nzval_S, bool bit)
{
	TITLE("HS_Matrix","save_HS_ccf");

	if(DRANK!=0)return;
	
	stringstream ssh;
	stringstream sss;
	if(bit)
	{
		ssh << global_out_dir << "H_bit.ccf";
		sss << global_out_dir << "S_bit.ccf";
	}
	else
	{
		ssh << global_out_dir << "H" << Local_Orbital_Elec::iter << "_" << iter+1 << ".ccf";
		sss << global_out_dir << "S" << Local_Orbital_Elec::iter << "_" << iter+1 << ".ccf";
	}

	if(bit)
	{
	    FILE *g1 = fopen(ssh.str().c_str(),"wb");
        FILE *g2 = fopen(sss.str().c_str(),"wb");

        fwrite(&NLOCAL,sizeof(int),1,g1);
        fwrite(&Hnnz,sizeof(int),1,g1);
        fwrite(&NLOCAL,sizeof(int),1,g2);
        fwrite(&Hnnz,sizeof(int),1,g2);

        fclose(g1);
        fclose(g2);
	}

		
	if(!bit)
	{
		ofstream g1(ssh.str().c_str());
        ofstream g2(sss.str().c_str());

        g1 << NLOCAL << " " << Hnnz << endl;
        g2 << NLOCAL << " " << Hnnz << endl;

		for(int i=0; i<NLOCAL+1; ++i)
		{
			g1 << colptr_H[i] << " ";
			g2 << colptr_H[i] << " ";
		}
		g1 << endl;
		g2 << endl;

		for(int i=0; i<Hnnz; ++i)
		{
			g1 << rowind_H[i] << " ";
			g2 << rowind_H[i] << " ";
		}
		g1 << endl;
		g2 << endl;

		for(int i=0; i<Hnnz; ++i)
		{
			g1 << nzval_H[i] << " ";
			g2 << nzval_S[i] << " ";
		}
		g1 << endl;
		g2 << endl;

        g1.close();
        g2.close();
	}
	return;
}

// mohan add 2010/3/20, output H and S matrix, convinence for diagonalization
// test or save the middle information for next start.
void HS_Matrix::save_HS(const double *H, const double *S, bool bit)
{
    TITLE("HS_Matrix","save_HS_bit");
    timer::tick("HS_Matrix","save_HS_bit");
    OUT(ofs_running,"Dimension of H and S",NLOCAL);

    stringstream ssh;
    stringstream sss;

	if(bit)
	{
    	ssh << global_out_dir << "data-H-bit";
    	sss << global_out_dir << "data-S-bit";
	}
	else 
	{
    	ssh << global_out_dir << "data-H";
    	sss << global_out_dir << "data-S";
	}

    if (bit)
    {
#ifdef __MPI
        FILE *g1;
        FILE *g2;

        if (DRANK==0)
        {
            g1 = fopen(ssh.str().c_str(),"wb");
            g2 = fopen(sss.str().c_str(),"wb");
            fwrite(&NLOCAL,sizeof(int),1,g1);
            fwrite(&NLOCAL,sizeof(int),1,g2);
        }

        int ir,ic;
        for (int i=0; i<NLOCAL; i++)
        {
            double* lineH = new double[NLOCAL-i];
            double* lineS = new double[NLOCAL-i];
            ZEROS(lineH, NLOCAL-i);
            ZEROS(lineS, NLOCAL-i);

            ir = ParaO.trace_loc_row[i];
            if (ir>=0)
            {
                // data collection
                for (int j=i; j<NLOCAL; j++)
                {
                    ic = ParaO.trace_loc_col[j];
                    if (ic>=0)
                    {
                        lineH[j-i] = H[ir*ParaO.ncol+ic];
                        lineS[j-i] = S[ir*ParaO.ncol+ic];
                    }
                }
            }
            else
            {
                //do nothing
            }

            Parallel_Reduce::reduce_double_all(lineH,NLOCAL-i);
            Parallel_Reduce::reduce_double_all(lineS,NLOCAL-i);

            if (DRANK==0)
            {
                for (int j=i; j<NLOCAL; j++)
                {
                    fwrite(&lineH[j-i],sizeof(double),1,g1);
                    fwrite(&lineS[j-i],sizeof(double),1,g2);
                }
            }
            delete[] lineH;
            delete[] lineS;

            MPI_Barrier(DIAG_WORLD);
        }

        if (DRANK==0)
        {
            fclose(g1);
            fclose(g2);
        }
#else
        FILE *g1 = fopen(ssh.str().c_str(),"wb");
        FILE *g2 = fopen(sss.str().c_str(),"wb");

        fwrite(&NLOCAL,sizeof(int),1,g1);
        fwrite(&NLOCAL,sizeof(int),1,g2);

        for (int i=0; i<NLOCAL; i++)
        {
            for (int j=i; j<NLOCAL; j++)
            {
                fwrite(&H[i*NLOCAL+j],sizeof(double),1,g1);
                fwrite(&S[i*NLOCAL+j],sizeof(double),1,g2);
            }
        }
        fclose(g1);
        fclose(g2);
#endif
    } //end bit
    else
    {
#ifdef __MPI
        ofstream g1;
        ofstream g2;

        if (DRANK==0)
        {
            g1.open(ssh.str().c_str());
            g2.open(sss.str().c_str());
            g1 << NLOCAL;
            g2 << NLOCAL;
        }

        int ir,ic;
        for (int i=0; i<NLOCAL; i++)
        {
            double* lineH = new double[NLOCAL-i];
            double* lineS = new double[NLOCAL-i];
            ZEROS(lineH, NLOCAL-i);
            ZEROS(lineS, NLOCAL-i);

            ir = ParaO.trace_loc_row[i];
            if (ir>=0)
            {
                // data collection
                for (int j=i; j<NLOCAL; j++)
                {
                    ic = ParaO.trace_loc_col[j];
                    if (ic>=0)
                    {
                        lineH[j-i] = H[ir*ParaO.ncol+ic];
                        lineS[j-i] = S[ir*ParaO.ncol+ic];
                    }
                }
            }
            else
            {
                //do nothing
            }

            Parallel_Reduce::reduce_double_all(lineH,NLOCAL-i);
            Parallel_Reduce::reduce_double_all(lineS,NLOCAL-i);

            if (DRANK==0)
            {
                for (int j=i; j<NLOCAL; j++)
                {
                    g1 << " " << lineH[j-i];
                    g2 << " " << lineS[j-i];
                }
                g1 << endl;
                g2 << endl;
            }
            delete[] lineH;
            delete[] lineS;
        }

        if (DRANK==0);
        {
            g1.close();
            g2.close();
        }
#else
        ofstream g1(ssh.str().c_str());
        ofstream g2(sss.str().c_str());

        g1 << NLOCAL;
        g2 << NLOCAL;

        for (int i=0; i<NLOCAL; i++)
        {
            for (int j=i; j<NLOCAL; j++)
            {
                g1 << " " << H[i*NLOCAL+j];
                g2 << " " << S[i*NLOCAL+j];
            }
            g1 << endl;
            g2 << endl;
        }
        g1.close();
        g2.close();
#endif
    }

    timer::tick("HS_Matrix","save_HS_bit");
    return;
}


