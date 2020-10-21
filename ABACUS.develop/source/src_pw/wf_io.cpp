#include "wf_io.h"
#include "../src_pw/global.h"
//#include "../../src_develop/src_wannier/eximport.h"

void WF_io::write_wfc(const string &fn, const ComplexMatrix *psi)
{
    if (test_wf) TITLE("WF_io","write_wfc");

    ofstream ofs( fn.c_str() );

    //eximport exi;
    //exi.out_input(ofs);
    //exi.out_kpoints(ofs);
    //exi.out_igk(ofs);
    //exi.out_planewave(ofs);

    ofs << "\n<WAVEFUNC>";
    ofs << "\n" << NBANDS << " Number of bands." << endl;
    ofs << setprecision(6);
    for (int i=0; i<NBANDS; i++)
    {
        for (int ik=0; ik<kv.nks; ik++)
        {
            ofs << "\n" << ik;
            for (int ig=0; ig<kv.ngk[ik]; ig++)
            {
                if (ig%4==0) ofs << "\n";
                ofs << setw(15) << psi[ik](i, ig).real()
                << setw(15) << psi[ik](i, ig).imag();
            }
            ofs << "\n";
        }
    }
    ofs << "\n<WAVEFUNC>";

    ofs.close();
    return;
}

void WF_io::write_wfc2(const string &fn, const ComplexMatrix *psi)
{
    if (test_wf) TITLE("WF_io","write_wfc2");

	ofstream ofs( fn.c_str() );
	
	if(MY_RANK==0)
	{
		ofs.precision(10);
		ofs << ucell.lat0 << endl;

		ofs << ucell.latvec.e11 << " " << ucell.latvec.e12 << " " << ucell.latvec.e13 << endl;
		ofs << ucell.latvec.e21 << " " << ucell.latvec.e22 << " " << ucell.latvec.e23 << endl;
		ofs << ucell.latvec.e31 << " " << ucell.latvec.e32 << " " << ucell.latvec.e33 << endl;

		ofs << ucell.ntype << " ntype" << endl;
		for (int it=0; it<ucell.ntype; it++)
		{
			ofs << ucell.atoms[it].label << " label" << endl; // mohan add 2009-07-23
			ofs << ucell.atoms[it].na << " na" << endl;
			for (int ia=0; ia<ucell.atoms[it].na; ia++)
			{
				ofs << ucell.atoms[it].tau[ia].x
					<< " " << ucell.atoms[it].tau[ia].y
					<< " " << ucell.atoms[it].tau[ia].z << endl;
			}
		}
		// ecutwfc_jlq determine the jlq corresponding to plane wave calculation.
		ofs << pw.ecutwfc << " ecutwfc" << endl; // mohan add 2009-09-08
        ofs << kv.nkstot << " nks" << endl;
        ofs	<< NBANDS << " nbands" << endl;
        ofs << "\n<WEIGHT_OF_KPOINTS>";
	}

	// (2)
    for (int ik=0; ik<kv.nkstot; ik++)
    {
        double kx, ky, kz, wknow;
#ifdef __MPI
        const int pool = Pkpoints.whichpool[ik];
        const int iknow = ik - Pkpoints.startk_pool[MY_POOL];
        if (RANK_IN_POOL==0)
        {
            if (MY_POOL==0)
            {
                if (pool==0)
                {
                    kx = kv.kvec_c[ik].x;
                    ky = kv.kvec_c[ik].y;
                    kz = kv.kvec_c[ik].z;
                    wknow = kv.wk[ik];
                }
                else
                {
                    MPI_Status ierror;
                    MPI_Recv(&kx, 1, MPI_DOUBLE, Pkpoints.startpro_pool[pool], ik*4, MPI_COMM_WORLD,&ierror);
                    MPI_Recv(&ky, 1, MPI_DOUBLE, Pkpoints.startpro_pool[pool], ik*4+1, MPI_COMM_WORLD,&ierror);
                    MPI_Recv(&kz, 1, MPI_DOUBLE, Pkpoints.startpro_pool[pool], ik*4+2, MPI_COMM_WORLD,&ierror);
                    MPI_Recv(&wknow, 1, MPI_DOUBLE, Pkpoints.startpro_pool[pool], ik*4+3, MPI_COMM_WORLD,&ierror);
                }
            }
            else
            {
                if (MY_POOL == pool)
                {
                    MPI_Send(&kv.kvec_c[iknow].x, 1, MPI_DOUBLE, 0, ik*4, MPI_COMM_WORLD);
                    MPI_Send(&kv.kvec_c[iknow].y, 1, MPI_DOUBLE, 0, ik*4+1, MPI_COMM_WORLD);
                    MPI_Send(&kv.kvec_c[iknow].z, 1, MPI_DOUBLE, 0, ik*4+2, MPI_COMM_WORLD);
                    MPI_Send(&kv.wk[iknow], 1, MPI_DOUBLE, 0, ik*4+3, MPI_COMM_WORLD);
                }
            }
        }
        // this barrier is very important
        MPI_Barrier(MPI_COMM_WORLD);
#else
		if (MY_RANK==0)
        {
            kx = kv.kvec_c[ik].x;
            ky = kv.kvec_c[ik].y;
            kz = kv.kvec_c[ik].z;
            wknow = kv.wk[ik];
        }
#endif

        if (MY_RANK==0)
        {
            ofs << "\n" << kx << " " << ky << " " << kz;
            ofs << " " << wknow * 0.5;
        }
    }

    // (3)
    if (MY_RANK==0)
    {
        ofs << "\n</WEIGHT_OF_KPOINTS>" << endl;
		ofs.close();
    }

	// out put the wave functions in plane wave basis.
	for(int ip=0; ip<NPOOL; ip++)
	{
		if( MY_POOL == ip )
		{
			for(int ik=0; ik<kv.nks; ik++)
			{
				for(int ib=0; ib<NBANDS; ib++)
				{
					for( int id=0; id<NPROC_IN_POOL; id++)
					{
						if (RANK_IN_POOL == id)
						{
							ofstream ofs2( fn.c_str(),ios::app);
							ofs2 << scientific;
							ofs2 << setprecision(6);
							for (int ig=0; ig<kv.ngk[ik]; ig++)
							{
								if (ig%4==0) ofs2 << "\n";
								ofs2 << setw(15) << psi[ik](ib, ig).real()
									<< setw(15) << psi[ik](ib, ig).imag();
							} // end ig
							ofs2.close();
						}
					}// end id
				}// end ib
			}// end ik
		}
	}// end ip



	return;
}
