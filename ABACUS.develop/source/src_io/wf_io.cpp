#include "wf_io.h"
#include "../src_pw/global.h"
//#include "../../src_develop/src_wannier/eximport.h"

//qianrui add 2020-10-15
/*inline int getink(const int &ik,const int &ipool,const int &nktot,const int &npool)
{
    int nkp=nktot/npool;
    int rem=nktot%npool;
    if(ipool<rem)
    {
        return ipool*nkp+ipool+ik;
    }
    else
    {
        return ipool*nkp+rem+ik;       
    }
    
}*/

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

void WF_io::write_wfc2(const string &fn, const ComplexMatrix *psi, const Vector3<double> *gkk)
{
    if (test_wf) TITLE("WF_io","write_wfc2"); 

    string * wfilename;
    wfilename=new string[kv.nkstot];
    for(int ik=0;ik<kv.nkstot;ik++)
    {
        stringstream wfss;
        if(wf.out_wf==1)
            wfss<<fn<<ik+1<<".txt";
        else
        {
            wfss<<fn<<ik+1<<".dat";
        }
        wfilename[ik]=wfss.str();
        if ( MY_RANK == 0 )
        {
            if(wf.out_wf==1)
            {
                ofstream ofs(wfss.str().c_str()); //clear all wavefunc files.
                ofs.close();
            }
            else
            {
                Rwstream wfs(wfss.str(),"w");
                wfs.close();
            }
        }
    }
    //if(MY_RANK!=0) cout.clear();
    //cout<<"Hello"<<endl;
    //if(MY_RANK!=0) cout.setstate(ios::failbit);
    MPI_Barrier(MPI_COMM_WORLD);
	
    // out put the wave functions in plane wave basis.
	for(int ip=0; ip<NPOOL; ip++)
	{
        if( MY_POOL == ip )
		{
			for(int ik=0; ik<kv.nks; ik++)
			{
                int ikngtot; //ikngtot: the total number of plane waves of ikpoint
                MPI_Allreduce(&kv.ngk[ik],&ikngtot,1,MPI_INT,MPI_SUM,POOL_WORLD);
                int ikstot=Pkpoints.startk_pool[ip]+ik;//ikstot : the index within all k-points
                //int ikstot=getink(ik,ip,kv.nkstot,NPOOL);
					for( int id=0; id<NPROC_IN_POOL; id++)
					{
						MPI_Barrier(POOL_WORLD);
                        if (RANK_IN_POOL == id)
						{
                        if(wf.out_wf==1)
                        {
                            ofstream ofs2( wfilename[ikstot].c_str(),ios::app);
                            if(id==0)
                            {
                                ofs2<<setprecision(6);
                                ofs2<<setw(10)<<"Kpoint"<<setw(10)<<"nKpoint"<<setw(10)<<"kv.x"<<setw(10)
                                <<"kv.y"<<setw(10)<<"kv.z"<<setw(10)<<"weight"<<setw(10)
                                <<"ngtot"<<setw(10)<<"nband"<<setw(10)<<"ecut"<<setw(10)<<"lat0"<<setw(10)<<"2pi/lat0"<<endl;
                                ofs2<<setw(10)<<ikstot+1<<setw(10)<<kv.nkstot<<setw(10)<<kv.kvec_c[ik].x<<setw(10)
                                <<kv.kvec_c[ik].y<<setw(10)<<kv.kvec_c[ik].z<<setw(10)<<kv.wk[ik]<<setw(10)
                                <<ikngtot<<setw(10)<<NBANDS<<setw(10)<<pw.ecutwfc<<setw(10)<<ucell.lat0<<setw(10)<<ucell.tpiba<<endl;
                                ofs2<<"\n<Reciprocal Lattice Vector>"<<endl;
                                ofs2<<setw(10)<<ucell.G.e11<<setw(10)<<ucell.G.e12<<setw(10)<<ucell.G.e13<<endl;
                                ofs2<<setw(10)<<ucell.G.e21<<setw(10)<<ucell.G.e22<<setw(10)<<ucell.G.e23<<endl;
                                ofs2<<setw(10)<<ucell.G.e31<<setw(10)<<ucell.G.e32<<setw(10)<<ucell.G.e33<<endl;
                                ofs2<<"<Reciprocal Lattice Vector>\n"<<endl;
                                ofs2<<"<G vectors>"<<endl;
                            }
                            for (int ig=0; ig<kv.ngk[ik]; ig++)
						    {
                                ofs2<<setw(10)<<gkk[ig].x<<setw(10)<<gkk[ig].y<<setw(10)<<gkk[ig].z<<endl;
							}
                            if(id==NPROC_IN_POOL-1)
                            {
                                ofs2<<"<G vectors>\n"<<endl;
                            }
							ofs2.close();
                        }
                        else
                        {
                            Rwstream wfs2( wfilename[ikstot],"a");
                            if(id==0)
                            {
                                wfs2<<int(72)<<ikstot+1<<kv.nkstot<<kv.kvec_c[ik].x
                                    <<kv.kvec_c[ik].y<<kv.kvec_c[ik].z<<kv.wk[ik]
                                    <<ikngtot<<NBANDS<<pw.ecutwfc<<ucell.lat0<<ucell.tpiba<<72; //4 int + 7 double is 72B
                                wfs2<<72<<ucell.G.e11<<ucell.G.e12<<ucell.G.e13
                                        <<ucell.G.e21<<ucell.G.e22<<ucell.G.e23
                                        <<ucell.G.e31<<ucell.G.e32<<ucell.G.e33<<72; //9 double is 72B
                            }
                            if(id==0)
                            {
                                wfs2<<ikngtot*8*3;
                            }
                             
                            for (int ig=0; ig<kv.ngk[ik]; ig++)
						    {
                                wfs2<<gkk[ig].x<<gkk[ig].y<<gkk[ig].z;
							}
                            if(id==NPROC_IN_POOL-1)
                            {
                                wfs2<<ikngtot*8*3;
                            } 
                            wfs2.close();
                        } 
						}
					}// end id
                for(int ib=0; ib<NBANDS; ib++)
				{
					for( int id=0; id<NPROC_IN_POOL; id++)
					{
						MPI_Barrier(POOL_WORLD); //qianrui add
                        if (RANK_IN_POOL == id)
						{
                        if(wf.out_wf==1)
                        {
                            ofstream ofs2( wfilename[ikstot].c_str(),ios::app);
                            if(id==0)   ofs2 << "\n< Band "<<ib+1 <<" >" <<endl; 
							ofs2 << scientific;
							
                            for (int ig=0; ig<kv.ngk[ik]; ig++)
							{
								if (ig%4==0&&(ig!=0||id!=0)) ofs2 << "\n";
								ofs2 << setw(15) << psi[ik](ib, ig).real()
									<< setw(15) << psi[ik](ib, ig).imag();
							} // end ig
                            if(id==NPROC_IN_POOL-1)   ofs2 << "\n< Band "<<ib+1 <<" >" <<endl; 
							ofs2.close();
                        }
                        else
                        {
                            Rwstream wfs2(wfilename[ikstot],"a");
                            if(id==0) wfs2<<ikngtot*16;
                            for (int ig=0; ig<kv.ngk[ik]; ig++)
							{
								wfs2 << psi[ik](ib, ig).real() << psi[ik](ib, ig).imag();
							}
                            if(id==NPROC_IN_POOL-1) wfs2<<ikngtot*16;
                            wfs2.close();
                        }
                        
						}
					}// end id
				}// end ib
			}// end ik
		}
	}// end ip


    delete [] wfilename;
	return;
}
