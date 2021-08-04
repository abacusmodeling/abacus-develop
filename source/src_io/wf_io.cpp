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
    if (GlobalV::test_wf) TITLE("WF_io","write_wfc");

    ofstream ofs( fn.c_str() );

    //eximport exi;
    //exi.out_input(ofs);
    //exi.out_kpoints(ofs);
    //exi.out_igk(ofs);
    //exi.out_planewave(ofs);

    ofs << "\n<WAVEFUNC>";
    ofs << "\n" << GlobalV::NBANDS << " Number of bands." << std::endl;
    ofs << setprecision(6);
    for (int i=0; i<GlobalV::NBANDS; i++)
    {
        for (int ik=0; ik<GlobalC::kv.nks; ik++)
        {
            ofs << "\n" << ik;
            for (int ig=0; ig<GlobalC::kv.ngk[ik]; ig++)
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
    if (GlobalV::test_wf) TITLE("WF_io","write_wfc2"); 

    string * wfilename;
    wfilename=new string[GlobalC::kv.nkstot];
    for(int ik=0;ik<GlobalC::kv.nkstot;ik++)
    {
        stringstream wfss;
        if(GlobalC::wf.out_wf==1)
            wfss<<fn<<ik+1<<".txt";
        else
        {
            wfss<<fn<<ik+1<<".dat";
        }
        wfilename[ik]=wfss.str();
        if ( GlobalV::MY_RANK == 0 )
        {
            if(GlobalC::wf.out_wf==1)
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
    //if(GlobalV::MY_RANK!=0) std::cout.clear();
    //std::cout<<"Hello"<<std::endl;
    //if(GlobalV::MY_RANK!=0) std::cout.setstate(ios::failbit);
#ifdef __MPI
    MPI_Barrier(MPI_COMM_WORLD);

	
    // out put the wave functions in plane wave basis.
	for(int ip=0; ip<GlobalV::NPOOL; ip++)
	{
        if( GlobalV::MY_POOL == ip )
		{
#endif
			for(int ik=0; ik<GlobalC::kv.nks; ik++)
			{
                int ikngtot=0; //ikngtot: the total number of plane waves of ikpoint
                int ikstot=0;//ikstot : the index within all k-points
#ifdef __MPI
                MPI_Allreduce(&GlobalC::kv.ngk[ik],&ikngtot,1,MPI_INT,MPI_SUM,POOL_WORLD);
                ikstot=GlobalC::Pkpoints.startk_pool[ip]+ik;
#else
                ikngtot=GlobalC::kv.ngk[ik];
                ikstot=ik;
#endif
#ifdef __MPI
                //int ikstot=getink(ik,ip,GlobalC::kv.nkstot,GlobalV::NPOOL);
				for( int id=0; id<GlobalV::NPROC_IN_POOL; id++)
				{
					MPI_Barrier(POOL_WORLD);
                    if (GlobalV::RANK_IN_POOL == id)
					{
#else
                    int id=0;               
#endif
                    if(GlobalC::wf.out_wf==1)
                    {
                        ofstream ofs2( wfilename[ikstot].c_str(),ios::app);
                        if(id==0)
                        {
                            ofs2<<setprecision(6);
                            ofs2<<setw(10)<<"Kpoint"<<setw(10)<<"nKpoint"<<setw(10)<<"GlobalC::kv.x"<<setw(10)
                            <<"GlobalC::kv.y"<<setw(10)<<"GlobalC::kv.z"<<setw(10)<<"weight"<<setw(10)
                            <<"ngtot"<<setw(10)<<"nband"<<setw(10)<<"ecut"<<setw(10)<<"lat0"<<setw(10)<<"2pi/lat0"<<std::endl;
                            ofs2<<setw(10)<<ikstot+1<<setw(10)<<GlobalC::kv.nkstot<<setw(10)<<GlobalC::kv.kvec_c[ik].x<<setw(10)
                            <<GlobalC::kv.kvec_c[ik].y<<setw(10)<<GlobalC::kv.kvec_c[ik].z<<setw(10)<<GlobalC::kv.wk[ik]<<setw(10)
                            <<ikngtot<<setw(10)<<GlobalV::NBANDS<<setw(10)<<GlobalC::pw.ecutwfc<<setw(10)<<GlobalC::ucell.lat0<<setw(10)<<GlobalC::ucell.tpiba<<std::endl;
                            ofs2<<"\n<Reciprocal Lattice Vector>"<<std::endl;
                            ofs2<<setw(10)<<GlobalC::ucell.G.e11<<setw(10)<<GlobalC::ucell.G.e12<<setw(10)<<GlobalC::ucell.G.e13<<std::endl;
                            ofs2<<setw(10)<<GlobalC::ucell.G.e21<<setw(10)<<GlobalC::ucell.G.e22<<setw(10)<<GlobalC::ucell.G.e23<<std::endl;
                            ofs2<<setw(10)<<GlobalC::ucell.G.e31<<setw(10)<<GlobalC::ucell.G.e32<<setw(10)<<GlobalC::ucell.G.e33<<std::endl;
                            ofs2<<"<Reciprocal Lattice Vector>\n"<<std::endl;
                            ofs2<<"<G vectors>"<<std::endl;
                        }
                        for (int ig=0; ig<GlobalC::kv.ngk[ik]; ig++)
					    {
                            ofs2<<setw(10)<<gkk[ig].x<<setw(10)<<gkk[ig].y<<setw(10)<<gkk[ig].z<<std::endl;
						}
                        if(id==GlobalV::NPROC_IN_POOL-1)
                        {
                            ofs2<<"<G vectors>\n"<<std::endl;
                        }
						ofs2.close();
                    }
                    else
                    {
                        Rwstream wfs2( wfilename[ikstot],"a");
                        if(id==0)
                        {
                            wfs2<<int(72)<<ikstot+1<<GlobalC::kv.nkstot<<GlobalC::kv.kvec_c[ik].x
                                <<GlobalC::kv.kvec_c[ik].y<<GlobalC::kv.kvec_c[ik].z<<GlobalC::kv.wk[ik]
                                <<ikngtot<<GlobalV::NBANDS<<GlobalC::pw.ecutwfc<<GlobalC::ucell.lat0<<GlobalC::ucell.tpiba<<72; //4 int + 7 double is 72B
                            wfs2<<72<<GlobalC::ucell.G.e11<<GlobalC::ucell.G.e12<<GlobalC::ucell.G.e13
                                    <<GlobalC::ucell.G.e21<<GlobalC::ucell.G.e22<<GlobalC::ucell.G.e23
                                    <<GlobalC::ucell.G.e31<<GlobalC::ucell.G.e32<<GlobalC::ucell.G.e33<<72; //9 double is 72B
                        }
                        if(id==0)
                        {
                            wfs2<<ikngtot*8*3;
                        }
                         
                        for (int ig=0; ig<GlobalC::kv.ngk[ik]; ig++)
					    {
                            wfs2<<gkk[ig].x<<gkk[ig].y<<gkk[ig].z;
						}
                        if(id==GlobalV::NPROC_IN_POOL-1)
                        {
                            wfs2<<ikngtot*8*3;
                        } 
                        wfs2.close();
                    } 
#ifdef __MPI
					}
				}// end id
#endif
                for(int ib=0; ib<GlobalV::NBANDS; ib++)
				{
#ifdef __MPI
					for( int id=0; id<GlobalV::NPROC_IN_POOL; id++)
					{
						MPI_Barrier(POOL_WORLD); //qianrui add
                        if (GlobalV::RANK_IN_POOL == id)
						{
#else
                    int id=0;
#endif
                        if(GlobalC::wf.out_wf==1)
                        {
                            ofstream ofs2( wfilename[ikstot].c_str(),ios::app);
                            if(id==0)   ofs2 << "\n< Band "<<ib+1 <<" >" <<std::endl; 
							ofs2 << scientific;
							
                            for (int ig=0; ig<GlobalC::kv.ngk[ik]; ig++)
							{
								if (ig%4==0&&(ig!=0||id!=0)) ofs2 << "\n";
								ofs2 << setw(15) << psi[ik](ib, ig).real()
									<< setw(15) << psi[ik](ib, ig).imag();
							} // end ig
                            if(id==GlobalV::NPROC_IN_POOL-1)   ofs2 << "\n< Band "<<ib+1 <<" >" <<std::endl; 
							ofs2.close();
                        }
                        else
                        {
                            Rwstream wfs2(wfilename[ikstot],"a");
                            if(id==0) wfs2<<ikngtot*16;
                            for (int ig=0; ig<GlobalC::kv.ngk[ik]; ig++)
							{
								wfs2 << psi[ik](ib, ig).real() << psi[ik](ib, ig).imag();
							}
                            if(id==GlobalV::NPROC_IN_POOL-1) wfs2<<ikngtot*16;
                            wfs2.close();
                        }
#ifdef __MPI                      
						}
					}// end id
#endif
				}// end ib
			}// end ik
#ifdef __MPI
		}
	}// end ip
#endif


    delete [] wfilename;
	return;
}
