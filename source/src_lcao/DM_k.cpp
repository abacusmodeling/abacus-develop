#include "local_orbital_charge.h"
#include "../src_pw/global.h"
#include "../src_io/wf_local.h"
#include "../module_base/blas_connector.h"
#include "../src_parallel/parallel_common.h"
#include "../module_base/memory.h"
#include "../module_base/timer.h"

void Local_Orbital_Charge::allocate_DM_k(void)
{
    ModuleBase::TITLE("Local_Orbital_Charge","allocate_k");

    this->nnrg_now = GlobalC::GridT.nnrg;
    //xiaohui add 'GlobalV::OUT_LEVEL' line, 2015-09-16
    if(GlobalV::OUT_LEVEL != "m") ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"nnrg_last",nnrg_last);
    if(GlobalV::OUT_LEVEL != "m") ModuleBase::GlobalFunc::OUT(GlobalV::ofs_running,"nnrg_now",nnrg_now);

    if(this->init_DM_R)
    {
        assert(nnrg_last > 0);
        for(int is=0; is<GlobalV::NSPIN; is++)
        {
            delete[] DM_R[is];
        }
        delete[] DM_R;
        init_DM_R=false;
    }

    if(nnrg_now>0)
    {
        this->DM_R = new double*[GlobalV::NSPIN];
        for(int is=0; is<GlobalV::NSPIN; is++)
        {
            this->DM_R[is] = new double[nnrg_now];
            ModuleBase::GlobalFunc::ZEROS(DM_R[is], nnrg_now);
        }
        this->nnrg_last = nnrg_now;
        this->init_DM_R = true;
        ModuleBase::Memory::record("LocalOrbital_Charge","Density_Matrix",GlobalV::NSPIN*nnrg_now,"double");
    }
    else if(nnrg_now==0)
    {
        this->init_DM_R = false;
    }
    else
    {
        ModuleBase::WARNING_QUIT("Local_Orbital_Charge::allocate_k","check init_DM_R.");
    }

	// Peize Lin test 2019-01-16 
    this->init_dm_2d();

    return;
}

#include "record_adj.h"
inline void cal_DM_ATOM(
	const Grid_Technique &gt, 
	const std::complex<double> fac, 
	Record_adj RA,
 	const int ia1, 
	const int iw1_lo, 
	const int nw1, 
    const int gstart,
    std::complex<double>*** wfc_k_grid,
    std::complex<double>* WFC_PHASE,
	std::complex<double> **DM_ATOM)
{

    const char transa='N';
	const char transb='T';  
    const std::complex<double> alpha=1;
	const std::complex<double> beta=1;

    for(int ik=0; ik<GlobalC::kv.nks; ik++)
    {
        std::complex<double> **wfc = wfc_k_grid[ik];
        const int ispin = GlobalC::kv.isk[ik];
        int atom2start=0;

        for (int ia2 = 0; ia2 < RA.na_each[ia1]; ++ia2)
        {
            std::complex<double> *DM=&DM_ATOM[ispin][atom2start];
            const int T2 = RA.info[ia1][ia2][3];
            const int I2 = RA.info[ia1][ia2][4];
            Atom* atom2 = &GlobalC::ucell.atoms[T2];
            const int start2 = GlobalC::ucell.itiaiw2iwt(T2,I2,0);
            const int iw2_lo=gt.trace_lo[start2];
            const int nw2=atom2->nw;
            std::complex<double> exp_R= exp( fac * (
                        GlobalC::kv.kvec_d[ik].x * RA.info[ia1][ia2][0] + 
                        GlobalC::kv.kvec_d[ik].y * RA.info[ia1][ia2][1] + 
                        GlobalC::kv.kvec_d[ik].z * RA.info[ia1][ia2][2]  
                        ) );
            
            //ModuleBase::GlobalFunc::ZEROS(WFC_PHASE, GlobalV::NBANDS*nw1);
            int ibStart=0;
            int nRow=0;
            for(int ib=0; ib<GlobalV::NBANDS; ++ib)
            {
                const double wg_local=GlobalC::wf.wg(ik,ib);
                if(wg_local>0)
                {
                    if(nRow==0) ibStart=ib;
                    const int iline=nRow*nw1;
                    std::complex<double> phase=exp_R*wg_local;
                    for(int iw1=0; iw1<nw1; ++iw1)
					{
                        WFC_PHASE[iline+iw1]=phase*conj(wfc[ib][iw1_lo+iw1]);
					}
                    ++nRow;
                }
                else
				{
                    break;
				}
            } // ib
            zgemm_(&transa, &transb, &nw2, &nw1, &nRow, &alpha,
                &wfc[ibStart][iw2_lo], &gt.lgd, 
                WFC_PHASE, &nw1,
                &beta, DM, &nw2);           

            atom2start+=nw1*nw2;
        } // ia2
    } // ik
    return;
}

//added by zhengdy-soc, for non-collinear case
inline void cal_DM_ATOM_nc(
	const Grid_Technique &gt, 
	const std::complex<double> fac, 
	Record_adj RA,
	const int ia1, 
	const int iw1_lo, 
	const int nw1, 
    const int gstart,
    std::complex<double>*** wfc_k_grid,
    std::complex<double>* WFC_PHASE,
	std::complex<double> **DM_ATOM)
{

    if(GlobalV::NSPIN !=4 ) 
	{
		ModuleBase::WARNING_QUIT("Local_Orbital_Charge","NSPIN not match!");
	}

    const char transa='N';
	const char transb='T';  
    const std::complex<double> alpha=1;
	const std::complex<double> beta=1;
    int ispin=0;

    for(int is1=0;is1<2;is1++)
    {
        for(int is2=0;is2<2;is2++)
        {
            for(int ik=0; ik<GlobalC::kv.nks; ik++)
            {
                std::complex<double> **wfc = wfc_k_grid[ik];
                int atom2start=0;

                for (int ia2 = 0; ia2 < RA.na_each[ia1]; ++ia2)
                {
                    std::complex<double> *DM=&DM_ATOM[ispin][atom2start];
                    const int T2 = RA.info[ia1][ia2][3];
                    const int I2 = RA.info[ia1][ia2][4];
                    Atom* atom2 = &GlobalC::ucell.atoms[T2];
                    const int start2 = GlobalC::ucell.itiaiw2iwt(T2,I2,0);
                    const int iw2_lo=gt.trace_lo[start2]/GlobalV::NPOL + gt.lgd/GlobalV::NPOL*is2;
                    const int nw2=atom2->nw;
                    std::complex<double> exp_R= exp( fac * (
                                GlobalC::kv.kvec_d[ik].x * RA.info[ia1][ia2][0] + 
                                GlobalC::kv.kvec_d[ik].y * RA.info[ia1][ia2][1] + 
                                GlobalC::kv.kvec_d[ik].z * RA.info[ia1][ia2][2]  
                                ) );
            
            		//ModuleBase::GlobalFunc::ZEROS(WFC_PHASE, GlobalV::NBANDS*nw1);
                    int ibStart=0;
                    int nRow=0;
                    for(int ib=0; ib<GlobalV::NBANDS; ++ib)
                    {
                        const double w1=GlobalC::wf.wg(ik,ib);
                        if(w1>0)
                        {
                            if(nRow==0) 
							{
								ibStart=ib;
							}
                            const int iline=nRow*nw1;
                            std::complex<double> phase=exp_R*w1;
                            for(int iw1=0; iw1<nw1; ++iw1)
							{
                                WFC_PHASE[iline+iw1]=phase*conj(wfc[ib][iw1_lo+iw1 + gt.lgd/GlobalV::NPOL*is1]);
							}
                            ++nRow;
                        }
                        else
                            break;
                    } // ib
                    zgemm_(&transa, &transb, &nw2, &nw1, &nRow, &alpha,
                        &wfc[ibStart][iw2_lo], &gt.lgd, 
                        WFC_PHASE, &nw1,
                        &beta, DM, &nw2);           

                    atom2start+=nw1*nw2;
                } // ia2
            } // ik
            ispin++;
        }//is2
    }//is1
    return;
}


void Local_Orbital_Charge::cal_dk_k(const Grid_Technique &gt)
{
    ModuleBase::TITLE("Local_Orbital_Charge","cal_dk_k");
    ModuleBase::timer::tick("LCAO_Charge","cal_dk_k");  
    //int nnrg = 0;
    ModuleBase::Vector3<double> tau1;
	ModuleBase::Vector3<double> dtau;
        
    Record_adj RA;
    RA.for_grid(gt);

    int ca = 0;
    std::complex<double> fac = ModuleBase::TWO_PI * ModuleBase::IMAG_UNIT;

    std::complex<double> *WFC_PHASE=new std::complex<double>[GlobalV::NLOCAL*GlobalC::ucell.nwmax];
    
    int DM_ATOM_SIZE=1; 
    std::complex<double> **DM_ATOM=new std::complex<double> *[GlobalV::NSPIN];

    for(int is=0; is<GlobalV::NSPIN; ++is)
    {
        DM_ATOM[is]=new std::complex<double>[DM_ATOM_SIZE];
        ModuleBase::GlobalFunc::ZEROS(DM_ATOM[is], DM_ATOM_SIZE);
    }
    for(int T1=0; T1<GlobalC::ucell.ntype; T1++)
    {
        Atom* atom1 = &GlobalC::ucell.atoms[T1];
        for(int I1=0; I1<atom1->na; I1++)
        {
            const int iat = GlobalC::ucell.itia2iat(T1,I1);
            if(gt.in_this_processor[iat])
            {
                const int start1 = GlobalC::ucell.itiaiw2iwt(T1,I1,0);
                const int gstart = GlobalC::GridT.nlocstartg[iat];
                const int ng = GlobalC::GridT.nlocdimg[iat];
                const int iw1_lo=gt.trace_lo[start1]/GlobalV::NPOL;
                const int nw1=atom1->nw;

                if(DM_ATOM_SIZE<ng)
                {
                    DM_ATOM_SIZE=ng;
                    for(int is=0; is<GlobalV::NSPIN; ++is)
					{
                        delete[] DM_ATOM[is];
					}
                    for(int is=0; is<GlobalV::NSPIN; ++is)
					{
                        DM_ATOM[is]=new std::complex<double>[DM_ATOM_SIZE];
					}
                }
                for(int is=0; is<GlobalV::NSPIN; ++is)
				{
                    ModuleBase::GlobalFunc::ZEROS(DM_ATOM[is], ng);
				}
                ModuleBase::GlobalFunc::ZEROS(WFC_PHASE, GlobalV::NBANDS*nw1);
                if(GlobalV::NSPIN!=4)
				{
					cal_DM_ATOM(gt, fac, RA, ca, iw1_lo, nw1, gstart, this->LOWF->wfc_k_grid, WFC_PHASE, DM_ATOM);
				}
                else 
				{
					cal_DM_ATOM_nc(gt, fac, RA, ca, iw1_lo, nw1, gstart, this->LOWF->wfc_k_grid, WFC_PHASE, DM_ATOM);
				}
                ++ca;

                if(GlobalV::NSPIN!=4)
                {
                    for(int is=0; is<GlobalV::NSPIN; ++is)
                    {
                        for(int iv=0; iv<ng; ++iv)
                        {
                            this->DM_R[is][gstart+iv]=DM_ATOM[is][iv].real();
                        }
                    }
                }
                else
                {//zhengdy-soc
                    for(int iv=0; iv<ng; ++iv)
                    {
						//note: storage nondiagonal term as Re[] and Im[] respectly;
						this->DM_R[0][gstart+iv]=DM_ATOM[0][iv].real() + DM_ATOM[3][iv].real();
						if(GlobalV::NONCOLIN)
						{//GlobalV::DOMAG
							this->DM_R[1][gstart+iv]=DM_ATOM[1][iv].real() + DM_ATOM[2][iv].real();
							this->DM_R[2][gstart+iv]=DM_ATOM[1][iv].imag() - DM_ATOM[2][iv].imag();
							this->DM_R[3][gstart+iv]=DM_ATOM[0][iv].real() - DM_ATOM[3][iv].real();
						}
						else if(!GlobalV::NONCOLIN)//GlobalV::DOMAG_Z
						{
							this->DM_R[1][gstart+iv]= 0.0;
							this->DM_R[1][gstart+iv]= 0.0;
							this->DM_R[3][gstart+iv]=DM_ATOM[0][iv].real() - DM_ATOM[3][iv].real();
						}
						else//soc with no mag 
						{
							this->DM_R[1][gstart+iv]= 0.0;
							this->DM_R[2][gstart+iv]= 0.0;
							this->DM_R[3][gstart+iv]= 0.0;
						}
					}
				}
            } // if gt.in_this_processor
        }// I1
    }// T1


    //------------
    // for test
    //------------
/*  std::cout << std::setprecision(3);
    for(int i=0; i<nnrg_now; i++)

    for(int ik=0; ik<GlobalC::kv.nkstot; ++ik)
    {
        for(int ib=0; ib<GlobalV::NBANDS; ++ib)
        {
            std::cout << " ik=" << ik << " ib=" << ib << " occ=" << GlobalC::wf.wg(ik,ib) << " e=" << GlobalC::wf.ekb[ik][ib] << std::endl;
        }
    }

    for(int i=0; i<10; i++)
    {
        if(DM_R[0][i]>1.0e-8)
        {
            std::cout << " i=" << i << " DM_R=" << DM_R[0][i] << std::endl;
        }
    }
*/
    for(int i=0; i<GlobalV::NSPIN; ++i)
	{
        delete[] DM_ATOM[i];
	}
    delete[] DM_ATOM;
    delete[] WFC_PHASE;

    RA.delete_grid();//xiaohui add 2015-02-04

    ModuleBase::timer::tick("LCAO_Charge","cal_dk_k");  
    return;
}


