#include "pw_basis_k.h"
PW_Basis_K::PW_Basis_K()
{
    nks = 1;
    kvec_d = NULL;
    ngk = NULL;
    GR_index = NULL;
}
PW_Basis_K::~PW_Basis_K()
{
    if(kvec_d != NULL) delete[] kvec_d;
    if(ngk != NULL) delete[] ngk;
    if(GR_index != NULL) delete[] GR_index;
}

void PW_Basis_K:: initparameters(
    bool gamma_only_in,
    double ecut_in,
    double gk_ecut_in,
    int nks_in, //number of k points in this pool
    ModuleBase::Vector3<double> *kvec_d_in, // Direct coordinates of k points
    int poolnproc_in, // Number of processors in this pool
    int poolrank_in, // Rank in this pool
    int distribution_type_in,
)
{
    initparameters(gamma_only_in,ecut_in,poolnproc_in,poolrank_in,distribution_type_in);
    this->nks = nks_in;
    this->kvec_d = kvec_d_in;
    this->gk_ecut = gk_ecut_in;
    for(int ik = 0 ; ik < this->nks ; ++ik)
    {
        kvec_c = 
    }
    return;
}

void PW_Basis_K::setupIndGk()
{
    ModuleBase::TITLE("PW_Basis_K","setupIndGk");
    ModuleBase::timer::tick("PW_Basis_K","setupIndGk");

    igk.create(nks, npw);
    ModuleBase::Memory::record("PW_Basis_K","igk",nks*npw,"int");

    //=============================================
    // Done this in each cpu ( use nks,not nkstot )
    // notice : using cartesian coordinate
    //=============================================
    for (int ik = 0; ik < this->nks; ik++)
    {
        int ng = 0;
        for (int ig = 0; ig < npw ; ig++)
        {
            const double gk2 = this->get_GPlusK_cartesian(ik, ig).norm2();       
            if (gk2 <= this->gk_ecut)
            {
                this->igk(ik, ng) = ig;
                ng++;
            }
        }
        this->npwk[ik] = ng;
        if ( this->npwk_max < ng)
        {
            this->npwk_max = ng;
        }
    }

	// bool out_gk =0; //DIY! mohan 2011-10-03
	// if(out_gk)
	// {
	// 	std::stringstream ss;
	// 	ss << GlobalV::global_out_dir << "PW_GK" << GlobalV::MY_RANK+1 << ".dat";
	// 	std::ofstream ofs( ss.str().c_str() );
	// 	ofs << GlobalC::pw.ggpsi << " (ggpsi, Ry)" << std::endl;
	// 	ofs << GlobalC::pw.ggwfc << " (ggwfc, Ry)" << std::endl;
	// 	ofs << GlobalC::pw.ggwfc2 << " (ggwfc2, Ry)" << std::endl;
	// 	ModuleBase::Vector3<double> f;
	// 	for(int ik=0; ik < nks; ++ik)
	// 	{
	// 		ofs << ik+1 << " (Index of k)" << std::endl;
	// 		ofs << pwb.ngmw << " (Number of plane waves)" << std::endl;	
	// 		for(int ig=0; ig < pwb.ngmw; ++ig)
	// 		{
    //             f = pwb.get_GPlusK_cartesian(ik, ig);
    //             ofs << f.x << " " << f.y << " " << f.z << " " << f.norm() << std::endl;
	// 		}
	// 	}
	// 	ofs.close();
	// }

    ModuleBase::timer::tick("PW_Basis_K","setupIndGk");

    return;
}
