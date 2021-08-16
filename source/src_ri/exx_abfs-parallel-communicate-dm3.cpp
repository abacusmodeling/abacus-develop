#include "exx_abfs-parallel-communicate-dm3.h"
#include "exx_abfs-parallel-communicate-dm3-template.h"
#include "../src_pw/global.h"
#include "../src_lcao/global_fp.h"
#include "abfs-template.h"

#include "../src_external/src_test/test_function.h"
#include "../src_external/src_test/src_global/complexmatrix-test.h"

/*
template<> vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>>
Exx_Abfs::Parallel::Communicate::DM3::K_to_R(const vector<matrix> &DK_2D, const double threshold_D) const
{

	TITLE("Exx_Abfs::Parallel::Communicate::DM3::K_to_R");
	assert(DK_2D.size()==GlobalV::NSPIN);
	const double SPIN_multiple = 0.5*GlobalV::NSPIN;
	
	vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> DR_a2D(GlobalV::NSPIN);
	for(int is=0; is!=GlobalV::NSPIN; ++is)
	{
		for(int iwt1_local=0; iwt1_local!=DK_2D[is].nr; ++iwt1_local)
		{
			const int iwt1 = GlobalC::ParaO.MatrixInfo.col_set[iwt1_local];
			const int iat1 = GlobalC::ucell.iwt2iat[iwt1];
			const int iw1 = GlobalC::ucell.iwt2iw[iwt1];
			for(int iwt2_local=0; iwt2_local!=DK_2D[is].nc; ++iwt2_local)
			{
				const int iwt2 = GlobalC::ParaO.MatrixInfo.row_set[iwt2_local];
				const int iat2 = GlobalC::ucell.iwt2iat[iwt2];
				const int iw2 = GlobalC::ucell.iwt2iw[iwt2];
				const double dm = DK_2D[is](iwt1_local,iwt2_local);
				if(abs(dm) > threshold_D)
				{
					matrix &DR_a2D_box2 = DR_a2D[is][iat1][iat2][{0,0,0}];
					if(!DR_a2D_box2.c)
						DR_a2D_box2.create(GlobalC::ucell.atoms[GlobalC::ucell.iat2it[iat1]].nw, GlobalC::ucell.atoms[GlobalC::ucell.iat2it[iat2]].nw);
					DR_a2D_box2(iw1,iw2) = DK_2D[is](iwt1_local,iwt2_local) * SPIN_multiple;
				}
			}
		}
	}
	return DR_a2D;
}
*/

/*
template<> vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>>
Exx_Abfs::Parallel::Communicate::DM3::K_to_R(const vector<ComplexMatrix> &DK_2D, const double threshold_D) const
{
	{
		static int istep=0;
		ofstream ofs("DK_2D_"+TO_STRING(istep++)+"_"+TO_STRING(GlobalV::MY_RANK));
		ofs<<DK_2D<<endl;
	}
	
	TITLE("Exx_Abfs::Parallel::Communicate::DM3::K_to_R");
	const double SPIN_multiple = 0.5*GlobalV::NSPIN;
	const Abfs::Vector3_Order<int> Born_von_Karman_period = Vector3<int>{GlobalC::kv.nmp[0],GlobalC::kv.nmp[1],GlobalC::kv.nmp[2]};
	const vector<Abfs::Vector3_Order<int>> supercell_boxes = Abfs::get_Born_von_Karmen_boxes(Born_von_Karman_period);
	
	vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> DR_a2D(GlobalV::NSPIN);
	for(int ik=0; ik!=DK_2D.size(); ++ik)
	{
		for(int iwt1_local=0; iwt1_local!=DK_2D[ik].nr; ++iwt1_local)
		{
			const int iwt1 = GlobalC::ParaO.MatrixInfo.col_set[iwt1_local];
			const int iat1 = GlobalC::ucell.iwt2iat[iwt1];
			const int iw1 = GlobalC::ucell.iwt2iw[iwt1];
			for(int iwt2_local=0; iwt2_local!=DK_2D[ik].nc; ++iwt2_local)
			{
				const int iwt2 = GlobalC::ParaO.MatrixInfo.row_set[iwt2_local];
				const int iat2 = GlobalC::ucell.iwt2iat[iwt2];
				const int iw2 = GlobalC::ucell.iwt2iw[iwt2];
				for(const Abfs::Vector3_Order<int> &box2 : supercell_boxes)
				{
					matrix &DR_a2D_box2 = DR_a2D[GlobalC::kv.isk[ik]][iat1][iat2][box2];
					if(!DR_a2D_box2.c)
						DR_a2D_box2.create(GlobalC::ucell.atoms[GlobalC::ucell.iat2it[iat1]].nw, GlobalC::ucell.atoms[GlobalC::ucell.iat2it[iat2]].nw);
					DR_a2D_box2(iw1,iw2) += real( DK_2D[ik](iwt1_local,iwt2_local) * exp(-TWO_PI*IMAG_UNIT*(GlobalC::kv.kvec_c[ik]*(box2*GlobalC::ucell.latvec))) ) * SPIN_multiple;
				}
			}
		}
	}
	for(auto &DR_a2D_is : DR_a2D)
		Abfs::delete_threshold_ptrs(DR_a2D_is,threshold_D);
	
	{
		static int istep=0;
		ofstream ofs("DR_a2D_"+TO_STRING(istep++)+"_"+TO_STRING(GlobalV::MY_RANK));
		ofs<<DR_a2D<<endl;
	}	
	return DR_a2D;
}
*/


void Exx_Abfs::Parallel::Communicate::DM3::cal_DM(const double threshold_D)
{
	TITLE("Exx_Abfs::Parallel::Communicate::DM3::cal_DM");
	vector<map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,matrix>>>> DR_a2D = GlobalV::GAMMA_ONLY_LOCAL
		? K_to_R(GlobalC::LOC.wfc_dm_2d.dm_gamma, threshold_D)
		: K_to_R(GlobalC::LOC.wfc_dm_2d.dm_k, threshold_D);
	DMr = allreduce.a2D_to_exx(DR_a2D);

	/*{
		static int istep=0;
		ofstream ofs("DMr_"+TO_STRING(istep++)+"_"+TO_STRING(GlobalV::MY_RANK));
		ofs<<DMr<<endl;
	}*/	
}
