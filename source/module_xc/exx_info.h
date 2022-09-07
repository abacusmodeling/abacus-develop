#ifndef EXX_INFO_H
#define EXX_INFO_H

#include "xc_functional.h"

struct Exx_Info
{
	enum class Hybrid_Type {No,HF,PBE0,SCAN0,HSE,Generate_Matrix};

	struct Exx_Info_Global
	{
		Exx_Info::Hybrid_Type hybrid_type;

		double hybrid_alpha = 0.25;
		double hse_omega = 0.11;
		
		bool separate_loop = true;
		size_t hybrid_step = 1;
	};	
	Exx_Info_Global info_global;



	struct Exx_Info_Lip
	{
		const Exx_Info::Hybrid_Type &hybrid_type;

		const double &hse_omega;
		double lambda;

		Exx_Info_Lip( const Exx_Info::Exx_Info_Global &info_global )
			:hybrid_type(info_global.hybrid_type),
			 hse_omega(info_global.hse_omega){}
	};
	Exx_Info_Lip info_lip;	



	struct Exx_Info_RI
	{
		const Exx_Info::Hybrid_Type &hybrid_type;
		
		const double &hse_omega;
		
		double pca_threshold = 0;
		std::vector<std::string> files_abfs;
		double C_threshold  = 0;
		double V_threshold  = 0;
		double dm_threshold = 0;
		double cauchy_threshold = 0;
		double ccp_threshold = 0;
		double ccp_rmesh_times = 10;
		double kmesh_times = 4;

		Exx_Info_RI( const Exx_Info::Exx_Info_Global &info_global )
			:hybrid_type(info_global.hybrid_type),
			 hse_omega(info_global.hse_omega){}
	};
	Exx_Info_RI info_ri;


	Exx_Info()
		:info_lip(this->info_global),
		 info_ri(this->info_global){}
};

#endif
