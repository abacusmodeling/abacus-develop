#ifndef BERRYPHASE
#define BERRYPHASE


#include "../src_pw/unkOverlap_pw.h"
#include "../src_lcao/unkOverlap_lcao.h"



class berryphase
{
public:
	unkOverlap_pw pw_method;
	unkOverlap_lcao lcao_method;
	int total_string;
	vector<vector<int>> k_index;
	int nppstr;
	int direction;
	int occ_nbands;
	int GDIR;
	
	berryphase();
	~berryphase();
	void get_occupation_bands();
	void lcao_init();
	void set_kpoints(const int direction);
	double stringPhase(int index_str, int nbands);
	void Berry_Phase(int nbands, double &pdl_elec_tot, int &mod_elec_tot);
	void Macroscopic_polarization();
	
};

#endif
