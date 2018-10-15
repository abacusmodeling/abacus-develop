#ifndef LOCAL_ORBTIAL_ELEC
#define LOCAL_ORBITAL_ELEC

#include "../src_pw/tools.h"
#include "../src_pw/threshold_elec.h"

class Local_Orbital_Elec: private Threshold_Elec
{
public:
	Local_Orbital_Elec(){};
	~Local_Orbital_Elec(){};

	void scf(const int &istep);
	void nscf(void);

	static int iter;
	static double avg_iter;
protected:
	static int istep;

private:
	void cal_bands(void);
	void init_mixstep_final_scf(void);


};

#endif
