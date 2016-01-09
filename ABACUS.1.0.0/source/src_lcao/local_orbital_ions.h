#ifndef LOCAL_ORBITAL_IONS_H
#define LOCAL_ORBITAL_IONS_H

#include "local_orbital_elec.h"
#include "../src_ions/ions_move_methods.h"
#include "../src_pw/charge_extra.h"
//#include "../src_develop/src_md/md.h"
//2014-06-06, xiaohui
#include "../src_pw/md.h"
//added by daye 2014/6/19
#include "../src_pw/mdNVT.h"
#include "../src_pw/mdNVE.h"

class Local_Orbital_Ions
{
public:
	Local_Orbital_Ions();
	~Local_Orbital_Ions();

	Local_Orbital_Elec LOE;

	void opt_ions(void);
	//2014-06-06, xiaohui
        mdnvt MDNVT ;
        mdNVE MDNVE ;

private:
	Ions_Move_Methods IMM;
	//MD md;

	bool force_stress(void);
	int istep;

	Charge_Extra CE;
	//choose md ensemble, zheng daye
	int mdtype;

};

#endif
