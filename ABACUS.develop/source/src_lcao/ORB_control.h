#ifndef ORB_CONTROL_H 
#define ORB_CONTROL_H 

#include "ORB_gen_tables.h"
#include "ORB_read.h"

class ORB_control 
{

	public:

    ORB_control();
    ~ORB_control();

    // Generate the S(overlap),T,NL matrix.
    void set_orb_tables(
		ORB_gen_tables &OGT, 
		LCAO_Orbitals &orb,
		const double &lat0);

    void clear_after_ions(
		ORB_gen_tables &OGT, 
		LCAO_Orbitals &orb);

};
#endif
