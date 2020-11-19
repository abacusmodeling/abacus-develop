#ifndef ATOM_ARRANGE_H
#define ATOM_ARRANGE_H

#include "../src_pw/tools.h"
#include "sltk_grid.h"
#include "sltk_atom_input.h"

class atom_arrange
{
public:

	atom_arrange();
	~atom_arrange();
	
	static void search( const double &search_radius_bohr);
	static void set_sr_OV(void);
	static void set_sr_NL(void);
	//2015-05-07
	static void delete_vector(const double &search_radius_bohr);

private:

};

#endif
