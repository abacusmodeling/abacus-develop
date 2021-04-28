#include "../src_io/output.h"
#include "pseudo_us.h"

pseudo_us::pseudo_us()
{
}

pseudo_us::~pseudo_us()
{
}

//---------------------------------------------------------------------
void pseudo_us::set_pseudo_us(const Pseudopot_upf &upf)
{
	this->set_pseudo_nc(upf);
} // end subroutine set_pseudo_upf

void pseudo_us::print_pseudo_us(ofstream &ofs)
{
	ofs << "\n\n internal format for pseudopotential ";
	output out;
	print_pseudo_nc(ofs);
	ofs << "\n pseudo_us : ";
	ofs << "\n nqf	" << nqf
	<< "\n nqlc	" << nqlc;
	out.printrm(ofs, " qqq : ", qqq);
	out.printr3_d(ofs, " qfunc : ", qfunc);
	out.printr4_d(ofs, " qfcoef : ", qfcoef);
	ofs << "\n ----------------------\n";
}


