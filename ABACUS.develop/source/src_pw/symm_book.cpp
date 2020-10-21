#include "tools.h"
#include "symm_book.h"

Symm_Book::Symm_Book()
{

}


Symm_Book::~Symm_Book()
{

}


void Symm_Book::Bravais()
{
	TITLE("Symm_Book","Bravais");

	write("There are only 7 crystal systems");
	write("Any other point groups are the subgroup of the 7");
	write("The 7 are: S2, C2h, D2h, D4h, D3d, D6h, Oh");

	write("There are 32 point groups:");
	ofs_running << setw(5) << " S2" << setw(2) << ":" <<
	setw(5) << "S2" << 
	setw(5) << "C1" << endl; 

	ofs_running << setw(5) << " C2h" << setw(2) << ":" <<
	setw(5) << "C2h" <<
	setw(5) << "C2" <<
	setw(5) << "C1h" << endl; 

	ofs_running << setw(5) << " D2h" << setw(2) << ":" <<
	setw(5) << "D2h" <<
	setw(5) << "D2" <<
	setw(5) << "C2v" << endl; 

	ofs_running << setw(5) << " D4h" << setw(2) << ":" <<
	setw(5) << "D4h" <<
	setw(5) << "C4" <<
	setw(5) << "S4" <<
	setw(5) << "D4" <<
	setw(5) << "C4v" <<
	setw(5) << "C4h" <<
	setw(5) << "D2d" << endl; 

	ofs_running << setw(5) << " D3d" << setw(2) << ":" <<
	setw(5) << "D3d" <<
	setw(5) << "S6" <<
	setw(5) << "C3" <<
	setw(5) << "C3v" <<
	setw(5) << "D3" << endl;

	ofs_running << setw(5) << " D6h" << setw(2) << ":" <<
	setw(5) << "D6h" <<
	setw(5) << "C6" <<
	setw(5) << "C3h" <<
	setw(5) << "C6h" <<
	setw(5) << "C6v" <<
	setw(5) << "D6" <<
	setw(5) << "D3h" << endl;

	ofs_running << setw(5) << " Oh" << setw(2) << ":" <<
	setw(5) << "Oh" <<
	setw(5) << "T" <<
	setw(5) << "O" <<
	setw(5) << "Th" <<
	setw(5) << "Td" << endl;

	write("There are two types of point groups");
	write("1: C1.C2.C3.C4.C6.D2.D3.D4.D6.T.O");
	write("2: S2.S4.S6.C1h.C2h.C3h.C4h.C6h.C2v.C3v.C4v.C6v.D2h.D3h.D4h.D6h.D2d.D3d.Th.Oh.Td");

	write("The 14 Bravais Lattices are");
	write("01. Cubic P (a=b=c) (sc)");
	write("02. Cubic I (a=b=c) (bcc)");
	write("03. Cubic F (a=b=c) (fcc)");
	write("04. Trigonal and hexagonal P");
	write("05. Tetragonal P");
	write("06. Tetragonal I");
	write("07. Trigonal R");
	write("08. Orthorhombic P (a!=b!=c)");
	write("09. Orthorhombic I");
	write("10. Orthorhombic F");
	write("11. Orthorhombic C");
	write("12. Monoclinic P");
	write("13. Monoclinic A");
	write("14. Triclinic P");

	write("The rank of the point groups");
	write("48: Oh");
	write("24: D6h. Th. O. Td");
	write("16: D4h");
	write("12: D3h. C6h. D3d. C6v. D6. T");
	write("08: D2h. C4h. D2d. C4v. D4");
	write("06: C6.  C3h. S6.  C3v. D3");
	write("04: C4.  C2h. S4.  C2v. D2");
	write("03: C3");
	write("02: C2,  C1h, S2");
	write("01: C1");

}

void Symm_Book::write(const string &name)
{
	ofs_running << " " << name << "." << endl;
}
