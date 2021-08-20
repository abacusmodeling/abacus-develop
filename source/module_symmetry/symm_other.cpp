#include "symm_other.h"

namespace ModuleSymmetry
{
void Symm_Other::print1(const int &ibrav, const double *cel_const, std::ofstream &ofs_running)
{
	TITLE("Symm_Other","print1");

	OUT(ofs_running,"IBRAV",ibrav);
	if(ibrav==1)
	{
		OUT(ofs_running,"BRAVAIS","SIMPLE CUBIC");
    	OUT(ofs_running,"LATTICE CONSTANT A",cel_const[0]);
	}
	else if(ibrav==2)
	{
		OUT(ofs_running,"BRAVAIS","BODY CENTERED CUBIC");
    	OUT(ofs_running,"LATTICE CONSTANT A",cel_const[0]);
	}
	else if(ibrav==3)
	{
		OUT(ofs_running,"BRAVAIS","FACE CENTERED CUBIC");
    	OUT(ofs_running,"LATTICE CONSTANT A",cel_const[0]);
	}
	else if(ibrav==4)
	{
		OUT(ofs_running,"BRAVAIS","HEXAGONAL CELL");
    	OUT(ofs_running,"LATTICE CONSTANT A",cel_const[0]);
		OUT(ofs_running,"C/A RATIO",cel_const[2]);
	}
	else if(ibrav==5)
	{
		OUT(ofs_running,"BRAVAIS","SIMPLE TETROGONAL CELL");
    	OUT(ofs_running,"LATTICE CONSTANT A",cel_const[0]);
		OUT(ofs_running,"C/A RATIO",cel_const[2]);
	}
	else if(ibrav==6)
	{
		OUT(ofs_running,"BRAVAIS","BODY CENTERED TETROGONAL CELL");
    	OUT(ofs_running,"LATTICE CONSTANT A",cel_const[0]);
		OUT(ofs_running,"C/A RATIO",cel_const[2]);
	}
	else if(ibrav==7)
	{
		OUT(ofs_running,"BRAVAIS","TRIGONAL (RHOMBOEDRIC) CELL");
    	OUT(ofs_running,"LATTICE CONSTANT A",cel_const[0]);
		OUT(ofs_running,"COS(ALPHA)",cel_const[3]);
	}
	else if(ibrav==8)
	{
		OUT(ofs_running,"BRAVAIS","SIMPLE ORTHORHOMBIC CELL");
    	OUT(ofs_running,"LATTICE CONSTANT A",cel_const[0]);
		OUT(ofs_running,"B/A RATIO",cel_const[1]);
		OUT(ofs_running,"C/A RATIO",cel_const[2]);
	}
	else if(ibrav==9)
	{
		OUT(ofs_running,"BRAVAIS","BODY CENTERED ORTHORHOMBIC CELL");
    	OUT(ofs_running,"LATTICE CONSTAddNT A",cel_const[0]);
		OUT(ofs_running,"B/A RATIO",cel_const[1]);
		OUT(ofs_running,"C/A RATIO",cel_const[2]);
	}
	else if(ibrav==10)
	{
		OUT(ofs_running,"BRAVAIS","FACE CENTERED ORTHORHOMBIC CELL");
    	OUT(ofs_running,"LATTICE CONSTAddNT A",cel_const[0]);
		OUT(ofs_running,"B/A RATIO",cel_const[1]);
		OUT(ofs_running,"C/A RATIO",cel_const[2]);
	}
	else if(ibrav==11)
	{
		OUT(ofs_running,"BRAVAIS","BASE CENTERED ORTHORHOMBIC CELL");
    	OUT(ofs_running,"LATTICE CONSTAddNT A",cel_const[0]);
		OUT(ofs_running,"B/A RATIO",cel_const[1]);
		OUT(ofs_running,"C/A RATIO",cel_const[2]);
	}
	else if(ibrav==12)
	{
		OUT(ofs_running,"BRAVAIS","SIMPLE MONOLINIC CELL");
    	OUT(ofs_running,"LATTICE CONSTAddNT A",cel_const[0]);
		OUT(ofs_running,"B/A RATIO",cel_const[1]);
		OUT(ofs_running,"C/A RATIO",cel_const[2]);
		OUT(ofs_running,"COS(BETA)",cel_const[4]);
	}
	else if(ibrav==13)
	{
		OUT(ofs_running,"BRAVAIS","BASED CENTERED MONOLINIC CELL");
    	OUT(ofs_running,"LATTICE CONSTAddNT A",cel_const[0]);
		OUT(ofs_running,"B/A RATIO",cel_const[1]);
		OUT(ofs_running,"C/A RATIO",cel_const[2]);
		OUT(ofs_running,"COS(BETA)",cel_const[4]);	
	}
	else if(ibrav==14)
	{
		OUT(ofs_running,"BRAVAIS","TRICLINIC CELL");
    	OUT(ofs_running,"LATTICE CONSTAddNT A",cel_const[0]);
		OUT(ofs_running,"B/A RATIO",cel_const[1]);
		OUT(ofs_running,"C/A RATIO",cel_const[2]);
		OUT(ofs_running,"COS(ALPHA)",cel_const[3]);	
		OUT(ofs_running,"COS(BETA)",cel_const[4]);	
		OUT(ofs_running,"COS(GAMMA)",cel_const[5]);	
	}
	else
	{
		WARNING_QUIT("Symm_Other::print1","ibrav is wrong.");
	}
	return;
}

bool Symm_Other::right_hand_sense(Vector3<double> &v1,Vector3<double> &v2,Vector3<double> &v3)
{
	double volume = Symm_Other::celvol(v1,v2,v3);
	//OUT(ofs_running,"volume = ",volume);
	if(volume < 0)
	{
		v1.reverse();
		v2.reverse();
		v3.reverse();
		return false;
	}
	return true;
}

//calculate the volume of the cell spanned by the vectors
double Symm_Other::celvol(const Vector3<double> &a, 
const Vector3<double> &b, const Vector3<double> &c)
{
	return a.x * ( b.y * c.z - b.z * c.y ) + a.y * ( b.z * c.x - b.x * c.z ) 
	+ a.z * ( b.x * c.y - b.y * c.x );
}
}

