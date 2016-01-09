//==========================================================
// AUTHOR : sltk, ywcui, mohan
// START DATE : 2007-04-11
// LAST DATE : 2008-11-22
//==========================================================
#ifndef GRID_H
#define GRID_H

#include <stdexcept>
#include <boost/functional/hash/hash.hpp>
#include "util.h"
#include "sltk_atom.h"
#include "sltk_atom_input.h"
#include "../src_pw/tools.h"
#include "../src_pw/pw_basis.h"
extern PW_Basis pw;

//==========================================================
// STRUCT :
// NAME : AtomLink
// NAME : CellSet
//==========================================================

struct AtomLink
{
	FAtom fatom;
	AtomLink* next_p;

	// Constructors and destructor
	AtomLink
	(
	    const FAtom& atom = FAtom(),
	    AtomLink* const pointNext = NULL //mohan fix bug 2011/09/26, from NullPtr->NULL
	);

};

struct CellSet
{
	AtomLink* address;
	int length;
	int in_grid[3];
	CellSet();
};

//==========================================================
// CLASS NAME :
// Atom_input : defined elsewhere
//==========================================================

class Atom_input;

//==========================================================
// CLASS NAME :
// Grid :
//==========================================================

class Grid
{
public:

	// Constructors and destructor
	Grid();
	~Grid();

	void init(const Atom_input &input);
	//2015-05-07
	void delete_vector(const Atom_input &input);


	//Static data
	static const double TOLERATE_ERROR;
	static const boost::hash<int> INT_HASHER;
	static const char* const ERROR[3];

	//Data
	int natom;// Total atoms.
	bool pbc; // periodic boundary condition
	bool expand_flag;
	double sradius;// searching radius
	double d_minX;// origin of all cells
	double d_minY;
	double d_minZ;
	int dx;
	int dy;
	int dz;
	int layer;
	double cell_x_length;
	double cell_y_length;
	double cell_z_length;
	CellSet ***Cell; //dx , dy ,dz is cell number in each direction,respectly.

	double grid_length[3];
	double vec1[3];
	double vec2[3];
	double vec3[3];
	double lat_now;
	bool init_cell_flag;

	// Inner Function
protected:
	AtomLink* const getHashCode(const FAtom &atom)const;
	AtomLink* atomlink;
	AtomLink* cordon_p;// Warning! A guard! Don't delete it!

private:

//==========================================================
// MEMBER FUNCTIONS :
// Three Main Steps:
// NAME : setMemberVariables (read in datas from Atom_input,
// 			init cells.)
// NAME : setAtomLinkArray( set the AtomLinkArray twice,
// 			first use Build_Hash,second use Fold_Hash)
// NAME : setBoundaryAdjacent( Consider different situations,
// 			if not_expand case : nature/periodic boundary
// 			condition , if expand_case)
//==========================================================
	void setMemberVariables(const Atom_input &input);
	void setAtomLinkArray(const Atom_input &input);
	void setBoundaryAdjacent(const Atom_input &input);

//==========================================================
//
//==========================================================
	AtomLink* const Build_Cache(const Atom_input &input) throw(std::out_of_range, std::logic_error);
	bool Push(const FAtom& atom);
	void In_Which_Cell(int &a, int &b, int &c, const FAtom &atom)const;
	void Build_Cell(void);
	void Build_Hash_Table(AtomLink* const pointCache);
	void Fold_Hash_Table(void) throw(std::logic_error);
	static int Hash_one_hit;

//==========================================================
//
//==========================================================
	void Construct_Adjacent_expand(const int i, const int j, const int k);
	void Construct_Adjacent_expand_periodic(
	    const int i, const int j, const int k, const int ia);

	void Construct_Adjacent_begin(void);
	void Construct_Adjacent_nature(
	    const int i, const int j, const int k, const int ia);
	void Construct_Adjacent_periodic(
	    const int i, const int j, const int k, const int ia);
	void Construct_Adjacent_final(
	    const int i, const int j, const int k, const int ia,
	    const int i2, const int j2, const int k2, const int ia2);
};

#endif // GRID_H
