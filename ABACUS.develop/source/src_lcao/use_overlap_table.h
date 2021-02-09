//=========================================================
//AUTHOR : Mohan 
//DATE : 2009-04-22
//=========================================================
#ifndef USE_OVERLAP_TABLE_H
#define USE_OVERLAP_TABLE_H

#include "../src_pw/tools.h"
#include "make_overlap_table.h"
#include "make_gaunt_table.h"
#include "ylm.h"

class Use_Overlap_Table
{
	public:

	friend class Hamilt_Linear;
	
	Use_Overlap_Table();
	~Use_Overlap_Table();

	void gen_tables( const int &job0 );
	void set_unit( const double &v ){lat0=v;}
	
	void snap_psipsi(
		double olm[],
		const int &job, // 0 for matrix element of either S or T, 1 for its derivatives
	    const char &dtype, // derivative type, 'S' for overlap, 'T' for kinetic energy
		const Vector3<double> &R1,
    	const int &I1,
    	const int &l1,
    	const int &m1,
    	const int &n1,
    	const Vector3<double> &R2,
    	const int &I2,
    	const int &l2,
    	const int &m2,
		const int &n2,
		complex<double> *olm1=NULL)const;
		
	//job = 0 for vnl matrix elements
	//job = 1 for its derivatives
	void snap_psibeta(
		double nlm[],
		const int& job,
		const Vector3<double> &R1,
		const int &I1,
		const int &l1,
		const int &m1,
		const int &n1,
		const Vector3<double> &R2,
		const int &I2,
		const int &l2,
		const int &m2,
		const int &n2,
		const Vector3<double> &Rnl,
		const int &type,
		complex<double> *nlm1=NULL,
		const int is=0)const;


	// set as public because in hamilt_linear, 
	// we need to destroy the tables: SR,TR,NR
	// after ionic optimization is done.
	Make_Overlap_Table MOT;

	// if we want to add table for descriptors,
	// we should consider here -- mohan 2021-02-09

	private:

	Make_Gaunt_Table MGT;

	double get_distance(const Vector3<double> &R1, const Vector3<double> &R2)const;

	double lat0;

};

extern Use_Overlap_Table UOT;

#endif
