#ifndef EXX_ABFS_UNITTEST_H
#define EXX_ABFS_UNITTEST_H

#include<iostream>
#include<sstream>
using namespace std;

#include "../../../src_ri/exx_abfs.h"

#include "../src_global/matrix-test.h"
#include "../src_ri/exx_abfs-unittest.h"

static void of_abfs_cpp( const string &file_name, const vector<vector<vector<Numerical_Orbital_Lm>>> &orb ) 
{
	for( int i=0; i!=orb.size(); ++i )
		for( int j=0; j!=orb[i].size(); ++j )
			for( int k=0; k!=orb[i][j].size(); ++k )
			{
				stringstream ss;	ss<<file_name<<"_"<<i<<"_"<<j<<"_"<<k;
				ofstream ofs(ss.str().c_str());
				for( const auto & c : orb[i][j][k].psi_uniform )
					ofs<<c<<endl;
				ofs.close();
			}
}

static void cout_matrix_4( const map<size_t,map<size_t,map<size_t,map<size_t,matrix>>>> &matrix )
{
	for( const auto & m1 : matrix )
	{
		const size_t TA = m1.first;
		{
			for( const auto & m2 : m1.second )
			{
				const size_t IA = m2.first; 
				for( const auto & m3 : m2.second )
				{
					const size_t TB = m3.first;
					for( const auto & m4 : m3.second )
					{
						const size_t IB = m4.first;
						const auto & m = m4.second;
						
						cout<<TA<<"\t"<<IA<<"\t"<<TB<<"\t"<<IB<<endl;
						cout<<m<<endl;
					}
				}
			}
		}
	}
}

static void ofs_ms( 
	const string& file_name, 
	const map<size_t,map<size_t,map<size_t,map<size_t,matrix>>>> &ms )
{
	ofstream ofs(file_name.c_str());
	ofs<<"0000"<<endl;
	ofs<<ms.at(0).at(0).at(0).at(0)<<endl;
	ofs<<"0001"<<endl;
	ofs<<ms.at(0).at(0).at(0).at(1)<<endl;
	ofs.close();
}

static void ofs_ms( 
	const string& file_name, 
	const map<size_t,map<size_t,map<size_t,map<size_t,vector<matrix>>>>> &ms )
{
	ofstream ofs(file_name.c_str());
	ofs<<"0000"<<endl;
	for( const auto m : ms.at(0).at(0).at(0).at(0) )
		ofs<<m<<endl;
	ofs<<"0001"<<endl;
	for( const auto m : ms.at(0).at(0).at(0).at(1) )
		ofs<<m<<endl;
	ofs.close();
}

static void ofs_ms( 
	const string& file_name, 
	const map<size_t,map<size_t,map<size_t,map<size_t,vector<vector<matrix>>>>>> &ms )
{
	ofstream ofs(file_name.c_str());
	ofs<<"0000"<<endl;
	for( const auto m1 : ms.at(0).at(0).at(0).at(0) )
		for( const auto m2 : m1 )
			ofs<<m2<<endl;
	ofs<<"0001"<<endl;
	for( const auto m1 : ms.at(0).at(0).at(0).at(1) )
		for( const auto m2 : m1 )
			ofs<<m2<<endl;
	ofs.close();	
}

#endif
