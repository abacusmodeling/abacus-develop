#pragma once

#include<iostream>
#include<fstream>
using namespace std;

template<typename T>
static void print_tmp(const T &S)
{
	for( const auto &s1 : S)
	{
		const auto i1=s1.first;
		for( const auto &s2 : s1.second )
		{
			const auto i2=s2.first;
			for( const auto &s3 : s2.second )
			{
				const auto i3=s3.first;
				cout<<i1<<"\t"<<i2<<"\t"<<i3<<endl;
			}
		}
	}
}

void print_matrixess(const map<size_t,map<size_t,map<Abfs::Vector3_Order<double>,weak_ptr<matrix>>>> &ms, ofstream &ofs)
{
	for( const auto & c1 : ms )
		for( const auto & c2 : c1.second )
			for( const auto & c3 : c2.second )
			{
				ofs<<c1.first<<"\t"<<c2.first<<"\t"<<c3.first<<endl;
				ofs<<*c3.second.lock()<<endl;
			}
}

void print_matrixess(const map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,shared_ptr<matrix>>>> &ms, ofstream &ofs)
{
	for( const auto & c1 : ms )
		for( const auto & c2 : c1.second )
			for( const auto & c3 : c2.second )
			{
				ofs<<c1.first<<"\t"<<c2.first<<"\t"<<c3.first<<endl;
				ofs<<*c3.second<<endl;
			}
}