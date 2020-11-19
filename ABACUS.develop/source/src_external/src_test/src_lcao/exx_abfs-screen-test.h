#ifndef EXX_ABFS_SCREEN_TEST_H
#define EXX_ABFS_SCREEN_TEST_H

#include <map>
#include <string>
#include "src_lcao/abfs.h"

static void test_screen( const string & file_name, const map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,double>>> & m )
{
	ofstream ofs(file_name);
	for( const auto m1 : m )
		for( const auto m2 : m1.second )
			for( const auto m3 : m2.second )
				ofs<<m1.first<<"\t"<<m2.first<<"\t"<<m3.first<<"\t"<<m3.second<<endl;
	ofs.close();
}

#endif