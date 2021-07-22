#ifndef ABFS_TEMPLATE_H
#define ABFS_TEMPLATE_H

#include "abfs.h"

template<typename T1,typename T2,typename T3,typename Tmatrix>
void Abfs::delete_empty_ptrs( map<T1,map<T2,map<T3,weak_ptr<Tmatrix>>>> &ptrs )
{
	TITLE("Abfs","delete_empty_ptrs");
	for( auto iter1=ptrs.begin(); iter1!=ptrs.end(); )
	{
		for( auto iter2=iter1->second.begin(); iter2!=iter1->second.end(); )
		{
			for( auto iter3=iter2->second.begin(); iter3!=iter2->second.end(); )
			{
				if( iter3->second.expired() )	iter2->second.erase(iter3++);
				else	++iter3;
			}
			if( iter2->second.empty() )	iter1->second.erase(iter2++);
			else	++iter2;
		}
		if( iter1->second.empty() )	ptrs.erase(iter1++);
		else	++iter1;
	}
}

template<typename Tkey,typename Tmatrix>
void Abfs::delete_threshold_ptrs( map<Tkey,Tmatrix> &ptrs, const double threshold )
{
	for( auto iter=ptrs.begin(); iter!=ptrs.end(); )
	{
		if( iter->second.absmax()<=threshold )	ptrs.erase(iter++);
		else	++iter;
	}
}

template<typename Tkey,typename Tmatrix>
void Abfs::delete_threshold_ptrs( map<Tkey,shared_ptr<Tmatrix>> &ptrs, const double threshold )
{
	for( auto iter=ptrs.begin(); iter!=ptrs.end(); )
	{
		if( iter->second->absmax()<=threshold )	ptrs.erase(iter++);
		else	++iter;
	}
}


template<typename Tkey1,typename Tkey2,typename Tvalue>
void Abfs::delete_threshold_ptrs( map<Tkey1,map<Tkey2,Tvalue>> &ptrs, const double threshold )
{
	for( auto iter=ptrs.begin(); iter!=ptrs.end(); )
	{
		delete_threshold_ptrs(iter->second, threshold);
		if( iter->second.empty() )	ptrs.erase(iter++);
		else	++iter;
	}
}

/*
template<typename T1,typename T2,typename T3,typename Tmatrix>
void Abfs::delete_threshold_ptrs( map<T1,map<T2,map<T3,shared_ptr<Tmatrix>>>> &ptrs, const double threshold)
{
	TITLE("Abfs","delete_threshold_ptrs");
	for( auto iter1=ptrs.begin(); iter1!=ptrs.end(); )
	{
		for( auto iter2=iter1->second.begin(); iter2!=iter1->second.end(); )
		{
			for( auto iter3=iter2->second.begin(); iter3!=iter2->second.end(); )
			{
				if( iter3->second->absmax()<=threshold )	iter2->second.erase(iter3++);
				else	++iter3;
			}
			if( iter2->second.empty() )	iter1->second.erase(iter2++);
			else	++iter2;
		}
		if( iter1->second.empty() )	ptrs.erase(iter1++);
		else	++iter1;
	}
}

template<typename T1,typename T2,typename T3,typename Tmatrix>
void Abfs::delete_threshold_ptrs( map<T1,map<T2,map<T3,Tmatrix>>> &ptrs, const double threshold)
{
	TITLE("Abfs","delete_threshold_ptrs");
	for( auto iter1=ptrs.begin(); iter1!=ptrs.end(); )
	{
		for( auto iter2=iter1->second.begin(); iter2!=iter1->second.end(); )
		{
			for( auto iter3=iter2->second.begin(); iter3!=iter2->second.end(); )
			{
				if( iter3->second.absmax()<=threshold )	iter2->second.erase(iter3++);
				else	++iter3;
			}
			if( iter2->second.empty() )	iter1->second.erase(iter2++);
			else	++iter2;
		}
		if( iter1->second.empty() )	ptrs.erase(iter1++);
		else	++iter1;
	}
}
*/

template<typename T1, typename T2, typename Tother>
vector<pair<T1,T2>> Abfs::get_atom_pair(const map<T1,map<T2,Tother>> &m)
{
	vector<pair<T1,T2>> atom_pairs;
	for(const auto &mA : m)
		for(const auto &mB : mA.second)
			atom_pairs.push_back({mA.first,mB.first});
	return atom_pairs;
}

#endif