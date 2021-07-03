#ifndef CONTAINER_OPERATOR_H
#define CONTAINER_OPERATOR_H

#include <cassert>
#include <vector>
#include <map>
using namespace std;



template< typename T>
vector<T> operator + ( const vector<T> & x1, const vector<T> & x2 )
{
	assert(x1.size()==x2.size());
	vector<T> x;
	for( size_t i=0; i!=x1.size(); ++i )
		x.push_back(x1[i]+x2[i]);
	return x;
}

template< typename T>
vector<T> operator - ( const vector<T> & x1, const vector<T> & x2 )
{
	assert(x1.size()==x2.size());
	vector<T> x;
	for( size_t i=0; i!=x1.size(); ++i )
		x.push_back(x1[i]-x2[i]);
	return x;
}

template< typename T1, typename T2 >
map<T1,T2> operator + ( const map<T1,T2> & x1, const map<T1,T2> & x2 )
{
	map<T1,T2> x;
	for( const auto &x1i : x1 )
		x.insert(make_pair( x1i.first, x1i.second + x2.at(x1i.first) ));
	return x;
}

template< typename T1, typename T2 >
map<T1,T2> operator - ( const map<T1,T2> & x1, const map<T1,T2> & x2 )
{
	map<T1,T2> x;
	for( const auto &x1i : x1 )
		x.insert(make_pair( x1i.first, x1i.second - x2.at(x1i.first) ));
	return x;
}



template< typename T1, typename T2 >
vector<T2> operator * ( const T1 & x1, const vector<T2> & x2 )
{
	vector<T2> x;
	for( size_t i=0; i!=x2.size(); ++i )
		x.push_back(x1*x2[i]);
	return x;
}

template< typename T1, typename T21, typename T22 >
map<T21,T22> operator * ( const T1 & x1, const map<T21,T22> & x2 )
{
	map<T21,T22> x;
	for( const auto & x2i : x2 )
		x.insert(make_pair( x2i.first, x1*x2i.second ));
	return x;
}

#endif	// CONTAINER_OPERATOR_H