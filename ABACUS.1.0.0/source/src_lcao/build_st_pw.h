#ifndef BUILD_ST_PW_H
#define BUILD_ST_PW_H

#include "../src_pw/tools.h"

class Build_ST_pw
{
	public:
	Build_ST_pw();
	~Build_ST_pw();

	void set_ST(const int &ik, const char& dtype);
	void set_local(const int &ik);

};

#endif
