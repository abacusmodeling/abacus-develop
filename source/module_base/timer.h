//==========================================================
// AUTHOR : fangwei , mohan
// DATE : 2008-11-06
// UPDATE : Peize Lin at 2019-11-21
//==========================================================
#ifndef TIMER_H
#define TIMER_H

#include <ctime>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <map>
namespace ModuleBase
{

//==========================================================
// CLASS :
// NAME : timer(calculate time)
//==========================================================
class timer
{
	public:
	
	struct Timer_One
	{
		double cpu_start;
		double cpu_second = 0.0;
		size_t calls = 0;
		size_t order = n_now++;
		bool start_flag = true;
	};
	
	static std::map<std::string,std::map<std::string,Timer_One>> timer_pool;
//==========================================================
// MEMBER FUNCTIONS :
// NAME : tick(use twice at a time)
// NAME : start
// NAME : finish
// NAME : enable
// NAME : disable
// NAME : print_all
//==========================================================
	static void tick(const std::string &class_name_in,const std::string &name_in);

	static void start(void);
	static void finish(std::ofstream &ofs,const bool print_flag = 1);

	static void enable(void){ disabled = false; }
	static void disable(void){ disabled = true; }

	static void print_all(std::ofstream &ofs);
	static long double print_until_now(void);

	private:

//==========================================================
// MEMBER VARIABLES : 
// NAME : disabled(if disabled , timer can't work)
// NAME : n_now(the index of clocks) 
//==========================================================
	static bool disabled;
	static size_t n_now;

//==========================================================
// MEMBER FUNCTIONS :
// NAME : cpu_time(calculate time)
//==========================================================
	static double cpu_time(void);

};

}
#endif

