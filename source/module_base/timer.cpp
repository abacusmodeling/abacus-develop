//==========================================================
// AUTHOR : fangwei , mohan
// DATE : 2008-11-06,
// UPDATE : Peize Lin at 2019-11-21
//==========================================================
#include "timer.h"
#include "chrono"
#include<vector>

#ifdef __MPI
#include "mpi.h"
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

namespace ModuleBase
{

//----------------------------------------------------------
// EXPLAIN :
//----------------------------------------------------------
bool timer::disabled = false;
size_t timer::n_now = 0;
std::map<std::string,std::map<std::string,timer::Timer_One>> timer::timer_pool;

void timer::finish(std::ofstream &ofs,const bool print_flag)
{
	timer::tick("","total");
	if(print_flag)
		print_all( ofs );
}

//----------------------------------------------------------
//
//----------------------------------------------------------
void timer::start(void)
{
	// first init ,then we can use tick
	timer::tick("","total");
}

double timer::cpu_time(void)
{
//----------------------------------------------------------
// EXPLAIN : here static is important !!
// only first call can let t0 = 0,clock begin
// when enter this function second time , t0 > 0
//----------------------------------------------------------
	// static clock_t t0 = clock();
	// const clock_t t1 = clock() - t0;
	// return (t1<0) ? 0 : (double)t1/CLOCKS_PER_SEC;

	// static time_t t0 = time(NULL);
	// const time_t t1 = time(NULL);
	// double res = difftime(t1, t0);
	// return (res<0) ? 0 : res;
	static auto t1 = std::chrono::system_clock::now();
	const auto t2 = std::chrono::system_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1);
	return double(duration.count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den;
		// mohan add, abandon the cross point time 2^32 ~ -2^32 .
}

void timer::tick(const std::string &class_name,const std::string &name)
{
//----------------------------------------------------------
// EXPLAIN : if timer is disabled , return
//----------------------------------------------------------
	if (disabled)
		return;

#ifdef _OPENMP
	if(!omp_get_thread_num())
#endif
	{
		Timer_One &timer_one = timer_pool[class_name][name];

//----------------------------------------------------------
// CALL MEMBER FUNCTION :
// NAME : cpu_time
//
// EXPLAIN :
// if start_flag == true,means a new clock counting begin,
// hence we record the start time of this clock counting.
// if start_flag == false, means it's the end of this counting,
// so we add the time during this two 'time point'  to the clock time storage.
//----------------------------------------------------------
		if(timer_one.start_flag)
		{
#ifdef __MPI
			timer_one.cpu_start = MPI_Wtime();
#else
			timer_one.cpu_start = cpu_time();
#endif
			++timer_one.calls;
			timer_one.start_flag = false;
		}
		else
		{
#ifdef __MPI
			timer_one.cpu_second += MPI_Wtime() - timer_one.cpu_start;
#else
			// if(class_name=="electrons"&&name=="c_bands")
			// {
			// 	cout<<"call times"<<timer_one.calls<<endl;
			// 	cout<<"electrons c_bands cost time:"<<endl;
			// 	cout<<cpu_time()<<"-"<<timer_one.cpu_start<<endl;
			// }

			timer_one.cpu_second += (cpu_time() - timer_one.cpu_start);
#endif
			timer_one.start_flag = true;
		}
	} // end if(!omp_get_thread_num())
}

long double timer::print_until_now(void)
{
	// stop the clock
	timer::tick("","total");
	// start again
	timer::tick("","total");
	return timer_pool[""]["total"].cpu_second;
}

void timer::print_all(std::ofstream &ofs)
{
	constexpr double small = 0.1; // cpu = 10^6
	// if want to print > 1s , set small = 10^6

	std::vector<std::pair<std::pair<std::string,std::string>,Timer_One>> timer_pool_order;
	for(auto &timer_pool_A : timer_pool)
	{
		const std::string class_name = timer_pool_A.first;
		for(auto &timer_pool_B : timer_pool_A.second)
		{
			const std::string name = timer_pool_B.first;
			const Timer_One timer_one = timer_pool_B.second;
			if(timer_pool_order.size() < timer_one.order+1)
				timer_pool_order.resize(timer_one.order+1);
			//timer_pool_order[timer_one.order] = {{class_name, name}, timer_one}; //qianrui change it to make it compatible with old compiler version
			timer_pool_order[timer_one.order] = std::pair<std::pair<std::string,std::string>, Timer_One> {std::pair<std::string,std::string >{class_name,name}, timer_one};
		}
	}

	std::cout << std::setprecision(2);
	ofs << std::setprecision(3);
	std::cout<<"\n  |CLASS_NAME---------|NAME---------------|TIME(Sec)-----|CALLS----|AVG------|PER%-------" << std::endl;
	ofs <<"\n\n\n\n  |CLASS_NAME---------|NAME---------------|TIME(Sec)-----|CALLS----|AVG------|PER%-------" << std::endl;
	for(auto &timer_pool_order_A : timer_pool_order)
	{
		const std::string &class_name = timer_pool_order_A.first.first;
		const std::string &name = timer_pool_order_A.first.second;
		const Timer_One &timer_one = timer_pool_order_A.second;

		if(timer_one.cpu_second < small)
			continue;

		ofs << std::resetiosflags(std::ios::scientific);
		ofs  << " "
			// << std::setw(2)  << timer_one.level
			 << std::setw(2)  << " "
			 << std::setw(20) << class_name
			 << std::setw(20) << name
			 << std::setw(15) << std::setprecision(5) << timer_one.cpu_second
			 << std::setw(10) << timer_one.calls
			 << std::setw(10) << std::setprecision(2) << timer_one.cpu_second/timer_one.calls
			 << std::setw(10) << timer_one.cpu_second / timer_pool_order[0].second.cpu_second * 100 << "%" << std::endl;

		std::cout << std::resetiosflags(std::ios::scientific);
		
		std::cout << " " 
			// << std::setw(2)  << timer_one.level
			 << std::setw(2)  << " "
			 << std::setw(20) << class_name
			 << std::setw(20) << name
			 << std::setw(15) << std::setprecision(5) << timer_one.cpu_second
			 << std::setw(10) << timer_one.calls
			 << std::setw(10) << std::setprecision(2) << timer_one.cpu_second/timer_one.calls
			 << std::setw(10) << timer_one.cpu_second / timer_pool_order[0].second.cpu_second * 100 << "%" << std::endl;
	}
	std::cout<<" ----------------------------------------------------------------------------------------"<<std::endl;
	ofs <<" ----------------------------------------------------------------------------------------"<<std::endl;
}
}

/*
void timer::print_all(std::ofstream &ofs)
{
//	std::cout<<"\n timer::print_all()"<<std::endl;
	const double small = 0.1; // cpu = 10^6
	// if want to print > 1s , set small = 10^6

	std::cout << std::setprecision(2);

	// prepare
	bool *print_flag = new bool[n_clock];
	for(int i=0; i<n_clock; i++)
	{
		print_flag[i] = false;
	}

	int type = 1; // 2:calls 1:total_time
	bool non_reorder = 1;

	std::cout<<"\n  |CLASS_NAME---------|NAME---------------|TIME(Sec)-----|CALLS----|AVG------|PER%-------" << std::endl;
	ofs <<"\n\n\n\n  |CLASS_NAME---------|NAME---------------|TIME(Sec)-----|CALLS----|AVG------|PER%-------" << std::endl;
	ofs << std::setprecision(3);
	for (int i=0; i<n_clock; i++)
	{
		int k = 0;
		double tmp = -1.0;

		if(non_reorder)
		{
			k = i;
		}
		else
		{
			// search in all clocks
			for(int j=0; j<n_clock; j++)
			{
				if(print_flag[j])
				{
					continue;
				}
				if(type==1)
				{
					if(tmp < cpu_second[j])
					{
						k = j;
						tmp = cpu_second[j];
					}
				}
				else if(type==2)
				{
					if(tmp < calls[j])
					{
						k = j;
						tmp = calls[j];
					}
				}
			}
		}
		print_flag[k]=true;

		if ((cpu_second[k] >= 0 && cpu_second[k] < small) ||
		        (cpu_second[k] <= 0 && cpu_second[k] > -small))
		{
			continue;
		}

		if( level[k] > 'X' ) continue;


		const long double spend_time = cpu_second[k];
		const double average_spend_time = spend_time/calls[k];


		ofs  << " "
			 << std::setw(2) << level[k]
			 << std::setw(20) << class_name[k]
			 << std::setw(20) << name[k]
			 << std::setw(15) << spend_time
			 << std::setw(10) << calls[k]
			 << std::setw(10) << std::setprecision(2) << average_spend_time
			 << std::setw(10) << spend_time / cpu_second[0] * 100 << "%" << std::endl;


		std::cout << std::resetiosflags(ios::scientific);

		std::cout  << " "
		     << std::setw(2) << level[k]
			 << std::setw(20) << class_name[k]
			 << std::setw(20) << name[k]
			 << std::setw(15) << spend_time
			 << std::setw(10) << calls[k]
			 << std::setw(10) << std::setprecision(2) << average_spend_time
			 << std::setw(10) << spend_time / cpu_second[0] * 100 << "%" << std::endl;

	}
	std::cout<<" ----------------------------------------------------------------------------------------"<<std::endl;
	ofs <<" ----------------------------------------------------------------------------------------"<<std::endl;
	delete[] print_flag;
	return;
}
*/
