#ifndef EXX_ABFS_PARALLEL_H
#define EXX_ABFS_PARALLEL_H

#include "exx_abfs.h"

class Exx_Abfs::Parallel
{
	public:
	
	class Distribute
	{
		public:
		class Htime;
		class Kmeans;
		class Order;
	};
	
	class Communicate
	{
		public:
		class Allreduce;
		class Hexx;
		#if EXX_DM==1
		class DM;
		#elif EXX_DM==2
		class DM2;
		#endif
	};
};

#endif