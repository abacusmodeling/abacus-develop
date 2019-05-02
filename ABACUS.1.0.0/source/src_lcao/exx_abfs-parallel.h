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
	};
	
	class Communicate
	{
		public:
		class Allreduce;
		class Hexx;
		class DM;
		class DM2;
	};
};

#endif