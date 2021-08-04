// =============================================================================
//                          C++ Header File
// Project:         
// File:            export.h 
// Author:          mohan
// Comment:         
// Warning:         
// Start time:      2008-9-3
// Last modified:   
// =============================================================================

#include "export.h"

void IF_MATCH(const string &name,const string &name2)
{
	if(name!=name2)
	{
		if(GlobalV::MY_RANK == 0)
		{
			std::cout<<"\n Can not match : "<<name<<"  "<<name2<<std::endl;
		}
#ifdef __MPI
		MPI_Finalize();
#endif
		exit(0);
	}
	
//	std::cout<<setw(12)<<name<<std::endl;
	return;
}



