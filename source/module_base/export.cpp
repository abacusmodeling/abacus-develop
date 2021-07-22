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
		if(MY_RANK == 0)
		{
			cout<<"\n Can not match : "<<name<<"  "<<name2<<endl;
		}
#ifdef __MPI
		MPI_Finalize();
#endif
		exit(0);
	}
	
//	cout<<setw(12)<<name<<endl;
	return;
}



