#include "pw_basis.h"

PW_Basis::PW_Basis()
{
    ig2isz = NULL;
    istot2ixy = NULL;  
    is2ixy = NULL;
    ixy2ip = NULL; 
    startis = NULL;
    nst_per = NULL;
    gdirect = NULL;		
    gcar = NULL; 
    gg = NULL;  
}
PW_Basis:: ~PW_Basis()
{
    if(ig2isz != NULL) delete[] ig2isz;
    if(istot2ixy != NULL) delete[] istot2ixy;
    if(is2ixy != NULL) delete[] is2ixy;
    if(ixy2ip != NULL) delete[] ixy2ip;
    if(startis != NULL) delete[] startis;
    if(nst_per != NULL) delete[] nst_per;
    if(gdirect != NULL) delete[] gdirect;
    if(gcar != NULL) delete[] gcar;
    if(gg != NULL) delete[] gg;
}