#include "pw_basis.h"

PW_Basis::PW_Basis()
{
    ig2fft = NULL;
    is2ir = NULL;
    gdirect = NULL;		
    gcar = NULL; 
    gg = NULL;
    startz = NULL;
    numz = NULL;  
}
PW_Basis:: ~PW_Basis()
{
    if(ig2fft != NULL) delete[] ig2fft;
    if(is2ir != NULL) delete[] is2ir;
    if(gdirect != NULL) delete[] gdirect;
    if(gcar != NULL) delete[] gcar;
    if(gg != NULL) delete[] gg;
}