#include "pdiag_common.h"

#ifdef __MPI
MPI_Comm DIAG_HPSEPS_WORLD;
#endif

void indxg2l(int ia,int ja,int nb,int nprow,int npcol,int *iia,int *jja)
{
    *iia=((ia/nb)/nprow)*nb+ia%nb;
    *jja=((ja/nb)/npcol)*nb+ja%nb;
}

#ifdef __MPI
void indxg2p(MPI_Comm comm2D,int nb,int ia,int ja,int *iarow,int *iacol)
{
    int dim[2],period[2],coord[2];
    MPI_Cart_get(comm2D,2,dim,period,coord);
    *iarow=(ia/nb)%dim[0];
    *iacol=(ja/nb)%dim[1];
    return;
}
#endif

#ifdef __MPI
void mpi_sub_row(MPI_Comm comm2D, MPI_Comm *comm_row)
{
    int blongs[2]={0,1};
    MPI_Cart_sub(comm2D,blongs,comm_row);
    return;
}
#endif

#ifdef __MPI
void mpi_sub_col(MPI_Comm comm2D, MPI_Comm *comm_col)
{
    int blongs[2]={1,0};
    MPI_Cart_sub(comm2D,blongs,comm_col);
    return;
}
#endif

