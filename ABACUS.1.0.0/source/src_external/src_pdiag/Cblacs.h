//------------------------------------>8======================================
//  Copyright (c) 2016, Yu Shen (shenyu@ustc.edu.cn)
//  All rights reserved.
//  
//  Redistribution and use in source and binary forms, with or without
//  modification, are permitted provided that the following conditions are met:
//      * Redistributions of source code must retain the above copyright
//        notice, this list of conditions and the following disclaimer.
//      * Redistributions in binary form must reproduce the above copyright
//        notice, this list of conditions and the following disclaimer in the
//        documentation and/or other materials provided with the distribution.
//      * Neither the name of the <organization> nor the
//        names of its contributors may be used to endorse or promote products
//        derived from this software without specific prior written permission.
//  
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT Yu Shen SHALL BE LIABLE FOR ANY
//  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
//  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
//  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//====================================8<----------------------------------------
// blacs
    // Initialization
void Cblacs_pinfo(int *myid, int *nprocs);
void Cblacs_get(int icontxt, int what, int *val);
void Cblacs_gridinit(int* icontxt, char *layout, int nprow, int npcol);
void Cblacs_gridmap(int* icontxt, int *usermap, int ldumap, int nprow, int npcol);
    // Destruction
void Cblacs_gridexit(int icontxt);
    // Informational and Miscellaneous
void Cblacs_gridinfo(int icontxt, int* nprow, int *npcol, int *myprow, int *mypcol);
int Cblacs_pnum(int icontxt, int prow, int pcol);
void Cblacs_pcoord(int icontxt, int pnum, int *prow, int *pcol);
void Cblacs_barrier(int icontxt, char *scope);
    // Point to  Point
void Cdgesd2d(int icontxt, int m, int n, double *a, int lda, int rdest, int cdest);
void Cdgerv2d(int icontxt, int m, int n, double *a, int lda, int rsrc, int csrc);
    // Combine
//void Cdgamx2d(int icontxt, int scope, int top, int m, int n, 
//              double *a, int lda, int *ra, int *ca, int rcflag, int rdest, int cdest);
