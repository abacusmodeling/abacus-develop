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

// pblas
void pdtran_(int *m , int *n ,
             double *alpha , double *a , int *ia , int *ja , int *desca ,
             double *beta ,  double *c , int *ic , int *jc , int *descc );
void pdgemm_(char *transa , char *transb , int *m , int *n , int *k ,
             double *alpha , double *a , int *ia , int *ja , int *desca ,
                             double *b , int *ib , int *jb , int *descb ,
             double *beta ,  double *c , int *ic , int *jc , int *descc );
void pdsymm_(char *side , char *uplo , int *m , int *n ,
             double *alpha , double *a , int *ia , int *ja , int *desca ,
                             double *b , int *ib , int *jb , int *descb ,
             double *beta ,  double *c , int *ic , int *jc , int *descc );
void pdtrmm_(char *side , char *uplo , char *transa , char *diag , int *m , int *n ,
             double *alpha , double *a , int *ia , int *ja , int *desca ,
             double *b , int *ib , int *jb , int *descb );

void pzgemm_(char *transa , char *transb , int *m , int *n , int *k ,
        double *alpha , double _Complex *a , int *ia , int *ja , int *desca ,
        double _Complex *b , int *ib , int *jb , int *descb ,
        double *beta ,  double _Complex *c , int *ic , int *jc , int *descc );
/*void pzsymm_(char *side , char *uplo , int *m , int *n ,
 *              double *alpha , double _Complex *a , int *ia , int *ja , int *desca ,
 *                                           double _Complex *b , int *ib , int *jb , int *descb ,
 *                                                        double *beta ,  double _Complex *c , int *ic , int *jc , int *descc );*/
void pztrmm_(char *side , char *uplo , char *transa , char *diag , int *m , int *n ,
        double *alpha , double _Complex *a , int *ia , int *ja , int *desca ,
        double _Complex *b , int *ib , int *jb , int *descb );

