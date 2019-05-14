# Makefile for plane wave basis code.
include Makefile.system
include Makefile.Objects

VPATH=./src_global\
:./src_parallel\
:./src_external/src_pdiag\
:./src_pw\
:./src_lcao\
:./src_ions\
:./src_external/src_pdiag/MRRR\
:./\

#==========================
# Define HONG
#==========================
#(1)no mpi
HONG_SER = -D__FP ${HONG_FFTW} 
HONG_SER_SELINV = -D__FP ${HONG_FFTW} -D__SELINV
HONG_GDB = -g -D__FP ${HONG_FFTW}
#(2)mpi
HONG_MPI = -D__FP ${HONG_FFTW} -D__MPI 
HONG_MPI_SELINV = -D__FP ${HONG_FFTW} -D__MPI -D__SELINV -DMETIS -DMKL_ILP64 -DEXX_DM=2 -DTEST_LIBXC=0 -DTEST_EXX_LCAO=0 -DTEST_EXX_RADIAL=1
#(3)memory
HONG_MEM = ${HONG_FFTW} -D__FP -D_MCD_CHECK -DWIN32 -DMCD_VERBOSE

#!!!!!!!!!!!!!!!!!!!! CHANE HERE IF YOU LIKE !!!!!!!!!!!!!!
#! change series version or parallel version~~~
HONG=${HONG_MPI_SELINV}
#HONG=${HONG_SER_SELINV}
#!!!!!!!!!!!!!!!!!!!! CHANE HERE IF YOU LIKE !!!!!!!!!!!!!!

#==========================
# OPTIMIZE OPTIONS
#==========================
#(1)gdb
#xiaohui modify 2013-12-27
OPTS_GDB = -g -W -I${BOOST_INCLUDE_DIR}
#(2)gprof
#POTS = -pg

#==========================
# OBJECTS NEEDED
#==========================
FP_OBJS=$(OBJS_OTHER)\
$(OBJS_FIRST_PRINCIPLE)\
$(OBJS_PARALLEL)\
$(OBJS_TOOLS) \
$(OBJS_ONSCALING)\
$(OBJS_SELINV)\
$(OBJS_CG)\
$(OBJS_WANNIER)\
$(OBJS_COMMON)\
$(OBJS_MLWF)\
$(OBJS_EXX)\
main.o \

PDIAG_OBJS=$(OBJS_PDIAG)\

#==========================
# MAKING OPTIONS
#==========================
fp_mpi : ${FP_OBJS} ${PDIAG_MR} ${PDIAG_OBJS} ${HEADERS}
	${CPLUSPLUS_MPI} ${OPTS} ${OPTS_MPI} $(FP_OBJS) ${PDIAG_OBJS} ${PDIAG_MR} ${LIBS} -o  ../bin/${VERSION}.mpi.1.0.0_exx

fp : ${FP_OBJS} ${HEADERS} 
	${CPLUSPLUS} ${OPTS} $(FP_OBJS) ${LIBS} -o ${VERSION}.fp.x 

fp_gdb : ${FP_OBJS} ${HEADERS}
	${CPLUSPLUS} ${OPTS_GDB} $(FP_OBJS) ${LIBS} -o ${VERSION}.fp_gdb.x

fp_gprof : ${FP_OBJS} ${HEADERS} 
	${CPLUSPLUS} ${OPTS} $(FP_OBJS) ${LIBS} -o ${VERSION}.fp_gprof.x 

fp_mem : ${FP_OBJS} ../src_parallel/mcd.o ${HEADERS}
	${CPLUSPLUS} $(FP_OBJS) ../src_parallel/mcd.o ${LIBS} -o ${VERSION}.fp_mem.x

#==========================
# rules
#==========================
.cpp.o:
	${CPLUSPLUS_MPI} ${OPTS} ${OPTS_MPI} -c ${HONG} $< -o $@

.PHONY:clean
clean:
	$(RM) *.o
