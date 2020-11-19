# Makefile for plane wave basis code.
include ../Makefile.system
include ../Makefile.Objects

VPATH=../src_tools\
:../src_algorithms\
:../src_parallel\
:../src_wannier\
:../src_dc\

#==========================
#Define HONG
#==========================
#(1)no mpi
HONG = -D__FFTW2 
#(2)mpi
#HONG = -D__FFTW2 -D__MPI 
#(3)memory
#HONG = -D__FFTW2 -D_MCD_CHECK -DWIN32 -DMCD_VERBOSE

#==========================
#OPTIMIZE OPTIONS
#==========================
#OPTS = -O3 -w -fpermissive
#OPTS= -v
OPTS= -O3
#(1)gdb
#OPTS = -g -W -O3 -w -fpermissive
#(2)gprof
#POTS = -pg

#==========================
# OBJECTS NEEDED
#==========================
PW_OBJS=$(OBJS_OTHER)\
$(OBJS_PARALLEL)\
$(OBJS_TOOLS) \
$(OBJS_FIRST_PRINCIPLE)\
$(OBJS_COMMON)\
$(OBJS_MLWF)\
$(OBJS_EXX)\
main.o \

#==========================
# MAKING OPTIONS
#==========================
pw_mpi : ${PW_OBJS} ${HEADERS}
	${CPLUSPLUS} ${OPTS} $(PW_OBJS) ${LIBS} ${LIBS_MPI} -o ${VERSION}.pw_mpi.x

pw : ${PW_OBJS} ${HEADERS}
	${CPLUSPLUS} ${OPTS} $(PW_OBJS) ${LIBS} -o ${VERSION}.pw.x

pw_gdb : ${PW_OBJS} ${HEADERS}
	${CPLUSPLUS} ${OPTS} $(PW_OBJS) ${LIBS} -o ${VERSION}.pw_gdb.x

pw_gprof : ${PW_OBJS} ${HEADERS}
	${CPLUSPLUS} $(OBJS) $(PW_OBJS) ${LIBS} -o ${VERSION}.pw_gprof.x

pw_mem : ${PW_OBJS} ./src_parallel/mcd.o ${HEADERS}
	${CPLUSPLUS} $(OBJS) ./src_parallel/mcd.o ${LIBS} -o ${VERSION}.pw_mem.x

#==========================
# rules
#==========================
.cpp.o:
	${CPLUSPLUS} ${OPTS} -c ${HONG}$ $< -o $@

.PHONY:clean
clean:
	$(RM) *.o *.x
