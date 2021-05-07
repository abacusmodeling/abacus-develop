# This is the Makefile of ABACUS-ORB API

CPLUSPLUS      = icpc
CPLUSPLUS_MPI = icpc 
FFTW_DIR = /home/mohan/1_Software/impi_fftw-3.3.8
OBJ_DIR = obj
NP      = 1

HONG_FFTW        = -D__FFTW3
FFTW_INCLUDE_DIR = ${FFTW_DIR}/include
FFTW_LIB_DIR     = ${FFTW_DIR}/lib
FFTW_LIB         = -L${FFTW_LIB_DIR} -lfftw3 -Wl,-rpath=${FFTW_LIB_DIR}

#==========================
# LIBS and INCLUDES
#==========================
LIBS = -lifcore -lm -lpthread ${FFTW_LIB}

#==========================
# OPTIMIZE OPTIONS
#==========================
INCLUDES = -I. -Icommands -I${FFTW_INCLUDE_DIR}

OPTS     = ${INCLUDES} -Ofast -traceback -std=c++11 -simd -march=native -xHost -m64 -qopenmp -Werror -Wall -pedantic -g
#OPTS_MPI = -cxx=${CPLUSPLUS}



include Makefile.Objects

VPATH=../../src_global\
:../../src_lcao\
:./\

#==========================
# Define HONG
#==========================
HONG= -DMETIS -DMKL_ILP64 -D__NORMAL ${HONG_FFTW} 

#==========================
# OPTIMIZE OPTIONS
#==========================
OPTS_GDB = -g -W

#==========================
# OBJECTS NEEDED
#==========================
#FP_OBJS_0=$(OBJS_ORBITAL)\
#$(OBJS_GLOBAL)\
#main.o\

FP_OBJS_0=main.o\
$(OBJS_TRY)\
$(OBJS_GLOBAL)\
$(OBJS_ORBITAL)\

FP_OBJS=$(patsubst %.o, ${OBJ_DIR}/%.o, ${FP_OBJS_0})
PDIAG_OBJS=$(patsubst %.o, ${OBJ_DIR}/%.o, ${OBJS_PDIAG})
PDIAG_MR=$(patsubst %.o, ${OBJ_DIR}/%.o, ${PDIAG_MR_0})

#==========================
# MAKING OPTIONS
#==========================
fp_mpi : 
	@ make init
	@ make -j $(NP) serial2 

init :
	@ if [ ! -d $(OBJ_DIR) ]; then mkdir $(OBJ_DIR); fi
	@ if [ ! -d $(OBJ_DIR)/README ]; then echo "This directory contains all of the .o files" > $(OBJ_DIR)/README; fi
	@ if [ ! -d ../bin ]; then mkdir ../bin; fi

serial : ${FP_OBJS} ${HEADERS} 
	${CPLUSPLUS} ${OPTS} $(FP_OBJS) ${LIBS} -o ${VERSION}.x 

serial2 : ${FP_OBJS} 
	${CPLUSPLUS} ${OPTS} $(FP_OBJS) ${LIBS} -o ${VERSION}.x 

#==========================
# rules
#==========================
${OBJ_DIR}/%.o:%.cpp
	${CPLUSPLUS_MPI} ${OPTS} ${OPTS_MPI} -c ${HONG} $< -o $@
${OBJ_DIR}/%.o:%.f
	${FORTRAN} -c ${HONG} $< -o $@	 

.PHONY:clean
clean:
	@ if [ -d $(OBJ_DIR) ]; then rm -rf $(OBJ_DIR); fi
