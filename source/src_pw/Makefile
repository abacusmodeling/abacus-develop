# This is the Makefile of ABACUS-PW API

#==========================
# Compiler information 
#==========================
CPLUSPLUS     = icpc
CPLUSPLUS_MPI = mpiicpc 
FFTW_DIR = /home/qianrui/intelcompile/impi_fftw
OBJ_DIR = obj
NP      = 1

#==========================
# FFTW package needed 
#==========================
HONG_FFTW        = -D__FFTW3
FFTW_INCLUDE_DIR = ${FFTW_DIR}/include
FFTW_LIB_DIR     = ${FFTW_DIR}/lib
FFTW_LIB         = -L${FFTW_LIB_DIR} -lfftw3 -Wl,-rpath=${FFTW_LIB_DIR}

#==========================
# LAPACK package needed 
#==========================
LAPACK_DIR    		= $(MKLROOT)
LAPACK_INCLUDE_DIR	= ${LAPACK_DIR}/include 
LAPACK_LIB_DIR		= ${LAPACK_DIR}/lib/intel64
LAPACK_LIB 			= -L${LAPACK_LIB_DIR} -Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread\
					  -lmkl_core -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -Wl,--end-group -Wl,-rpath=${LAPACK_LIB_DIR}

#==========================
# LIBS and INCLUDES
#==========================
LIBS = -lifcore -lm -lpthread ${FFTW_LIB} ${LAPACK_LIB}

#==========================
# OPTIMIZE OPTIONS
#==========================
INCLUDES = -I. -Icommands -I${FFTW_INCLUDE_DIR} -I${LAPACK_INCLUDE_DIR}

# -pedantic turns off more extensions and generates more warnings
# -xHost generates instructions for the highest instruction set available on the compilation host processor
OPTS = ${INCLUDES} -Ofast -std=c++11 -simd -march=native -xHost -m64 -qopenmp -Werror -Wall -pedantic -g

include Makefile.Objects

VPATH=../src_global\
:../module_base\
:../src_parallel\
:../module_cell\
:../src_io\
:../src_ions\
:../\

#==========================
# Define HONG
#==========================
#HONG= -DMETIS -DMKL_ILP64 -D__NORMAL ${HONG_FFTW} 
HONG= ${HONG_FFTW}  -D__MPI
FP_OBJS_0=main.o\
$(OBJS_BASE)\
$(OBJS_PW)\
$(OBJS_PARALLEL)\
${OBJS_IO}\
${OBJS_CELL}\
$(OBJS_MAIN)

FP_OBJS=$(patsubst %.o, ${OBJ_DIR}/%.o, ${FP_OBJS_0})

#==========================
# MAKING OPTIONS
#==========================
ORB : 
	@ make init
	@ make -j $(NP) parallel

init :
	@ if [ ! -d $(OBJ_DIR) ]; then mkdir $(OBJ_DIR); fi
	@ if [ ! -d $(OBJ_DIR)/README ]; then echo "This directory contains all of the .o files" > $(OBJ_DIR)/README; fi

parallel : ${FP_OBJS} 
	${CPLUSPLUS_MPI} ${OPTS} $(FP_OBJS) ${LIBS} -o ${VERSION}.x 

#==========================
# rules
#==========================
${OBJ_DIR}/main.o:../main.cpp
	${CPLUSPLUS_MPI} ${OPTS} ${OPTS_MPI} -c ${HONG} ../main.cpp -o ${OBJ_DIR}/main.o
${OBJ_DIR}/%.o:%.cpp
	${CPLUSPLUS_MPI} ${OPTS} ${OPTS_MPI} -c ${HONG} $< -o $@

.PHONY:clean
clean:
	@ if [ -d $(OBJ_DIR) ]; then rm -rf $(OBJ_DIR); fi
	@ if [ -e ${VERSION}.x ]; then rm -f ${VERSION}.x; fi
