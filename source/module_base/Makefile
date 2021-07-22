# This is the Makefile of ABACUS-base API

#==========================
# Compiler information 
#==========================
CPLUSPLUS     = mpiicpc
CPLUSPLUS_MPI = mpiicpc 
FFTW_DIR = /home/lrx/ABACUS/fftw3.3.9/fftw
LAPACK_DIR = /opt/intel/compilers_and_libraries_2017.8.262/linux/mkl
HONG_LAPACK        = -D__MKL -DMKL_ILP64
LAPACK_INCLUDE_DIR = ${LAPACK_DIR}/include
LAPACK_LIB_DIR     = ${LAPACK_DIR}/lib/intel64
LAPACK_LIB         = -L${LAPACK_LIB_DIR} -Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -Wl,--end-group -Wl,-rpath=${LAPACK_LIB_DIR}
OBJ_DIR = obj
NP      = 4

#==========================
# LIBS and INCLUDES
#==========================
LIBS = -lifcore -lm -lpthread -lfftw3 -L${FFTW_DIR}/lib ${LAPACK_LIB}

#==========================
# OPTIMIZE OPTIONS
#==========================
INCLUDES = -I. -Icommands -I${FFTW_DIR}/include -I${LAPACK_INCLUDE_DIR}

# -pedantic turns off more extensions and generates more warnings
# -xHost generates instructions for the highest instruction set available on the compilation host processor
OPTS = ${INCLUDES} -Ofast -std=c++11 -simd -march=native -xHost -m64 -qopenmp -Werror -Wall -pedantic -g

include Makefile.Objects

VPATH=../src_global\
:./\

#==========================
# Define HONG
#==========================
HONG= -DMETIS -DMKL_ILP64 -D__NORMAL ${HONG_LAPACK}

FP_OBJS_0=main.o\
$(OBJS_BASE)\

FP_OBJS=$(patsubst %.o, ${OBJ_DIR}/%.o, ${FP_OBJS_0})

#==========================
# MAKING OPTIONS
#==========================
base : 
	@ make init
	@ make -j $(NP) serial

init :
	@ if [ ! -d $(OBJ_DIR) ]; then mkdir $(OBJ_DIR); fi
	@ if [ ! -d $(OBJ_DIR)/README ]; then echo "This directory contains all of the .o files" > $(OBJ_DIR)/README; fi

serial : ${FP_OBJS} 
	${CPLUSPLUS} ${OPTS} $(FP_OBJS) ${LIBS} -o ${VERSION}.x 

#==========================
# rules
#==========================
${OBJ_DIR}/%.o:%.cpp
	${CPLUSPLUS_MPI} ${OPTS} ${OPTS_MPI} -c ${HONG} $< -o $@

.PHONY:clean
clean:
	@ if [ -d $(OBJ_DIR) ]; then rm -rf $(OBJ_DIR); fi
