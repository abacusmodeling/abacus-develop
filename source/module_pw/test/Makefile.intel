
#==========================
# Compiler information 
#==========================
CPLUSPLUS = mpiicpc
CUDA_COMPILE = nvcc
OBJ_DIR = pw_obj
NP      = 12
#==========================
# Objects
#==========================
VPATH=../../src_parallel\
:../../module_base\
:../

PW_OBJS_0=matrix.o\
matrix3.o\
tool_quit.o\
mymath3.o\
timer.o\
global_variable.o\
parallel_global.o\
pw_basis.o\
pw_distributer.o\
pw_init.o\
pw_transform.o\
pw_distributeg.o\
pw_distributeg_method1.o\
pw_distributeg_method2.o\
fft.o\
pw_basis_k.o\
pw_operation.o\
pw_transform_k.o

DOUBLEFILE=test1-2.o\
test1-3.o\
test1-4.o\
test2-2.o\
test2-3.o

FLOATFILE=test1-2f.o\
test1-3f.o\
test1-4f.o

TESTFILE0 = ${DOUBLEFILE}

#==========================
# Options
#==========================
#No MPI
# HONG = -D__NORMAL
# CPLUSPLUS = icpc

#Mix Precision
# HONG = -D__MIX_PRECISION -D__NORMAL
# TESTFILE0 = ${DOUBLEFILE} ${FLOATFILE}
# CPLUSPLUS = icpc

#Only MPI
HONG = -D__MPI -D__NORMAL

#MPI + Mix Precision
# HONG = -D__MPI -D__MIX_PRECISION -D__NORMAL
# TESTFILE0 = ${DOUBLEFILE} ${FLOATFILE}

#Cuda
#HONG = -D__MPI -D__CUDA -D__NORMAL

#Cuda & Mix Precision
#HONG = -D__MPI -D__CUDA -D__MIX_PRECISION -D__NORMAL
# TESTFILE0 = ${DOUBLEFILE} ${FLOATFILE}


PW_OBJS=$(patsubst %.o, ${OBJ_DIR}/%.o, ${PW_OBJS_0})
TESTFILE=$(patsubst %.o, ${OBJ_DIR}/%.o, ${TESTFILE0})

##==========================
## FFTW package needed 
##==========================
#Use fftw package
#  FFTW_DIR = /home/qianrui/intelcompile/impi_fftw
#  FFTW_LIB_DIR     = ${FFTW_DIR}/lib
#  FFTW_INCLUDE_DIR = ${FFTW_DIR}/include
#  FFTW_LIB         = -L${FFTW_LIB_DIR} -lfftw3 -Wl,-rpath=${FFTW_LIB_DIR}
# FFTW_LIB         = -L${FFTW_LIB_DIR} -lfftw3 -lfftw3f -Wl,-rpath=${FFTW_LIB_DIR}

#Use mkl_fftw
FFTW_INCLUDE_DIR = ${MKLROOT}/include/fftw
FFTW_LIB = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core


##==========================
## CUDA needed 
##==========================
# CUDA_DIR = /usr/local/cuda-11.0
# CUDA_INCLUDE_DIR	= ${CUDA_DIR}/include 
# CUDA_LIB_DIR		= ${CUDA_DIR}/lib64
# CUDA_LIB			= -L${CUDA_LIB_DIR} -lcufft -lcublas -lcudart

##==========================
## GTEST needed 
##==========================
GTEST_DIR = /home/qianrui/intelcompile/impi_gtest
GTESTOPTS = -I${GTEST_DIR}/include -L${GTEST_DIR}/lib -lgtest -lpthread


#LIBS = ${FFTW_LIB} ${CUDA_LIB} -ltcmalloc -lprofiler
LIBS = ${FFTW_LIB} ${CUDA_LIB}
OPTS = -I${FFTW_INCLUDE_DIR} ${HONG} -Ofast -std=c++11 -simd -m64 -qopenmp -Wall -pedantic -g

#==========================
# MAKING OPTIONS
#==========================
pw : 
	@ make init
	@ make -j $(NP) utest.exe

init :
	@ if [ ! -d $(OBJ_DIR) ]; then mkdir $(OBJ_DIR); fi

utest.exe: ${PW_OBJS} ${TESTFILE}
	${CPLUSPLUS} ${OPTS} ${TESTFILE} utest.cpp test_tool.cpp ${PW_OBJS}  ${LIBS} -o utest.exe ${GTESTOPTS}
${OBJ_DIR}/%.o:%.cpp
	${CPLUSPLUS} ${OPTS} -c ${HONG} $< -o $@ ${GTESTOPTS}

.PHONY:clean
clean:
	@ if [ -d $(OBJ_DIR) ]; then rm -rf $(OBJ_DIR); fi
	@ if [ -e utest.exe ]; then rm -f utest.exe; fi
