# This is the Makefile of  ABACUS-CMD API

#==========================
# Compiler information 
#==========================
CPLUSPLUS     = icpc
CPLUSPLUS_MPI = mpiicpc
OBJ_DIR = md_obj
NP      = 8

##==========================
## FFTW package needed 
##==========================
#FFTW_DIR = /home/qianrui/intelcompile/impi_fftw
#FFTW_DIR = /home/mohan/1_Software/impi_fftw-3.3.8
HONG_FFTW        = -D__FFTW3
#FFTW_INCLUDE_DIR = ${FFTW_DIR}/include
FFTW_INCLUDE_DIR = ${MKLROOT}/include/fftw
#FFTW_LIB_DIR     = ${FFTW_DIR}/lib
#FFTW_LIB         = -L${FFTW_LIB_DIR} -lfftw3 -Wl,-rpath=${FFTW_LIB_DIR}

##==========================
## FFTW package needed 
##==========================
#BOOST_DIR = /home/mohan/1_Software/impi_boost-1.70.0
BOOST_DIR = /home/qianrui/intelcompile/impi_boost
BOOST_INCLUDE_DIR = ${BOOST_DIR}/include

#==========================
# LIBS and INCLUDES
#==========================
LIBS = -lifcore -lm -lpthread ${FFTW_LIB}

#==========================
# OPTIMIZE OPTIONS
#==========================
INCLUDES = -I. -Icommands -I${FFTW_INCLUDE_DIR} -I${BOOST_INCLUDE_DIR}

# -pedantic turns off more extensions and generates more warnings
# -xHost generates instructions for the highest instruction set available on the compilation host processor
OPTS = ${INCLUDES} -Ofast -std=c++11 -simd -march=native -xHost -m64 -qopenmp -Werror -Wall -pedantic -g

include Makefile.Objects

VPATH=../module_base\
:../module_cell\
:../src_parallel\
:./\
:../\

#==========================
# Define HONG
#==========================
HONG= ${HONG_FFTW} -D__MPI -DMETIS -DMKL_ILP64 -D__NORMAL -D__CMD

FP_OBJS_0=main.o\
$(OBJS_MAIN)\
$(OBJS_BASE)\
$(OBJS_MD)\

FP_OBJS=$(patsubst %.o, ${OBJ_DIR}/%.o, ${FP_OBJS_0})

#==========================
# MAKING OPTIONS
#==========================
cmd : 
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
${OBJ_DIR}/%.o:%.cpp
	${CPLUSPLUS_MPI} ${OPTS} ${OPTS_MPI} -c ${HONG} $< -o $@

.PHONY:clean
clean:
	@ if [ -d $(OBJ_DIR) ]; then rm -rf $(OBJ_DIR); fi
	@ if [ -e ${VERSION}.x ]; then rm -f ${VERSION}.x; fi
