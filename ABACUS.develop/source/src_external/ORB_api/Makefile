# This is the Makefile of ABACUS

include Makefile.system
include Makefile.Objects

VPATH=../../src_global\
:../../src_lcao\
:./\

#==========================
# Define HONG
#==========================
HONG=-D__MPI -DMETIS -DMKL_ILP64

#==========================
# OPTIMIZE OPTIONS
#==========================
OPTS_GDB = -g -W 

#==========================
# OBJECTS NEEDED
#==========================
FP_OBJS_0=$(OBJS_ORBITAL)\
$(OBJS_GLOBAL)\
main.o\

FP_OBJS=$(patsubst %.o, ${OBJ_DIR}/%.o, ${FP_OBJS_0})
PDIAG_OBJS=$(patsubst %.o, ${OBJ_DIR}/%.o, ${OBJS_PDIAG})
PDIAG_MR=$(patsubst %.o, ${OBJ_DIR}/%.o, ${PDIAG_MR_0})

#==========================
# MAKING OPTIONS
#==========================
fp_mpi : 
	@ make init
	@ make -j $(NP) parallel

init :
	@ if [ ! -d $(OBJ_DIR) ]; then mkdir $(OBJ_DIR); fi
	@ if [ ! -d $(OBJ_DIR)/README ]; then echo "This directory contains all of the .o files" > $(OBJ_DIR)/README; fi
	@ if [ ! -d ../bin ]; then mkdir ../bin; fi

parallel : ${FP_OBJS} ${PDIAG_OBJS} ${PDIAG_MR} ${HEADERS}
	${CPLUSPLUS_MPI} ${OPTS} ${OPTS_MPI} $(FP_OBJS) ${PDIAG_OBJS} ${PDIAG_MR} ${LIBS} -o  ../bin/${VERSION}.mpi

serial : ${FP_OBJS} ${HEADERS} 
	${CPLUSPLUS} ${OPTS} $(FP_OBJS) ${LIBS} -o ${VERSION}.fp.x 

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
	@ if [ -d ../bin ]; then rm -rf ../bin; fi
