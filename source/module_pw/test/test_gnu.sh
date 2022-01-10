#!/bin/bash

GTEST_DIR=/home/qianrui/gnucompile/g_gtest
FFTW_DIR=/home/qianrui/gnucompile/g_fftw-3.3.8-mpi

GTESTOPTS="-I$GTEST_DIR/include -L$GTEST_DIR/lib -lgtest -lpthread"

make clean > /dev/null 2>&1
mv Makefile Makefile.bak
for((i=0;i<4;++i))
do
head -n 50 Makefile.gnu > Makefile
if ((i==0)) ;then
cat >>Makefile<<EOF
HONG = -D__NORMAL
CPLUSPLUS = icpc
EOF
elif ((i==1)) ;then
cat >>Makefile<<EOF
HONG = -D__MIX_PRECISION -D__NORMAL
TESTFILE0 = \${DOUBLEFILE} \${FLOATFILE}
CPLUSPLUS = icpc
EOF
elif ((i==2)) ;then
cat >>Makefile<<EOF
HONG = -D__MPI -D__CUDA -D__NORMAL
EOF
elif ((i==3)) ;then
cat >>Makefile<<EOF
HONG = -D__MPI -D__MIX_PRECISION -D__NORMAL
TESTFILE0 = \${DOUBLEFILE} \${FLOATFILE}
EOF
fi
cat >>Makefile<<EOF
GTESTOPTS = $GTESTOPTS
FFTW_DIR = ${FFTW_DIR}
FFTW_LIB_DIR     = \${FFTW_DIR}/lib
FFTW_INCLUDE_DIR = \${FFTW_DIR}/include
FFTW_LIB         = -L\${FFTW_LIB_DIR} -lfftw3 -lfftw3f -Wl,-rpath=\${FFTW_LIB_DIR}
EOF
tail -n 36 Makefile.gnu >>Makefile
make > /dev/null 2>&1

if ((i==0||i==1)) ;then
    ./pw_test.exe
else
    ./pw_test.exe
    sleep 1
    echo "====================="
    mpirun -np 2 ./pw_test.exe |grep PASSED
    sleep 1
    echo "====================="
    mpirun -np 8 ./pw_test.exe |grep PASSED
fi
make clean > /dev/null 2>&1
done
mv Makefile.bak Makefile

exit 0