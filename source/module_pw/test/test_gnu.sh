#!/bin/bash

GTEST_DIR=/home/qianrui/gnucompile/g_gtest
FFTW_DIR=/home/qianrui/gnucompile/g_fftw-3.3.8-mpi

GTESTOPTS="-I$GTEST_DIR/include -L$GTEST_DIR/lib -lgtest -lpthread"
headn=`grep -n "FFTW package needed" Makefile.gnu |cut -d ':' -f1`
((headn-=2))
filen=`wc -l Makefile.gnu |cut -d ' ' -f1`
tailn=`grep -n "PW_OBJS=" Makefile.gnu |cut -d ':' -f1`
((tailn=filen-tailn+2))
make clean > /dev/null 2>&1
mv Makefile Makefile.bak
for((i=0;i<4;++i))
do
head -n $headn Makefile.gnu > Makefile
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
tail -n $tailn Makefile.gnu >>Makefile
make > /dev/null 2>&1

if ((i==0)) ;then
    ./pw_test.exe
elif ((i==1)) ;then
    ./pw_test.exe
    echo "valgrind test:(1 processors)"
    valgrind ./pw_test.exe >_tmp.txt 2>&1
    cat _tmp.txt|egrep "(ERROR SUMMARY)|(lost)";
else
    echo "1 processor:"
    ./pw_test.exe
    sleep 1
    echo "2 processors:"
    mpirun -np 2 ./pw_test.exe >_tmp.txt 2>&1
    cat _tmp.txt|grep PASSED
    cat _tmp.txt|grep FAILED
    sleep 1
    echo "8 processors:"
    mpirun -np 8 ./pw_test.exe >_tmp.txt 2>&1
    cat _tmp.txt|grep PASSED
    cat _tmp.txt|grep FAILED
fi
#for mix compile
if ((i==3)) ;then
    echo "valgrind test:(1 processors)"
    valgrind ./pw_test.exe >_tmp.txt 2>&1
    cat _tmp.txt|egrep "(ERROR SUMMARY)|(lost)";
    echo "valgrind test:(2 processors)"
    mpirun -np 2 valgrind ./pw_test.exe >_tmp.txt 2>&1
    cat _tmp.txt|egrep "(ERROR SUMMARY)|(lost)";
    echo "valgrind test:(8 processors)"
    mpirun -np 8 valgrind ./pw_test.exe >_tmp.txt 2>&1
    cat _tmp.txt|egrep "(ERROR SUMMARY)|(lost)";
fi
make clean > /dev/null 2>&1
done
mv Makefile.bak Makefile
test -e _tmp.txt && rm -f _tmp.txt

exit 0