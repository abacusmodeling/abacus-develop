#!/bin/bash

GTEST_DIR=/home/qianrui/intelcompile/impi_gtest


GTESTOPTS="-I$GTEST_DIR/include -L$GTEST_DIR/lib -lgtest -lpthread"
make clean > /dev/null 2>&1
mv Makefile Makefile.bak
for((i=0;i<4;++i))
do
head -n 50 Makefile.intel > Makefile
if ((i==0)) ;then
cat >>Makefile<<EOF
HONG = -D__NORMAL
CPLUSPLUS = icpc
GTESTOPTS = $GTESTOPTS
EOF
elif ((i==1)) ;then
cat >>Makefile<<EOF
HONG = -D__MIX_PRECISION -D__NORMAL
TESTFILE0 = \${DOUBLEFILE} \${FLOATFILE}
CPLUSPLUS = icpc
GTESTOPTS = $GTESTOPTS
EOF
elif ((i==2)) ;then
cat >>Makefile<<EOF
HONG = -D__MPI -D__CUDA -D__NORMAL
GTESTOPTS = $GTESTOPTS
EOF
elif ((i==3)) ;then
cat >>Makefile<<EOF
HONG = -D__MPI -D__MIX_PRECISION -D__NORMAL
TESTFILE0 = \${DOUBLEFILE} \${FLOATFILE}
GTESTOPTS = $GTESTOPTS
EOF
fi
tail -n 50 Makefile.intel >>Makefile
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