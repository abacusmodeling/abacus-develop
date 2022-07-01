#!/bin/bash

GTEST_DIR=/home/qianrui/intelcompile/impi_gtest


GTESTOPTS="-I$GTEST_DIR/include -L$GTEST_DIR/lib -lgtest -lpthread"
headn=`grep -n "GTEST needed" Makefile.intel |cut -d ':' -f1`
((headn-=2))
filen=`wc -l Makefile.intel |cut -d ' ' -f1`
tailn=`grep -n "PW_OBJS=" Makefile.intel |cut -d ':' -f1`
((tailn=filen-tailn+2))

make clean > /dev/null 2>&1
mv Makefile Makefile.bak
for((i=0;i<4;++i))
do
head -n $headn Makefile.intel > Makefile
if ((i==0)) ;then
cat >>Makefile<<EOF
HONG = -D__NORMAL
CPLUSPLUS = icpc
GTESTOPTS = $GTESTOPTS
EOF
elif ((i==1)) ;then
cat >>Makefile<<EOF
HONG = -D__MIX_PRECISION -D__NORMAL
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
GTESTOPTS = $GTESTOPTS
EOF
fi
tail -n $tailn Makefile.intel >>Makefile
make > /dev/null

if ((i==0||i==1)) ;then
    ./pw_test.exe
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
make clean > /dev/null 2>&1
done
mv Makefile.bak Makefile
test -e _tmp.txt && rm -f _tmp.txt
exit 0