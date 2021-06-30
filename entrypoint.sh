cd ABACUS.develop
cmake -B build cmake
cmake --build build -j4
cmake --install build
cd tests
ABACUS_PATH=ABACUS NP=2 bash Autotest.sh
bash clean.sh
