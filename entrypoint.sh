cmake -B build
cmake --build build -j4
cmake --install build
cd build
ctest --verbose
