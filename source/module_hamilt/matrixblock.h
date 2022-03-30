#ifndef MATRIXBLOCK
#define MATRIXBLOCK

template<typename T>
struct MatrixBlock
{
    /* this is a simple template block of a matrix
       would change to Eigen in the future */
    T* p;
    size_t row;
    size_t col;
};

#endif