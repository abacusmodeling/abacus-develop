#ifndef HAMILTLCAO_H
#define HAMILTLCAO_H

#include "hamilt.h"


namespace ModuleHamilt
{

template<typename T>
class HamiltLCAO : public Hamilt
{
    public:
    //HamiltLCAO();
    //~HamiltLCAO();

    //construct Hamiltonian matrix with inputed electonic density 
    void constructHamilt(const int iter, const MatrixBlock<double> rho) override
    {
        this->ch_mock();
    };
    
    //for target K point, update consequence of hPsi() and matrix()
    void updateHk(int ik) override
    {
        this->hk_mock();
    };
    
    //core function: for solving eigenvalues of Hamiltonian with iterative method
    virtual void hPsi(const ModulePsi::Psi<std::complex<double>>& psi, ModulePsi::Psi<std::complex<double>>& hpsi) const override
    {
        //should be updated for iterative diagonalization method
        return;
    };
    
    //core function: return H(k) and S(k) matrixs for direct solving eigenvalues.
    //not used in PW base
    void matrix(const MatrixBlock<std::complex<double>> hk_in, const MatrixBlock<std::complex<double>> sk_in) override
    {
        return;
    };

    private:
    void ch_mock();
    void hk_mock();

    LocalMatrix<T> localM;

};

//memory for storing Hamiltonian matrix and overlap matrix for one k point
template<typename T>
class LocalMatrix
{
    std::vector<T> Hloc;
    std::vector<T> Sloc;
    size_t size;
    void resize(size_t size_in);
    T* getH();
    T* getS();

};

}


#endif