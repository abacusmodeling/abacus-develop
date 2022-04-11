#ifndef HAMILTLCAO_H
#define HAMILTLCAO_H

#include "hamilt.h"


namespace ModuleHamilt
{

//memory for storing Hamiltonian matrix and overlap matrix for one k point
template<typename T>
class LocalMatrix
{
    public:
    std::vector<T> hloc;
    std::vector<T> sloc;
    size_t size=1;
    void resize(size_t size_in){size = size_in; hloc.resize(size); sloc.resize(size);}
    T* getH(){return hloc.data();}
    T* getS(){return sloc.data();}

};

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
    void matrix(const MatrixBlock<std::complex<double>> hk_in, const MatrixBlock<std::complex<double>> sk_in) override;

    private:
    void ch_mock();
    void hk_mock();

    LocalMatrix<T> localM;

};



}


#endif