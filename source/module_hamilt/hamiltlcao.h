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

template<typename T, typename T1>
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
    void updateHk(const int ik) override
    {
        this->hk_fixed_mock(ik);
        this->hk_update_mock(ik);
    };
    
    //core function: for solving eigenvalues of Hamiltonian with iterative method
    virtual void hPsi(const ModulePsi::Psi<std::complex<double>>& psi, ModulePsi::Psi<std::complex<double>>& hpsi) const override
    {
        //should be updated for iterative diagonalization method
        return;
    };
    
    //core function: return H(k) and S(k) matrixs for direct solving eigenvalues.
    //not used in PW base
    void matrix(MatrixBlock<std::complex<double>> hk_in, MatrixBlock<std::complex<double>> sk_in)const override;
    void matrix(MatrixBlock<double> hk_in, MatrixBlock<double> sk_in)const override;

    private:
    void ch_mock();
    void hk_fixed_mock(const int ik);
    void hk_update_mock(const int ik);

    //specific code for matrix()
    void getMatrix(MatrixBlock<T> hk_in, MatrixBlock<T> sk_in)const;

    //there are H and S matrix for each k point in reciprocal space
    //type double for gamma_only case, type complex<double> for multi-k-points case
    LocalMatrix<T> kM;
    //there are H and S matrix for fixed T+VNL and overlap terms in real space
    //type double for nspin<4, type complex<double> for nspin==4
    LocalMatrix<T1> fixedRealM;

};



}


#endif