//==========================================================
// Author: Xin Qu
// DATE : 2019-12-10
//==========================================================
#ifndef DFTU_H
#define DFTU_H

#include <string>

#include "../src_global/complexmatrix.h"
#include "../src_pw/charge_broyden.h"

using namespace std;

//==========================================================
// CLASS :
// NAME : DTFU (DFT+U)
//==========================================================
class DFTU
{
public:
    DFTU();                      // constructor 
    ~DFTU();                     // deconstructor

    //called at Run_Frag::frag_LCAO_line(void)
    void init();                        // initialize the input terms of  U, J, double_counting etc
    
    //called at Local_Orbital_Elec::scf(const int &istep)
    //calculate the local occupation number matrix
    void cal_occup_m_k(const int iter);
    void cal_occup_m_gamma(const int iter);

    void write_occup_m(const string &fn);
    void read_occup_m(const string &fn);
    void local_occup_bcast();
    
    //calculate the energy correction: en_cor
    //called at energy::calculate_etot(void)
    void cal_energy_correction( const int istep);

    //calculate the effective potential
    void cal_eff_pot_mat(const int ik, const int istep);
    double get_onebody_eff_pot
    (
        const int T, const int iat,
	    const int L, const int N, const int spin, 
	    const int m0,
        const int m1,
        const int type, const bool newlocale
    );


    //calculate force and stress
    void cal_force_stress_gamma();
    void force_stress();
    void folding_dSm_soverlap();
    void allocate_force_stress();
    void erase_force_stress();
    void cal_force_k(vector<vector<complex<double>>> &VU);
    void cal_force_gamma(vector<vector<double>> &VU);
    void cal_stress_k(vector<vector<complex<double>>> &VU);
    void cal_stress_gamma(vector<vector<double>> &VU);
    
    //LSCC
    void cal_slater_Fk(const int L, const int T); //L:angular momnet, T:atom type
    void cal_unscreened_slater_Fk(const int L, const int T); //L:angular momnet, T:atom type
    void cal_yukawa_lambda();
    void cal_slater_UJ(const int istep, const int iter);
    void cal_slater_Vsc(const int T, const int L);


   // void print(const int T, const int iat, const int L, const int N, const int iter);

    void output();
    
    //Sm_k[ik][irc]: for k_points algorithm, calculated at LCAO_nnr::folding_fixedH(const int &ik)
    vector<vector<complex<double>>> Sm_k;


    // effective potential matrix for k algorithm: pot_eff_k[ik][irc]
    // effective potential matrix for gamma only algorithm: pot_eff_gamma[is][irc]
    vector<vector<complex<double>>> pot_eff_k;
    vector<vector<double>> pot_eff_gamma;



    double EU;
    int iter_dftu;

    //force and stress
    vector<vector<double>> force_dftu;      //force_dftu[iat][dim] 
    vector<vector<double>> stress_dftu;


private:

    int cal_type;        //1:dftu_tpye=1, dc=1; 2:dftu_type=1, dc=2; 3:dftu_tpye=2, dc=1; 4:dftu_tpye=2, dc=2; 
    int dc;              //dc (type of double_counting)
    bool Yukawa;         //1:use Yukawa potential; 0: do not use Yukawa potential 
    double lambda;       //the parameter in Yukawa potential

    double *U;           //J (Hund parameter J)
    double *J;           //J (Hund parameter J)
    double Nval;         //Total nmuber of valence electrons of the system 
    double Nc;           //Total nmuber of correlated electrons of the system 
    
    vector<vector<matrix>> Vsc; //Vsc[T][N](i,j)
    vector<vector<vector<vector<double>>>> Fk; //slater integral:Fk[T][L][N][k]
    vector<vector<vector<double>>> U_Yukawa;   //U_Yukawa[T][L][N]
    vector<vector<vector<double>>> J_Yukawa;   //J_Yukawa[T]{L][N]

    //transform between iwt index and it, ia, L, N and m index
    vector<vector<vector<vector<vector<int>>>>> iatlnmipol2iwt;   //iatlnm2iwt[iat][l][n][m][ipol]
    vector<int> iwt2it;                               //iwt2it[iwt]
    vector<int> iwt2iat;                              //iwt2iat[iwt]
    vector<int> iwt2l;                                //iwt2l[iwt]
    vector<int> iwt2n;                                //iwt2n[iwt]
    vector<int> iwt2m;                                //iwt2m[iwt]
    vector<int> iwt2ipol;                             //iwt2ipol[iwt]
    vector<int> iat2it;

    //local occupancy matrix of the correlated subspace
    //locale: the out put local occupation number matrix of correlated electrons in the current electronic step
    //locale_save: the input local occupation number matrix of correlated electrons in the current electronic step
    //locale_save2: the input local occupation number matrix of correlated electrons in the last electronic step
    vector<vector<vector<vector<matrix>>>> locale;            // locale[iat][l][n][spin](m1,m2)
    vector<vector<vector<vector<matrix>>>> locale_save;       // locale_save[iat][l][n][spin](m1,m2)

    complex<double> ***dSm_k;                                   //dSm_k[ik][dim][irc]
    complex<double> ***soverlap_k;                              //soverlap_k[ik][xy][irc]
    double **soverlap_gamma;                                    //soverlap_gamma[xy][irc]
   
};

extern DFTU dftu;

#endif