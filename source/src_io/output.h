//==========================================================
// AUTHOR : Lixin He, Mohan Chen
// DATA : 2006-11 ~ 2008-11
//==========================================================
#ifndef OUTPUT_H
#define OUTPUT_H

//#include "../src_pw/tools.h"
#include "../module_base/realarray.h"
#include "../module_base/matrix3.h"
#include "../module_base/complexmatrix.h"
class output
{
public:
    //============================
    // Print realArray (3D or 4D)
    //============================
    void printr3_d(std::ofstream &ofs,const std::string &s,const ModuleBase::realArray &u) const;
    void printr4_d(std::ofstream &ofs,const std::string &s,const ModuleBase::realArray &u) const;

    //===========================
    // print matrix3
    //===========================
    void printM3(std::ofstream &ofs,const std::string& description, const ModuleBase::Matrix3 &m)const;
    void printM3(const std::string &description, const ModuleBase::Matrix3 &m)const;

    //===============================
    // print matrix
    //===============================
    void printrm(std::ofstream &ofs,const std::string &s, const ModuleBase::matrix &m, const double &limit = 1.0e-15) const;
    void printrm(const std::string &s, const ModuleBase::matrix &m, const double &limit = 1.0e-15) const;


    //===============================
    // print ModuleBase::ComplexMatrix
    //===============================
    void printcm(std::ofstream &ofs,const std::string &s, const ModuleBase::ComplexMatrix &m) const;

    void printcm(const std::string &s, const ModuleBase::ComplexMatrix &m) const;

    void printcm_real(const std::string &s, const ModuleBase::ComplexMatrix &m,const double &limit = 1.0e-15) const;

    void printcm_real_limit_hermit(const std::string &s, const ModuleBase::ComplexMatrix &m,const double &limit) const;

    void printcm_imag(const std::string &s, const ModuleBase::ComplexMatrix &m,const double &limit = 1.0e-15) const;
    void printcm_norm(const std::string &s, const ModuleBase::ComplexMatrix &m, const double &limit)const;
    void printcm_norm(std::ofstream &ofs, const std::string &s, const ModuleBase::ComplexMatrix &m, const double &limit)const;


    //***************
    // Template
    //***************
public:

    template <class T>
    void printr1_d(std::ofstream &ofs, const std::string &s,T *u, int n1) const
    {
        ofs<<"\n\n "<<s<< "  n1 = "<< n1;
        if (n1<=0)return;
        if (u==0) return;
        //	ofs.setf(ios::scientific,ios::floatfield);
        for (int i=0;i<n1;i++)
        {
            if (i%8==0)ofs<<"\n";
            ofs<< std::setw(12)<<u[i];
        }
        //	ofs.unsetf(ios::scientific);
    }

    //===================================================
    // print one dimension array (double,int,std::string ...)
    //===================================================
    template <class T>
    void printr1_d(const std::string &s, T *u,const int n1) const
    {
        std::cout << "\n " << s << "  Dimension = " << n1;
        if (n1 <= 0)return;
        if (u == 0) return;
        //std::cout.setf(ios::scientific, ios::floatfield);
        for (int i = 0;i < n1;i++)
        {
            if (i % 8 == 0) std::cout << "\n";
            std::cout << std::setw(12) << u[i];
        }
        std::cout<<std::endl;
        //std::cout.unsetf(ios::scientific);
    }

    //================================
    // output ModuleBase::Vector3 : printV3
    // output ModuleBase::Vector3[] : printv31_d
    //================================
    template <class T>
    void printV3(std::ofstream &ofs, const ModuleBase::Vector3 <T> v)const
    {
        ofs << " ";
        ofs << std::setw(18) << v.x << std::setw(18) << v.y << std::setw(18) << v.z << std::endl;
    }

    template <class T>
    void printV3(const ModuleBase::Vector3 <T> v)const
    {
        std::cout << " ";
        std::cout << std::setw(18) << v.x << std::setw(18) << v.y << std::setw(18) << v.z << std::endl;
    }

    template <class T>
    void printv31_d(std::ofstream &ofs, const std::string &s, ModuleBase::Vector3<T> *u, int n1) const
    {
        ofs << " " << s << " Dimension = " << n1 << std::endl;
        if (n1 <= 0)return;
        if (u == 0) return;
        for (int i = 0;i < n1;i++)
        {
            printV3(ofs, u[i]);
        }
    }

    template <class T>
    void printv31_d(const std::string &s, ModuleBase::Vector3<T> *u, int n1) const
    {
        std::cout << "\n " << s << "  dimension = " << n1;
        if (n1 <= 0)return;
        if (u == 0) return;
        for (int i = 0;i < n1;i++)
        {
            printV3(u[i]);
        }
    }

};

#endif // OUTPUT_H
