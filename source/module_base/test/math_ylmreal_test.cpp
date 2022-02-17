#include"../math_ylmreal.h"
#include"../vector3.h"
#include"../matrix.h"
#include"gtest/gtest.h"

#define PI 3.141592653589793238462643383279502884197169399
#define doublethreshold 1e-12

/************************************************
*  unit test of class YlmReal and Ylm
***********************************************/

/**
 * For lmax <5 cases, the reference values are calculated by the formula from 
 * https://formulasearchengine.com/wiki/Table_of_spherical_harmonics. Note, these
 * formula lack of the Condon–Shortley phase (-1)^m, and in this unit test, item 
 * (-1)^m is multiplied.
 * For lmax >=5, the reference values are calculated by YlmReal::Ylm_Real.
 *
 * - Tested functions of class YlmReal
 *      - Ylm_Real
 *      - Ylm_Real2
 *      - rlylm 
 */



//mock functions of WARNING_QUIT and WARNING
namespace ModuleBase
{
    void WARNING_QUIT(const std::string &file,const std::string &description) {return ;}
    void WARNING(const std::string &file,const std::string &description) {return ;}
}


class YlmRealTest : public testing::Test
{
    protected:

    int lmax = 7;
    int ng = 4; //test the 4 selected points on the sphere
    int nylm ; // total Ylm number;
    ModuleBase::matrix ylm;
    ModuleBase::Vector3<double> *g;
    double *ref;
    double *rly;

    //Ylm function
    //https://formulasearchengine.com/wiki/Table_of_spherical_harmonics
    //multipy the Condon–Shortley phase (-1)^m
    inline double norm(const double &x, const double &y, const double &z) {return sqrt(x*x + y*y + z*z);}
    double y00(const double &x, const double &y, const double &z)  {return 1.0/2.0/sqrt(PI);}
    double y10(const double &x, const double &y, const double &z)  {double r=norm(x,y,z); return sqrt(3.0/(4.0*PI)) * z / r;}
    double y11(const double &x, const double &y, const double &z)  {double r=norm(x,y,z); return -1.0*sqrt(3.0/(4.*PI)) * x / r;}
    double y1m1(const double &x, const double &y, const double &z) {double r=norm(x,y,z); return -1.0*sqrt(3./(4.*PI)) * y / r;} // y1m1 means Y1,-1
    double y20(const double &x, const double &y, const double &z)  {double r=norm(x,y,z); return 1./4. * sqrt(5./PI) * (-1.*x*x - y*y + 2.*z*z) / (r*r);}
    double y21(const double &x, const double &y, const double &z)  {double r=norm(x,y,z); return -1.0*1./2. * sqrt(15./PI) * (z*x) / (r*r);}
    double y2m1(const double &x, const double &y, const double &z) {double r=norm(x,y,z); return -1.0*1./2. * sqrt(15./PI) * (z*y) / (r*r);}
    double y22(const double &x, const double &y, const double &z)  {double r=norm(x,y,z); return 1./4. * sqrt(15./PI) * (x*x - y*y) / (r*r);}
    double y2m2(const double &x, const double &y, const double &z) {double r=norm(x,y,z); return 1./2. * sqrt(15./PI) * (x*y) / (r*r);}
    double y30(const double &x, const double &y, const double &z)  {double r=norm(x,y,z); return 1./4. * sqrt(7./PI) * z*(2.*z*z-3.*x*x-3.*y*y) / (r*r*r);}
    double y31(const double &x, const double &y, const double &z)  {double r=norm(x,y,z); return -1.0*1./4. * sqrt(21./2./PI) * x*(4.*z*z-x*x-y*y) / (r*r*r);}
    double y3m1(const double &x, const double &y, const double &z) {double r=norm(x,y,z); return -1.0*1./4. * sqrt(21./2./PI) * y*(4.*z*z-x*x-y*y) / (r*r*r);}
    double y32(const double &x, const double &y, const double &z)  {double r=norm(x,y,z); return 1./4. * sqrt(105./PI) * (x*x - y*y)*z / (r*r*r);}
    double y3m2(const double &x, const double &y, const double &z) {double r=norm(x,y,z); return 1./2. * sqrt(105./PI) * x*y*z / (r*r*r);}
    double y33(const double &x, const double &y, const double &z)  {double r=norm(x,y,z); return -1.0*1./4. * sqrt(35./2./PI) * x*(x*x - 3.*y*y) / (r*r*r);}
    double y3m3(const double &x, const double &y, const double &z) {double r=norm(x,y,z); return -1.0*1./4. * sqrt(35./2./PI) * y*(3.*x*x - y*y) / (r*r*r);}
    double y40(const double &x, const double &y, const double &z)  {double r=norm(x,y,z); return 3./16.*sqrt(1./PI) * (35.*z*z*z*z - 30.*z*z*r*r + 3*r*r*r*r) / (r*r*r*r);}
    double y41(const double &x, const double &y, const double &z)  {double r=norm(x,y,z); return -1.0*3./4.*sqrt(5./2./PI) * x*z*(7.*z*z - 3*r*r) / (r*r*r*r);}
    double y4m1(const double &x, const double &y, const double &z) {double r=norm(x,y,z); return -1.0*3./4.*sqrt(5./2./PI) * y*z*(7.*z*z - 3.*r*r) / (r*r*r*r);}
    double y42(const double &x, const double &y, const double &z)  {double r=norm(x,y,z); return 3./8.*sqrt(5./PI) * (x*x-y*y)*(7.*z*z-r*r) / (r*r*r*r);}
    double y4m2(const double &x, const double &y, const double &z) {double r=norm(x,y,z); return 3./4.*sqrt(5./PI) * x*y*(7.*z*z - r*r) / (r*r*r*r);}
    double y43(const double &x, const double &y, const double &z)  {double r=norm(x,y,z); return -1.0*3./4.*sqrt(35./2./PI) * x*z*(x*x - 3.*y*y) / (r*r*r*r);}
    double y4m3(const double &x, const double &y, const double &z) {double r=norm(x,y,z); return -1.0*3./4.*sqrt(35./2./PI) * y*z*(3.*x*x - y*y) / (r*r*r*r);}
    double y44(const double &x, const double &y, const double &z)  {double r=norm(x,y,z); return 3./16.*sqrt(35./PI) * (x*x*(x*x - 3.*y*y) - y*y*(3.*x*x-y*y)) / (r*r*r*r);}
    double y4m4(const double &x, const double &y, const double &z) {double r=norm(x,y,z); return 3./4.*sqrt(35./PI) * x*y*(x*x - y*y) / (r*r*r*r);}

    void SetUp()
    {
        nylm = (lmax + 1) * (lmax + 1);
        ylm.create(nylm,ng);
        g = new ModuleBase::Vector3<double>[ng];
        g[0].set(1.0,0.0,0.0);
        g[1].set(0.0,1.0,0.0);
        g[2].set(0.0,0.0,1.0);
        g[3].set(-1.0,-1.0,-1.0);

        rly = new double[nylm];
        ref = new double[nylm*ng]{
            y00(g[0].x, g[0].y, g[0].z),  y00(g[1].x, g[1].y, g[1].z),  y00(g[2].x, g[2].y, g[2].z),  y00(g[3].x, g[3].y, g[3].z),  
            y10(g[0].x, g[0].y, g[0].z),  y10(g[1].x, g[1].y, g[1].z),  y10(g[2].x, g[2].y, g[2].z),  y10(g[3].x, g[3].y, g[3].z),  
            y11(g[0].x, g[0].y, g[0].z),  y11(g[1].x, g[1].y, g[1].z),  y11(g[2].x, g[2].y, g[2].z),  y11(g[3].x, g[3].y, g[3].z),  
            y1m1(g[0].x, g[0].y, g[0].z), y1m1(g[1].x, g[1].y, g[1].z), y1m1(g[2].x, g[2].y, g[2].z), y1m1(g[3].x, g[3].y, g[3].z), 
            y20(g[0].x, g[0].y, g[0].z),  y20(g[1].x, g[1].y, g[1].z),  y20(g[2].x, g[2].y, g[2].z),  y20(g[3].x, g[3].y, g[3].z),  
            y21(g[0].x, g[0].y, g[0].z),  y21(g[1].x, g[1].y, g[1].z),  y21(g[2].x, g[2].y, g[2].z),  y21(g[3].x, g[3].y, g[3].z),  
            y2m1(g[0].x, g[0].y, g[0].z), y2m1(g[1].x, g[1].y, g[1].z), y2m1(g[2].x, g[2].y, g[2].z), y2m1(g[3].x, g[3].y, g[3].z),                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
            y22(g[0].x, g[0].y, g[0].z),  y22(g[1].x, g[1].y, g[1].z),  y22(g[2].x, g[2].y, g[2].z),  y22(g[3].x, g[3].y, g[3].z),  
            y2m2(g[0].x, g[0].y, g[0].z), y2m2(g[1].x, g[1].y, g[1].z), y2m2(g[2].x, g[2].y, g[2].z), y2m2(g[3].x, g[3].y, g[3].z), 
            y30(g[0].x, g[0].y, g[0].z),  y30(g[1].x, g[1].y, g[1].z),  y30(g[2].x, g[2].y, g[2].z),  y30(g[3].x, g[3].y, g[3].z),  
            y31(g[0].x, g[0].y, g[0].z),  y31(g[1].x, g[1].y, g[1].z),  y31(g[2].x, g[2].y, g[2].z),  y31(g[3].x, g[3].y, g[3].z),  
            y3m1(g[0].x, g[0].y, g[0].z), y3m1(g[1].x, g[1].y, g[1].z), y3m1(g[2].x, g[2].y, g[2].z), y3m1(g[3].x, g[3].y, g[3].z), 
            y32(g[0].x, g[0].y, g[0].z),  y32(g[1].x, g[1].y, g[1].z),  y32(g[2].x, g[2].y, g[2].z),  y32(g[3].x, g[3].y, g[3].z),  
            y3m2(g[0].x, g[0].y, g[0].z), y3m2(g[1].x, g[1].y, g[1].z), y3m2(g[2].x, g[2].y, g[2].z), y3m2(g[3].x, g[3].y, g[3].z), 
            y33(g[0].x, g[0].y, g[0].z),  y33(g[1].x, g[1].y, g[1].z),  y33(g[2].x, g[2].y, g[2].z),  y33(g[3].x, g[3].y, g[3].z),  
            y3m3(g[0].x, g[0].y, g[0].z), y3m3(g[1].x, g[1].y, g[1].z), y3m3(g[2].x, g[2].y, g[2].z), y3m3(g[3].x, g[3].y, g[3].z), 
            y40(g[0].x, g[0].y, g[0].z),  y40(g[1].x, g[1].y, g[1].z),  y40(g[2].x, g[2].y, g[2].z),  y40(g[3].x, g[3].y, g[3].z),  
            y41(g[0].x, g[0].y, g[0].z),  y41(g[1].x, g[1].y, g[1].z),  y41(g[2].x, g[2].y, g[2].z),  y41(g[3].x, g[3].y, g[3].z),  
            y4m1(g[0].x, g[0].y, g[0].z), y4m1(g[1].x, g[1].y, g[1].z), y4m1(g[2].x, g[2].y, g[2].z), y4m1(g[3].x, g[3].y, g[3].z), 
            y42(g[0].x, g[0].y, g[0].z),  y42(g[1].x, g[1].y, g[1].z),  y42(g[2].x, g[2].y, g[2].z),  y42(g[3].x, g[3].y, g[3].z),  
            y4m2(g[0].x, g[0].y, g[0].z), y4m2(g[1].x, g[1].y, g[1].z), y4m2(g[2].x, g[2].y, g[2].z), y4m2(g[3].x, g[3].y, g[3].z), 
            y43(g[0].x, g[0].y, g[0].z),  y43(g[1].x, g[1].y, g[1].z),  y43(g[2].x, g[2].y, g[2].z),  y43(g[3].x, g[3].y, g[3].z),  
            y4m3(g[0].x, g[0].y, g[0].z), y4m3(g[1].x, g[1].y, g[1].z), y4m3(g[2].x, g[2].y, g[2].z), y4m3(g[3].x, g[3].y, g[3].z), 
            y44(g[0].x, g[0].y, g[0].z),  y44(g[1].x, g[1].y, g[1].z),  y44(g[2].x, g[2].y, g[2].z),  y44(g[3].x, g[3].y, g[3].z),  
            y4m4(g[0].x, g[0].y, g[0].z), y4m4(g[1].x, g[1].y, g[1].z), y4m4(g[2].x, g[2].y, g[2].z), y4m4(g[3].x, g[3].y, g[3].z), 
              0.000000000000000,    0.000000000000000,    0.935602579627389,    0.090028400200397, 
             -0.452946651195697,   -0.000000000000000,   -0.000000000000000,   -0.348678494661834, 
             -0.000000000000000,   -0.452946651195697,   -0.000000000000000,   -0.348678494661834, 
             -0.000000000000000,    0.000000000000000,    0.000000000000000,   -0.000000000000000, 
             -0.000000000000000,   -0.000000000000000,    0.000000000000000,   -0.000000000000000, 
              0.489238299435250,    0.000000000000000,   -0.000000000000000,   -0.376615818502422, 
              0.000000000000000,   -0.489238299435250,   -0.000000000000000,    0.376615818502422, 
              0.000000000000000,    0.000000000000000,    0.000000000000000,    0.532615198330370, 
              0.000000000000000,    0.000000000000000,    0.000000000000000,   -0.000000000000000, 
             -0.656382056840170,   -0.000000000000000,   -0.000000000000000,   -0.168427714314628, 
             -0.000000000000000,   -0.656382056840170,   -0.000000000000000,   -0.168427714314628, 
             -0.317846011338142,   -0.317846011338142,    1.017107236282055,    0.226023830284901, 
             -0.000000000000000,   -0.000000000000000,   -0.000000000000000,    0.258942827786103, 
             -0.000000000000000,   -0.000000000000000,   -0.000000000000000,    0.258942827786103, 
              0.460602629757462,   -0.460602629757462,    0.000000000000000,   -0.000000000000000, 
              0.000000000000000,    0.000000000000000,    0.000000000000000,   -0.409424559784410, 
             -0.000000000000000,   -0.000000000000000,   -0.000000000000000,    0.136474853261470, 
             -0.000000000000000,    0.000000000000000,   -0.000000000000000,   -0.136474853261470, 
             -0.504564900728724,   -0.504564900728724,    0.000000000000000,   -0.598002845308118, 
             -0.000000000000000,   -0.000000000000000,    0.000000000000000,    0.000000000000000, 
             -0.000000000000000,   -0.000000000000000,   -0.000000000000000,    0.350610246256556, 
             -0.000000000000000,   -0.000000000000000,   -0.000000000000000,    0.350610246256556, 
              0.683184105191914,   -0.683184105191914,    0.000000000000000,   -0.000000000000000, 
              0.000000000000000,    0.000000000000000,    0.000000000000000,   -0.202424920056864, 
              0.000000000000000,    0.000000000000000,    1.092548430592079,   -0.350435072502801, 
              0.451658037912587,    0.000000000000000,   -0.000000000000000,    0.046358202625865, 
              0.000000000000000,    0.451658037912587,   -0.000000000000000,    0.046358202625865, 
              0.000000000000000,   -0.000000000000000,    0.000000000000000,    0.000000000000000, 
              0.000000000000000,    0.000000000000000,    0.000000000000000,    0.492067081245654, 
             -0.469376801586882,   -0.000000000000000,   -0.000000000000000,    0.187354445356332, 
             -0.000000000000000,    0.469376801586882,   -0.000000000000000,   -0.187354445356332, 
              0.000000000000000,    0.000000000000000,    0.000000000000000,    0.355076798886913, 
              0.000000000000000,    0.000000000000000,    0.000000000000000,   -0.000000000000000, 
              0.518915578720260,    0.000000000000000,   -0.000000000000000,   -0.443845998608641, 
              0.000000000000000,    0.518915578720260,   -0.000000000000000,   -0.443845998608641, 
              0.000000000000000,   -0.000000000000000,    0.000000000000000,    0.000000000000000, 
              0.000000000000000,    0.000000000000000,    0.000000000000000,    0.452635881587108, 
             -0.707162732524596,    0.000000000000000,   -0.000000000000000,    0.120972027847095, 
             -0.000000000000000,    0.707162732524596,   -0.000000000000000,   -0.120972027847095  
         } ; 
        
           
    }

    void TearDown()
    {
        delete [] g;
        delete [] ref;
        delete [] rly;
    }
};

TEST_F(YlmRealTest,YlmReal)
{
    ModuleBase::YlmReal::Ylm_Real(nylm,ng,g,ylm);
    for(int i=0;i<nylm;++i)
    {
        for(int j=0;j<ng;++j) 
        {
            EXPECT_NEAR(ylm(i,j),ref[i*ng+j],doublethreshold)  << "i=" << i << " ,j=" << j;
        }
    } 
}


TEST_F(YlmRealTest,YlmReal2)
{
    ModuleBase::YlmReal::Ylm_Real2(nylm,ng,g,ylm);
    for(int i=0;i<nylm;++i)
    {
        for(int j=0;j<ng;++j) 
        {
            EXPECT_NEAR(ylm(i,j),ref[i*ng+j],doublethreshold) << "i=" << i << " ,j=" << j;
        }
    } 
}


TEST_F(YlmRealTest,rlylm)
{    
    for(int j=0;j<ng;++j)
    {
        ModuleBase::YlmReal::rlylm(lmax,g[j].x,g[j].y,g[j].z,rly);
        for(int i=0;i<nylm;++i)
        {
            EXPECT_NEAR(rly[i],ref[i*ng+j],doublethreshold) << "i=" << i << " ,j=" << j;
        }
    }
}

