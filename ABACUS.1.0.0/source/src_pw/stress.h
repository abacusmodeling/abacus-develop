#ifndef STRESS_H
#define STRESS_H

#include "tools.h"

using namespace std;
class Stress
{
public: 
         Stress(){};
         ~Stress(){};
         void cal_stress();



       void dvloc_of_g (const int& msh,
                         const double* rab,
                         const double* r,
                         const double* vloc_at,
                         const double& zp,
                         double*  dvloc);
       void dvloc_coul (const double& zp, double* dvloc);

       void deriv_drhoc (
                         const bool &numeric,
                         const int mesh,
                         const double *r,
                         const double *rab,
                         const double *rhoc,
                         double *drhocg);

        void print_stress(const string &name, double f[][3], const bool screen, bool ry)const;
        void printstress_total (bool ry);
	double sigmatot[3][3];

private:
       void stres_knl();
       void stres_har();
       void stres_ewa();
       void stres_loc();
       void stres_cc();
       void stres_gradcorr();

       void stres_nl();
       void dylmr2 (
                         const int nylm,
                         const int ngy,
                         Vector3<double> *gk,
                         matrix &dylm,
                         const int ipol);
       double Polynomial_Interpolation_nl(
                         const realArray &table,
                         const int &dim1,
                         const int &dim2,
                         const double &table_interval,
                         const double &x);
       void get_dvnl2(
                         ComplexMatrix &vkb,
                         const int ik);
       void get_dvnl1(
                         ComplexMatrix &vkb,
                         const int ik,
                         const int ipol);



       double sigmaxc[3][3];
       double sigmakin[3][3];
       double sigmahar[3][3];
       double sigmaloc[3][3];
       double sigmanlc[3][3];
       double sigmaewa[3][3];
       double sigmaxcc[3][3];
}; 







#endif
