/* mymath1.cpp file */

#include <iostream>
#include <fstream>
#include <iomanip>

#include <cmath>
#include <cstdlib>
#include "mymath.h"

typedef double real;
#ifdef __FFTW3
void fftw_zeros(fftw_complex *data,int n)
{
    for (int i=0;i<n;i++)
    {
        data[i][0]= 0.0;
        data[i][1]= 0.0;
    }
}
#endif

double rndm()
{
    //generate a random number between [-1, 1]
    double rndNum;
    rndNum = 0.0;
    rndNum = (rand() / (1.0 * RAND_MAX) - 0.5) * 2.0;
    return rndNum;
}

void simpson(const int mesh,const double *func,const double *rab, double &asum)
{
    //     simpson's rule integrator for function stored on the
    //     radial logarithmic mesh

    // integer :: i, mesh

    double f1=0.0;
    double f2=0.0;
    double f3=0.0;
    double r12=1.00 / 12.00;
    int i=0;

    // routine assumes that mesh is an odd number so run check
    // (mesh + 1) ??

    if (mesh%2==0) 	//( mesh+1 - ( (mesh+1) / 2 ) * 2 != 1 )
    {
        cout << "\n error in subroutine simpson ";
        cout << "\n routine assumes mesh is odd but mesh = "
             << mesh << endl;
        // write(*,*) '***error in subroutine radlg';
        // write(*,*) 'routine assumes mesh is odd but mesh =',mesh+1,??;
        exit(0);
        //return ;	//stop;
    }//  endif

    asum = 0.00;

    f3 = func [0] * rab [0] * r12;

    for (i = 1;i < mesh;i+=2)  // do i = 2, mesh - 1, 2
    {
        f1 = f3;
        f2 = func [i] * rab [i] * r12;
        f3 = func [i + 1] * rab [i + 1] * r12;
        asum = asum + 4.00 * f1 + 16.00 * f2 + 4.00 * f3;
    }
    return;
} // end subroutine simpson

void simpson_cp90(int mesh, double *func, double *rab, double intg)
{
    // implicit none;
    // integer mesh;
    // real(kind=8)  func(mesh), rab(mesh), intg;

    // real(kind=8) c(4);
    double c[4];
    int i;

    if (mesh < 8)
        cout << "\n simpson , few mesh points,8 " << endl;

    c[0] = 109.0 / 48.0;

    c[1] = -5.0 / 48.0;

    c[2] = 63.0 / 48.0;

    c[3] = 49.0 / 48.0;

    intg = (func[0] * rab[0] + func[mesh-1] * rab[mesh-1]) * c[0]
           + (func[1] * rab[1] + func[mesh-2] * rab[mesh-2]) * c[1]
           + (func[2] * rab[2] + func[mesh-3] * rab[mesh-3]) * c[2]
           + (func[3] * rab[3] + func[mesh-4] * rab[mesh-4]) * c[3];

    for (i = 4;i < mesh - 4;i++)  //do i=5,mesh-4
    {
        intg = intg + func[i] * rab[i];
    } //  end do

    return;
} //end subroutine simpson_cp90


/*===============================================================
!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
*/

void simpson_fpmd(int n, double *func, double dx, double s)
{
    // INTEGER, INTENT(IN) ::  N
    // REAL(dbl), INTENT(IN)  ::  func(N), dx
    // REAL(dbl), INTENT(OUT) ::  S

    double c0, c1, c2, c3;
    int i;
    //  REAL(dbl) :: C1,C2,C3,C4
    //  PARAMETER(C1=109.d0/48.d0,C2=-5.d0/48.d0, C3=63.d0/48.d0,C4=49.d0/48.d0)
    //  INTEGER I

    c0 = 109.0 / 48.0;
    c1 = -5.0 / 48.0;
    c2 = 63.0 / 48.0;
    c3 = 49.0 / 48.0;

    s =     func[0] * c0;
    s = s + func[1] * c1;
    s = s + func[2] * c2;
    s = s + func[3] * c3;

    for (i = 4;i < n - 5;i++) 	// DO I = 5, (N-5)
    {
        s = s + func[i];
    } //  END DO

    s = s + func[n-4] * c3;

    s = s + func[n-3] * c2;

    s = s + func[n-2] * c1;

    s = s + func[n-1] * c0;

    s = s * dx;

    return;
}//  END SUBROUTINE simpson_fpmd

/*****************************************************
! Copyright (C) 2002-2003 PWSCF+CP group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
********************************************************/

// Error function
double my_erf(double x)
{
    //==================================================================
    // Error function - computed from the rational approximations of
    // W. J. Cody, Math. Comp. 22 (1969), pages 631-637.
    //
    // for abs(x) le 0.47 erf is calculated directly
    // for abs(x) gt 0.47 erf is calculated via erf(x)=1-erfc(x)
    //==================================================================

    // data p1 / 2.42667955230532d2, 2.19792616182942d1, &
    // 6.99638348861914d0, -3.56098437018154d-2 /
    // data q1 / 2.15058875869861d2, 9.11649054045149d1, &
    // 1.50827976304078d1, 1.00000000000000d0 /

    double p1[4] = { 242.667955230532, 21.9792616182942,
                     6.99638348861914, -0.0356098437018154
                   };

    double q1[4] = {215.058875869861, 91.1649054045149,
                    15.0827976304078, 1.00000000000000
                   };
    double erf0;

    if (fabs(x) > 6.0)
    {
        //  erf(6)=1-10^(-17) cannot be distinguished from 1 with 16-byte words
        erf0 = (x < 0) ? -1.0 : 1.0;// sign (1.0, x) ;
    }
    else if (fabs(x) <= 0.470)
    {
        const double x2 = x * x;
        erf0 = x * (p1 [0] + x2 * (p1 [1] + x2 * (p1 [2] + x2 * p1 [3])))
               / (q1 [0] + x2 * (q1 [1] + x2 * (q1 [2] + x2 * q1 [3])));
    }
    else
    {
        erf0 = 1.0 - erfc(x);
    }

    return erf0;
}

// complementary error function
double my_erfc(double x)
{
    double ax, x2, xm2, erfc0;
    // real ax, x2, xm2, p2 [8], q2 [8], p3 [5], q3 [5], pim1;
    double p2[8] = {
        300.459261020162, 451.918953711873, 339.320816734344,
        152.989285046940, 43.1622272220567, 7.21175825088309,
        0.564195517478974, -0.000000136864857382717
    };

    double p3[5] = {
        -0.00299610707703542, -0.0494730910623251, -0.226956593539687,
        -0.278661308609648, -0.0223192459734185
    };

    double q2[8] = {
        300.459260956983, 790.950925327898, 931.354094850610,
        638.980264465631, 277.585444743988, 77.0001529352295,
        12.7827273196294, 1.00000000000000
    };

    double q3[5] = { 0.0106209230528468, 0.191308926107830,
                     1.05167510706793, 1.98733201817135, 1.00000000000000
                   };

    const double pim1 = 0.564189583547756 ;
    ax = fabs(x);

    if (ax > 26.0)
    {
        //  erfc(26.0)=10^(-296); erfc( 9.0)=10^(-37);
        erfc0 = 0.0 ;
    }
    else if (ax > 4.0)
    {
        x2 = x * x;
        xm2 = (1.0 / ax) * (1.0 / ax);
        erfc0 = (1.0 / ax) * exp(- x2) * (pim1 + xm2 * (p3 [0]
                                          + xm2 * (p3 [1] + xm2 * (p3 [2] + xm2 * (p3 [3]
                                                                   + xm2 * p3 [4] )))) / (q3 [0] + xm2 * (q3 [1] + xm2 * (q3 [2] + xm2 *
                                                                                                          (q3 [3] + xm2 * q3 [4])))));
    }
    else if (ax > 0.470)
    {
        x2 = x * x;
        erfc0 = exp(- x2) * (p2 [0] + ax * (p2 [1] +
                                            ax * (p2 [2] + ax * (p2 [3] +
                                                                 ax * (p2 [4] + ax * (p2 [5] +
                                                                                      ax * (p2 [6] + ax * p2 [7]))))))) /
                (q2 [0] + ax * (q2 [1] + ax * (q2 [2] +
                                               ax * (q2 [3] + ax * (q2 [4] +
                                                                    ax * (q2 [5] +ax *(q2 [6] +
                                                                                       ax * q2 [7])))))));
    }
    else
    {
        erfc0 = 1.0 - my_erf(ax) ;
    }

    // erf(-x)=-erf(x)  =>  erfc(-x) = 2-erfc(x)
    if (x < 0.0) erfc0 = 2.0 - erfc0;

    return erfc0 ;
} //end function erfc

double gauss_freq(double x)
{
    const double g_freq = 0.5 * my_erfc(- x * 0.707106781186548);
    return g_freq;
} //end function gauss_freq

