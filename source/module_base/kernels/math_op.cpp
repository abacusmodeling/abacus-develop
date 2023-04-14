#include "module_base/kernels/math_op.h"
#include "module_base/libm/libm.h"

namespace ModuleBase {

template <typename FPTYPE>
__inline__
FPTYPE __fact(const int n) {
    FPTYPE f = 1.0;
    for (int i = n; i > 1; i--) {
        f *= i;
    }
    return f;
}

__inline__
int __semi_fact(const int n)
{
    int semif = 1;
    for (int i = n; i > 2; i -= 2)
    {
        semif *= i;
    }
    return semif;
}

template <typename FPTYPE>
struct cal_ylm_real_op<FPTYPE, psi::DEVICE_CPU> {
    void operator()(
            const psi::DEVICE_CPU *ctx,
            const int &ng,
            const int &lmax,
            const FPTYPE &SQRT2,
            const FPTYPE &PI,
            const FPTYPE &PI_HALF,
            const FPTYPE &FOUR_PI,
            const FPTYPE &SQRT_INVERSE_FOUR_PI,
            const FPTYPE *g,
            FPTYPE * p,
            FPTYPE * ylm)
    {
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int ig = 0; ig < ng; ig++) {
            //----------------------------------------------------------
            // EXPLAIN : if lmax = 1,only use Y00 , output result.
            //----------------------------------------------------------
            if (lmax == 0) {
                ylm[0 * ng + ig] = SQRT_INVERSE_FOUR_PI;
                continue;
            }
            //----------------------------------------------------------
            // LOCAL VARIABLES :
            // NAME : cost = cos(theta),theta and phi are polar angles
            // NAME : phi
            //----------------------------------------------------------
            const FPTYPE gmod = sqrt(g[ig * 3 + 0] * g[ig * 3 + 0] + g[ig * 3 + 1] * g[ig * 3 + 1] + g[ig * 3 + 2] * g[ig * 3 + 2]);
            FPTYPE cost = gmod < 1.0e-9 ? 0.0 : g[ig * 3 + 2] / gmod;
            FPTYPE phi;
            //  beware the arc tan, it is defined modulo pi
            if (g[ig * 3 + 0] > 1.0e-9) {
                phi = atan(g[ig * 3 + 1] / g[ig * 3 + 0]);
            }
            else if (g[ig * 3 + 0] < -1.e-9) {
                phi = atan(g[ig * 3 + 1] / g[ig * 3 + 0]) + PI;
            }
            else {
                phi = PI_HALF * ((g[ig * 3 + 1] >= 0.0) ? 1.0 : -1.0); //HLX: modified on 10/13/2006
            } // end if
            //==========================================================
            // NAME : p(Legendre Polynomials) (0 <= m <= l)
            //==========================================================
            int lm = -1;
            for (int l = 0; l <= lmax; l++) {
                const FPTYPE c = sqrt((2 * l + 1) / FOUR_PI);
                if (l == 0) {
                    p[0 * (lmax + 1) * ng + 0 * ng + ig] = 1.0;
                }
                else if (l == 1) {
                    p[0 * (lmax + 1) * ng + 1 * ng + ig] = cost;
                    p[1 * (lmax + 1) * ng + 1 * ng + ig] = -sqrt(std::max(0.0, 1.0 - cost * cost));
                }
                else {
                    const int l1 = l - 1,
                            l2 = l - 2,
                            l3 = 2 * l - 1;
                    //  recursion on l for P(:,l,m)
                    for (int m = 0; m <= l2; m++) {  // do m = 0, l - 2//mohan modify 2007-10-13
                        p[m * (lmax + 1) * ng + l * ng + ig] =
                                (cost * l3 * p[m * (lmax + 1) * ng + l1 * ng + ig] -
                                 (l1 + m) * p[m * (lmax + 1) * ng + l2 * ng + ig]) / (l - m);
                    } // end do
                    p[l1 * (lmax + 1) * ng + l * ng + ig] =
                            cost * l3 * p[l1 * (lmax + 1) * ng + l1 * ng + ig];
                    FPTYPE x2 = std::max(0.0, 1.0 - cost * cost);
                    p[l * (lmax + 1) * ng + l * ng + ig] = __semi_fact(l3) * pow(x2, static_cast<FPTYPE>(l) / 2.0);//mohan modify 2007-10-13
                    if (l % 2 == 1) {
                        p[l * (lmax + 1) * ng + l * ng + ig] *= -1;
                    }
                } // end if

                // Y_lm, m = 0
                ++lm;
                ylm[lm * ng + ig] = c * p[0 * (lmax + 1) * ng + l * ng + ig];

                for (int m = 1; m <= l; m++) {
                    // Y_lm, m > 0
                    const FPTYPE same =
                            c * sqrt(__fact<FPTYPE>(l - m) /
                                     __fact<FPTYPE>(l + m)) * SQRT2;
                    FPTYPE sinp, cosp;
                    ModuleBase::libm::sincos(m * phi, &sinp, &cosp);
                    ++lm;
                    ylm[lm * ng + ig] = same * p[m * (lmax + 1) * ng + l * ng + ig] * cosp;

                    // Y_lm, m < 0
                    ++lm;
                    ylm[lm * ng + ig] = same * p[m * (lmax + 1) * ng + l * ng + ig] * sinp;
                }
           }// end do
        }
    }
};

template struct cal_ylm_real_op<float, psi::DEVICE_CPU>;
template struct cal_ylm_real_op<double, psi::DEVICE_CPU>;

}  // namespace ModuleBase

