#ifndef IONS_MOVE_CG_H
#define IONS_MOVE_CG_H

#include "../src_pw/tools.h"
class Ions_Move_CG
{
public:
    Ions_Move_CG();
    ~Ions_Move_CG();

    void allocate(void);
    void start(const matrix& force, const double &etot);
        
        static double CG_THRESHOLD;
	int sd_step;
	int cg_step;

private:
	double* pos0;
	double* grad0;
        double* cg_grad0;
	double* move0;
	double  e0;
        // setup gradients.
	void setup_cg_grad(double *grad, const double *grad0, double *cg_grad, const double *cg_grad0, const int &ncggrad);
	void setup_move( double *move, double *cg_gradn, const double &trust_radius );
	void setup_etot_cg(const double &energy_in, const bool sd, const bool trial);
        void Brent(double &fa,double &fb,double &fc,double &xa,double &xb,double &xc,double &best_x,double &xpt);
        void f_cal(const double *g0,const double *g1,const int &dim,double &f_value);
        void third_order(const double &e0,const double &e1,const double &fa, const double &fb, const double x, double &best_x);
};


#endif
