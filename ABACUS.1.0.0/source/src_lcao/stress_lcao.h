#ifndef STRESS_LCAO_H
#define STRESS_LCAO_H

#include "../src_pw/tools.h"

#include "../src_pw/stress.h" //test

class Stress_LCAO
{
	public :
	
	Stress_LCAO ();
	~Stress_LCAO ();

	void allocate (void);
	void destroy (void);

	void cal_stress_loc (void);
	void cal_stress_ew (void);
	void cal_stress_scc (void);
	void cal_stress_cc (void);
        void cal_stress_har(void);
	void cal_stress_gradcorr(void);
	
	void print_stress(const string &name, double f[][3], const bool screen, bool ry)const;
	void printstress_total (bool ry);
	
	void start_stress(double overlap[][3],double tvnl_dphi[][3],double vnl_dbeta[][3],double vl_dphi[][3]);
	
	// total stress
	double scs[3][3];

        static double stress_invalid_threshold_ev;
	private:
	
	//each part of stress
	double sigmadvl[3][3];
	double sigmaewa[3][3]; 
	double sigmacc[3][3];
	double sigmaxc[3][3];
        double sigmahar[3][3];
        //each part of stress
        //only calculated in force_lcao part
        double soverlap[3][3];
        double stvnl_dphi[3][3];
        double svnl_dbeta[3][3];
        double svl_dphi[3][3];

        bool allocate_flag;
	double output_acc; // control the accuracy
   
        Stress str;    //test
};
#endif
