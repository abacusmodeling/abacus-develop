#ifndef CHARGE_EXTRA_H
#define CHARGE_EXTRA_H

using namespace std;

class Charge_Extra
{
	public:
	Charge_Extra();
	~Charge_Extra();

	void allocate(void);
	void record_rho(void);
	void extrapolate_charge(void);
	//xiaohui add 2014-05-07, to use "istep = ions.istep"
	int istep;
	//xiaohui add 2014-05-10, for second-order extrapolation
	double* pos_old1;
	double* pos_old2;
	double* pos_now;
	double* pos_next;

	private:

	double*** rho_ion; //(dim, nspin, pw.nrxx)
	bool init_rho;
	int dim;
	//xiaohui add 2014-05-07, to use first-order extrapolation
	double** delta_rho1;
	double** delta_rho2;
	double** delta_rho;
	//xiaohui add 2014-05-10, to use second-order extrapolation
	double** delta_rho3;
	int pos_dim;
	double alpha,beta;


	void find_alpha_and_beta(void);


};

#endif
