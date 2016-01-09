#ifndef EXX_DIV_H
#define EXX_DIV_H

class Exx_Divergence
{
	public:

	Exx_Divergence();
	~Exx_Divergence();

	void init(void);
	void get_factor(const int &ik, const int &iq);
	

	double* factor;

	private:

	bool gamma_extra;
	double grid_fac;


};

#endif
