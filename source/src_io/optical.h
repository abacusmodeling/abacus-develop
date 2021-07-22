#ifndef OPTICAL_H
#define OPTICAL_H

class Optical
{
	public:

	Optical();
	~Optical();
	
	static bool opt_epsilon2;
	static int  opt_nbands;
	
	void cal_epsilon2(const int &nbands);

	private:
	double element_cvk(const int &ik, const int &iv, const int &ic);


};

#endif
