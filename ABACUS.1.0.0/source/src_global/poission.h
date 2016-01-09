#ifndef POISSION_H
#define POISSION_H

class Poission
{
	public:
	
	Poission(){};
	~Poission(){};

	static void SolPoissonEq
	(
    	const double* rho,
    	const double* r,
    	const int& mesh,
    	double* pot
	);

};

#endif
