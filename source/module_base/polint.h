#ifndef POLINT_H
#define POLINT_H

namespace ModuleBase
{

class Polint
{
	public:

		Polint();
		~Polint();

		static double Lagrange3
		(
		 	const double* xa,
			const double* ya,
			const int& n,
			const double& x
		);
			
		static double RadialF
		(
		 	const double* rad,
			const double* rad_f,
			const int& msh,
			const int& l,
			const double& R
		);

};

}

#endif
