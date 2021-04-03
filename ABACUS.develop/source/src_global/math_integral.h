#ifndef MATH_INTEGRAL_H
#define MATH_INTEGRAL_H

// mohan add 2021-04-03

class Integral
{

	public:

    Integral();
    ~Integral();

	// Peize Lin accelerate 2017-10-02
    static void Simpson_Integral
    (
        const int mesh,
        const double *func,
        const double *rab,
        double &asum
    );

	// Peize Lin accelerate 2017-10-02
	static void Simpson_Integral
	(
		const int mesh,
		const double *func,
		const double dr,
		double &asum
	);

    // Peize Lin add 2016-02-14
    static void Simpson_Integral_0toall
    (
        const int mesh,
        const double *func,
        const double *rab,
        double *asum
    );

    // Peize Lin add 2016-02-14
    static void Simpson_Integral_alltoinf
    (
        const int mesh,
        const double *func,
        const double *rab,
        double *asum
    );     

};
#endif
