#ifndef SPH_BESSEL_UNITTEST_H
#define SPH_BESSEL_UNITTEST_H

#include<sys/time.h>
#include "../../../module_base/global_function.h"
#include "../../../module_base/sph_bessel.h"

static void Sph_Bessel_timetest(
	const int Lmax,
	const double dR,
	const size_t Rmesh,
	const double dk,
	const size_t kmesh
)
{
	std::cout<<Lmax<<"\t"<<dR<<"\t"<<Rmesh<<"\t"<<dk<<"\t"<<kmesh<<std::endl;
	
	auto generate_point = []( const double d, const size_t mesh ) -> std::vector<double>
	{
		std::vector<double> point(mesh);
		for( size_t i=0; i!=mesh; ++i )
			point[i] = i*d;
		return point;
	};

	const std::vector<double> rpoint = generate_point( dR, Rmesh );
	const std::vector<double> kpoint = generate_point( dk, kmesh );
	
	ModuleBase::Sph_Bessel SB;
	for (int l = 0; l < Lmax+1; l++)
	{
		std::vector<std::vector<double>> jlx( Rmesh, std::vector<double>(kmesh) );
		timeval t_start;	gettimeofday( &t_start, NULL);
		for (int ir = 0; ir < Rmesh; ir++)
		{
			SB.jlx( kmesh, ModuleBase::GlobalFunc::VECTOR_TO_PTR(kpoint), rpoint[ir], l, ModuleBase::GlobalFunc::VECTOR_TO_PTR(jlx[ir]) );
		}
		timeval t_end;	gettimeofday( &t_end, NULL);
		std::cout<<l<<"\t"<<(double)(t_end.tv_sec-t_start.tv_sec) + (double)(t_end.tv_usec-t_start.tv_usec)/1000000.0<<std::endl;
	}
}

#endif // SPH_BESSEL_UNITTEST_H