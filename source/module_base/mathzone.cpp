//#include "../src_pw/global.h"
#include "mathzone.h"
#include "constants.h"
#include "global_variable.h"
#include "global_function.h"
#include "timer.h"

namespace ModuleBase
{

Mathzone::Mathzone()
{
}

Mathzone::~Mathzone()
{
}


void Mathzone::To_Polar_Coordinate
(
    const double &x_cartesian,
    const double &y_cartesian,
    const double &z_cartesian,
    double &r,
    double &theta,
    double &phi
)
{
//	OUT("x",x_cartesian);
//	OUT("y",y_cartesian);
//	OUT("z",z_cartesian);
    r = sqrt(x_cartesian*x_cartesian
             + y_cartesian*y_cartesian
             + z_cartesian*z_cartesian);
//----------------------------------------------------------
// Calculate theta
//----------------------------------------------------------
    double small = 1.0e-9;
    if (r < small)
    {
        theta = ModuleBase::PI/2.0;
    }
    else
    {
        theta = acos(z_cartesian / r);
    }
//----------------------------------------------------------
// Calculate phi
//----------------------------------------------------------
    if (x_cartesian > small && y_cartesian > small)
    {
        phi = atan( y_cartesian / x_cartesian );
    }
    else if (  x_cartesian < small  )
    {
        phi = atan( y_cartesian / x_cartesian ) + ModuleBase::PI;
    }
    else if ( x_cartesian > small && y_cartesian < -small)
    {
        phi = atan( y_cartesian/x_cartesian ) + ModuleBase::TWO_PI;
    }
    else
    {
        phi = ModuleBase::PI/2.0 * ( (y_cartesian >= 0.0) ? 1.0 : -1.0);
    }
//----------------------------------------------------------
// degress = radians / ModuleBase::PI * 180
//----------------------------------------------------------
    theta = theta/ModuleBase::PI*180;
    phi = phi/ModuleBase::PI*180;

//	OUT("r",r);
//	OUT("theta(degree)",theta);
//	OUT("phi(degree)",phi);
    return;
}

}
