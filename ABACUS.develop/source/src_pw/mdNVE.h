#ifndef MDNVE_H
#define MDNVE_H
#include"md.h"

using namespace std;
class mdNVE:public md{
public:
	mdNVE(int n=1):md(1){};
	void runNVE(int step1);
	void RescaleVelocity(double KE,double newKE);
	double Conserved(double KE, double PE);
	double MAXVALF();

private:
	                      //>> INTERNAL VARIABLES <<//
//  int nfrozen;               // Number of frozen atoms
   int ii,k,it,ion;
   int i;                       // Counters to parse through the ion table.
   int step;                    // step number we are on
   //int a.nat;                  // total number of atoms
   int nStep;                   // Total number of steps
   int errorFlag;

  //vector3<double> velocity;           // Velocity of ions (atomic units)
  Vector3<double> fracStep;           // Step (due to velocity) in fractional coordinates
  double mass;   // atom mass, temperay var
   double maxStep;     // Largest movement by an ion
   double hamiltonian,conservedE,maxForce;
   double msd;                  //  mean square displacement
   double diffuCoeff;           // diffusion coefficient
   double newKE;
   double twiceKE;     // Kinetic energy x 2
   double oldEtot;        // old energy

  //TYPE(stopwatch):: &           // Timer objects
  //  watch, &
  //  watch2
  
  //string parameterInfo;      // Line containing values for Qmass, dt, temperature for output to header

  double tempNow;

};
#endif
