#ifndef MDNVT_H
#define MDNVT_H

#include "md.h"

using namespace std;
class mdnvt:public md{
public:	mdnvt(int n=1):md(1){};
    void runnvt(int step1);
	void NHIntegrator();
	double NHhamiltonian(double KE,double PE);
	double MAXVALF();
private:
//	int  nfrozen;//                 ! Number of frozen atoms
  int  ii,k,it,ion;//
  int  i;//                       ! Counters to parse through the ion table.
  int  step;//                    ! step number we are on
  //int  numIon;//                  ! total number of atoms
  int  nStep;//                   ! Total number of steps
  int  errorFlag;  
  double  mass;            // atom mass, temperay var
  double  maxStep;       // Largest movement by an ion
  //double  xLogS;         // position of thermostat
  //double  vLogS;         // vel of thermostat
  double  hamiltonian;   // conserved energy
  double  maxForce;      
  //double  msd;           // mean square displacement
  //double  diffuCoeff;      // diffusion coefficient
  double  twiceKE;         // Kinetic energy x 2
  double  oldEtot;           // old energy
};
#endif
