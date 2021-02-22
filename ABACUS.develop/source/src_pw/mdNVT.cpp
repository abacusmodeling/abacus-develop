#include "mdNVT.h"

const double EXP=2.71828;
/*
using namespace std;
class mdnvt:public md{
public:	mdnvt(int n=1):md(1){};
    void runnvt();
	void NHIntegrator();
	double NHhamiltonian(double KE,double PE);
	double MAXVALF();
private:
	int  nfrozen;//                 ! Number of frozen atoms
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
*/
void mdnvt::runnvt(int step1){
//------------------------------------------------------------------------------
// DESCRIPTION:
//  This subroutine integrates motion for molecular dynamics using the
//  Nose-Hoover thermostat coupled with a predictor-corrector integrator.
// REFERENCES:
//  Martyna et. al., "Explicit reversible integrators for extended system
//     dynamics," Mol. Phys. 87 (1996), 1117.
//  Allen and Tildesley, "Computer Simulation of Liquids" (predictor-corrector)
//  Frenkel and Smit, "Understanding Molecular Simulation".
//------------------------------------------------------------------------------
//REVISION LOG:
//    10/20/2008 : Adapted from subroutines in IonOptimizers and symplectic
//                 algorithm from Martyna et. al. with help from Frenkel/Smit.
//
//------------------------------------------------------------------------------
//  Vector4 rho;   // Electron density, realspace
//  Vector energy;      // Energy from electronic optimization
  //Vector3 forces;// Total forces.  Final index is 0 for total force,
  // First is ion number, second direction (1,2,3 for x,y,z)

  
  //REAL(KIND=DP), ALLOCATABLE, DIMENSION(:,:) :: vel // vel of ions (atomic units)
  //double intCoeff[];
  //REAL(KIND=DP), DIMENSION(3) :: fracStep // Step (due to vel) in fractional coordinates

  //string message;
  //string parameterInfo;      // Line containing values for Qmass, dt, temperature
                                          // for output to header


                      //  !>> INITIALIZATION <<!
  TITLE("mdnvt","runnvt");
   timer::tick("mdnvt","runnvt",'C');
  step=step1+step_rst;
  if (!MY_RANK)out<<"step: " << step <<endl;
  if (step!=1)ReadNewTemp( step );
  intCoeff=new double[nYosh];

  //xiaohui modify 2015-09-15
  //if(step==1)cout<<" (NVT) Start."<<endl;
  //else {
  //       cout<<" (NVT) step: "<<step<<endl;
  //}
  //// Calculate # of frozen atoms 
  //cout<<" (NVT) Fronzen Atom Number : "<<nfrozen<<endl;


  MakeIntCoeff();

  //xiaohui modify 2015-09-15
  //// Set up extra output to ion optimizer / MD header
  //cout<<" Qmass for NVT       : "<< Qmass<<" (a.u.)"<<endl;
  //cout<< " Time interval       : "<< dt*fundamentalTime*1e15<< " (fs)"<<endl;

  //xiaohui add 2015-09-15
  cout << " ---------------------------------------------------------" << endl;

  cout<< " Target temperature  : "<< temperature/K_BOLTZMAN_AU<< " (K)"<<endl;
 // cout<< " Number of NVT steps : "<< nStep<<endl;

  //if (rankGlobal==0)  GeometryMinimizerReportHeader(parameterInfo);

  if(step==1||step%fixTemperature==1){
      for(k=0;k<numIon;k++){
        if(ionmbl[k].x==0)vel[k].x=0;
        if(ionmbl[k].y==0)vel[k].y=0;
        if(ionmbl[k].z==0)vel[k].z=0;
      }
      if(nfrozen==0)RemoveMovementOfCenterOfMass(); //added 20170626
      scalevel();
  }

  // get the kinetic energy
  twiceKE = GetAtomKE();
  twiceKE = twiceKE * 2;
 /* vLogS=0;
  xLogS=0;*/
   //Set up forces 
  //Optimizer;                        //Optimize the density
  callforce();
	if(STRESS)
	{
		cout<<"output Pressure for check!"<<endl;
		double press;
		for(int i=0;i<3;i++)
			press += stress_lcao(i,i)/3;
		press += twiceKE/3/ucell.omega; //output virtual press = 2/3 *Ek/V + sum(sigma[i][i])/3
		double unit_transform = RYDBERG_SI / pow(BOHR_RADIUS_SI,3) * 1.0e-8 ;
		cout<<"Virtual Pressure is "<<press*unit_transform<<" Kbar "<<endl;
	}
  for(k=0;k<numIon;k++){
     if(ionmbl[k].x==0)force[k].x=0;
     if(ionmbl[k].y==0)force[k].y=0;
     if(ionmbl[k].z==0)force[k].z=0;
  }
  maxForce = MAXVALF();
  mdenergy=en.etot/2;
  hamiltonian = NHhamiltonian(twiceKE/2, mdenergy);
  //xiaohui move this line 2015-09-15
  //cout<<"ready part is ready!"<<endl;

  // Print intial output
  //GeometryMinimizerReportSteps(startStep, TimerStop(watch), forces, maxForce,0,1, hamiltonian,1,temperature)
  //----------------------------------------------
  // big loop
  //-----------------------------------------------
  cout<<" "<<std::left<<setw(12)<<"MD_STEP"<<std::left<<setw(12)<< "SystemE"<<std::left<<setw(12)<< "Conserved"<<std::left<<setw(12)<< "DeltaE"<<std::left<<setw(12)<< "Temperature"<<endl;

  cout<<" "<<std::left<<setw(12)<<step<<std::left<<setw(12)<< mdenergy<<std::left<<setw(12)<< hamiltonian<<std::left<<setw(12)<< mdenergy-oldEtot<<std::left<<setw(12)<<twiceKE/(3*double(numIon-nfrozen))/K_BOLTZMAN_AU<<endl;
 
  oldEtot=mdenergy;

 
  // Loop for MD step 

    /*xiaohui move 2015-09-15
    //cout<<" --------------------------------------------------"<<endl;
    cout << " ---------------------------------------------------------" << endl;
    cout<<" Molecular Dynamics (NVT) STEP "<< step<<endl;
    */

    //cout<< " --------------------------------------------------"<<endl;
    //cout << " ---------------------------------------------------------" << endl;
    if (!MY_RANK){
          out<<" --------------------------------------------------"<<endl;
          out<<" Molecular Dynamics (NVT) STEP "<< step<<endl;
          out<< " --------------------------------------------------"<<endl;
          ofs_running<<" --------------------------------------------------"<<endl;
          ofs_running<<" Molecular Dynamics (NVT) STEP "<< step<<endl;
          ofs_running<< " --------------------------------------------------"<<endl;
    }
    // (1) Calculate the Mean-Square-Displacement.
    if(step==1&&rstMD==0){ 
    // (2) First themorstat step: updates vel, twiceKE, xLogS, vLogS
   // MonitorMeanSquareDisplacement(step);
         NHIntegrator();
       
    // (3) New vel obtained (Verlet-like step)
        for( ii=0;ii<numIon;ii++){ 
            mass = allmass[ii];
            vel[ii] = vel[ii] + force[ii]/mass*dt/2.0;
	}
    }
    else {
        for( ii=0;ii<numIon;ii++){ 
            mass = allmass[ii];
            vel[ii] = vel[ii] + force[ii]/mass*dt/2.0;
        }
        twiceKE=GetAtomKE();
        twiceKE = 2 * twiceKE;
        NHIntegrator();
        twiceKE=GetAtomKE();
        twiceKE = 2 * twiceKE;
        hamiltonian = NHhamiltonian(twiceKE/2, mdenergy);
        NHIntegrator();
        for( ii=0;ii<numIon;ii++){ 
            mass = allmass[ii];
            vel[ii] = vel[ii] + force[ii]/mass*dt/2.0;
        }
    }
    // (4) Update the Non-Wrapped cartesion coordinates
        if(nfrozen==0)RemoveMovementOfCenterOfMass();  
        if(mdtype==2)scalevel();//choose Nose Hoover plus velocity scaling method
	for( ii=0;ii<numIon;ii++){ 
            cartNoWrap[ii] = vel[ii]*dt +cartNoWrap[ii];
	}
    
        MonitorMeanSquareDisplacement(step);
    // Calculate the maximal velocities.
    // The step in fractional coordinates
    maxStep = 0;
    Vector3<double> fracStep;
   	for( ii=0;ii<numIon;ii++){ 
	    	if((pow(vel[ii].x,2)+pow(vel[ii].y,2)+pow(vel[ii].z,2))>maxStep)
                maxStep = pow(vel[ii].x,2)+pow(vel[ii].y,2)+pow(vel[ii].z,2);
                Mathzone::Cartesian_to_Direct(vel[ii].x*dt/ucell.lat0,vel[ii].y*dt/ucell.lat0,vel[ii].z*dt/ucell.lat0,
                                    ionlatvec.e11,ionlatvec.e12,ionlatvec.e13,
                                    ionlatvec.e21,ionlatvec.e22,ionlatvec.e23,
                                    ionlatvec.e31,ionlatvec.e32,ionlatvec.e33,
                                    fracStep.x,fracStep.y,fracStep.z);

//            fracStep = ionlatvec.Inverse()*vel[ii]*dt/ucell.lat0; 
            taudirac[ii] = taudirac[ii] + fracStep;
	}

    moveatoms(step);
    string t("md_pos_");
    t=intTurnTostring(step,t);
    connection2();
    printpos(t,step);
    maxStep = sqrt(maxStep)*dt;
/*   if(step1){
         PDF(step1);
         printRDF(step1);
    }*/

    // (5) OFDFT optimizer
    // Refresh the ionpositions, including the Ewald term
    // and the local pseudopotentials term.
    //RefreshIonPositions(grids); // This must be called every time the 
                                    // ion positions are changed

    // Get the new charge density in the new ion positions.
     //Optimizer;                        //Optimize the density
    //CalculateStress( energy, stress);
    //PrintStress(stress);


    // (8)


    
    // (9) Second thermostat step: updates vel, twiceKE, xLogS, vLogS
// cout<<"testNVT: x "<<xLogS<<" v "<<vLogS<<endl;

    // For monitoring MD or restarting MD   
  //   OutputMDGeom(step);

     //xiaohui move this line, 2015-09-15
     //cout<<" ekin " << setprecision (9)<<twiceKE/2<<" etot "<< setprecision (9)<<hamiltonian<<endl;

     if (doMSD ==1){
     OutputIonVelocityAndStats(step);
}

    //OutputDen(step, grids(0)%rho);

    // (10) Conserved quantity during MD 

// kinetic, external, coulombic, exchange-correlation,
// ion-ion, Thomas-Fermi, von Weiszacker and
// third term Wang-Teter, WGC,
    if (!MY_RANK){
         out<<setw(15)<<"maxForce=     "<<setw(15)<<"maxstep=      "<<setw(15)<<"step=     "<<endl;
         out<<setw(15)<<maxForce<<setw(15)<<maxStep<<setw(15)<<step<<endl;

    //watch = TimerStart()
    
    // Change temperature if needed

    // Output the message to the screen.
         out<<step<<" "<< mdenergy<<" "<< hamiltonian<<" "<< mdenergy-oldEtot<<" "<< twiceKE/(3*double(numIon-nfrozen))/K_BOLTZMAN_AU<<endl;
    }
    oldEtot=mdenergy;
    mstout(step);
    double thisRDF=0;
    if(step%fixTemperature==0){
        for(i=42;i<48;i++){
           thisRDF+=rdf[i];
        }
  //      cout<<"RDFthirdpeak "<<step/fixTemperature<<" "<<thisRDF<<endl;
    }

 //----------------------------------------
  // End of NVT MD loop
  //----------------------------------------

  // GeometryMinimizerReportFooter(TimerStop(watch2));

  //xiaohui move these two lines, 2015-09-15
  //cout<<"(NVT): finished."<<endl;
  //if (!MY_RANK)  ofs_running<<"(NVT): finished."<<endl;

  //delete(cartNoWrap);
  //delete(vel);
  delete []intCoeff;
//2015-09-25, xiaohui
#ifdef __MPI
  if(step%100==0){
    MPI_Bcast(vel,numIon*3,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(cartNoWrap,numIon*3,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(taudirac,numIon*3,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&xLogS,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&vLogS,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  }
#endif
  //IF (ALLOCATED(msd_coords0)) DEALLOCATE(msd_coords0)
  //IF (ALLOCATED(msd_coordst)) DEALLOCATE(msd_coordst)
   timer::tick("mdnvt","runnvt",'C');
  return;
}


void mdnvt::NHIntegrator(){
//---------------------------------------------------------------------------
// DESCRIPTION:
//   This function propagates the Nose-Hoover thermostat extended-system
//   variables using multiple step technique.
//   Note that only thermostat variables are updated, including
//   thermostat positions (xLogS),
//   thermostat velocities (vLogS),
//   thermostat acceleration (gLogS),
//   but not including any variables for atoms. 
//--------------------------------------------------------------------

  double freedom,scaling,gLogS;

  int iResn, iYosh;
  double wdt2, wdt4;

                    //    !>> INITIALIZATION <<!

  

  freedom = 3*double(numIon-nfrozen);
  scaling = 1.0;

  // calculate the force of thermostat
  gLogS = (twiceKE - freedom*temperature) / Qmass;
//  cout<<"testNVT: x "<<xLogS<<" v "<<vLogS<<" g "<<gLogS<<" ekin "<<twiceKE<<" temperature " <<temperature/K_BOLTZMAN_AU<<endl;
                        //!>> FUNCTION BODY <<!

//  WRITE(outputUnit,*) " freedom=",freedom
//  WRITE(outputUnit,*) " gLogS=",gLogS
//  WRITE(outputUnit,*) " twiceKE=",twiceKE
//  WRITE(outputUnit,*) " temperature=",temperature
  
  for( iResn = 0;iResn< nResn;iResn++){
	  for( iYosh = 0;iYosh< nYosh;iYosh++){
      wdt2 = intCoeff[iYosh]*dt/(2*double(nResn));
      wdt4 = intCoeff[iYosh]*dt/(4*double(nResn));

      // calculate the first half of integration.
      vLogS = vLogS + gLogS*wdt4;
      scaling = scaling * pow(EXP,(-wdt2*vLogS));

      // In some cases, the exponential part will be NaN,
      // Try to give warning at this case, added by mohan
     // if( check_isnan(scaling) ){
      //  out<< " iResn="<< iResn<< " iYosh="<< iYosh<<endl;
      //  out<< " gLogS="<< gLogS<<endl;
      //  out<< " vLogS="<< vLogS<<endl;
      //  out<<"Scaling in NHIntegrator is NaN !"<<endl;
      //  STOP
      //}

      gLogS = (scaling*scaling*twiceKE - freedom*temperature) / Qmass;
      // Calculate xLogS only as a check of code
      xLogS = xLogS + vLogS*wdt2;

      // calculate the second half of integration
      vLogS = vLogS + gLogS*wdt4;
	  }
  }
  
  if (!MY_RANK)out<< " scaling factor of vel : "<< scaling<<endl;
  for(i=0;i<numIon;i++){
    vel[i] = vel[i]*scaling;
  }
  twiceKE = twiceKE*scaling*scaling;
  return;
}


double mdnvt::NHhamiltonian(double KE,double PE){
//---------------------------------------------------------------------------
// DESCRIPTION:
//   This function calculates the conserved quantity for the Nose-Hoover
//   system.
//----------------------------------------------------------------------------
 

//  REAL(kind=DP), INTENT(IN) :: &
//    KE, &    ! Kinetic energy of particles (due to their vel)
//    PE, &    ! Potential energy of particles (energy calculated using OFDFT)
//    xLogS, & ! "position" of constraint variable
//    vLogS    ! "vel" of constraint variable

  double NHhamiltonian0; // The conserved quantity

                       // !>> INTERIOR VARIABLES <<!

  NHhamiltonian0 = KE + PE + 0.5*vLogS*vLogS*Qmass +(3*double(numIon-nfrozen))*temperature*xLogS;
   
   if (!MY_RANK){
   out<< " "<<endl;
   out<< " --------------------------------------------------"<<endl;
   out<< " SUMMARY OF NVT CALCULATION"<<endl;
   out<<" --------------------------------------------------"<<endl;
   out<<" NVT Conservation     : "<<setw(10)<< NHhamiltonian0*2<<" (Rydberg)"<<endl;
   out<<" NVT Temperature      : "<<setw(10)<< twiceKE/(3*double(numIon-nfrozen))/K_BOLTZMAN_AU<<" (K)"<<endl;
   out<<" NVT Kinetic energy   : "<<setw(10)<< KE*2<<" (Rydberg)"<<endl;
   out<<" NVT Potential energy : "<<setw(10)<< PE*2<<" (Rydberg)"<<endl;
   out<< " Thermostat Position : "<<setw(10)<< xLogS<<endl;
   out<< " Thermostat vel : "<< vLogS<<endl;
   ofs_running<< " "<<endl;
   ofs_running<< " --------------------------------------------------"<<endl;
   ofs_running<< " SUMMARY OF NVT CALCULATION"<<endl;
   ofs_running<<" --------------------------------------------------"<<endl;
   ofs_running<<" NVT Conservation     : "<<setw(10)<< NHhamiltonian0*2<<" (Rydberg)"<<endl;
   ofs_running<<" NVT Temperature      : "<<setw(10)<< twiceKE/(3*double(numIon-nfrozen))/K_BOLTZMAN_AU<<" (K)"<<endl;
   ofs_running<<" NVT Kinetic energy   : "<<setw(10)<< KE*2<<" (Rydberg)"<<endl;
   ofs_running<<" NVT Potential energy : "<<setw(10)<< PE*2<<" (Rydberg)"<<endl;
   ofs_running<< " Thermostat Position : "<<setw(10)<< xLogS<<endl;
   ofs_running<< " Thermostat vel : "<< vLogS<<endl;

   }
   return NHhamiltonian0;
}

// check 
//LOGICAL FUNCTION check_isnan(a) 
//REAL(KIND=DP) a 
//if (a.ne.a) then 
//check_isnan = .true. 
//else 
//check_isnan = .false. 
//end if 
//return 
//END FUNCTION check_isnan
double mdnvt::MAXVALF(){
	double max,force0;
	int i,j;
	max=0;
	for(i=0;i<numIon;i++){
		force0=0;
		for(j=0;j<3;j++){
			force0+=pow(force[i].x,2)+pow(force[i].y,2)+pow(force[i].z,2);
		}
		if(max<force0)max=force0;
	}
	return max;
}

