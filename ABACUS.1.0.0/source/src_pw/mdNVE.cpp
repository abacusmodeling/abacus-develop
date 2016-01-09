#include"mdNVE.h"

/*
using namespace std;
class mdNVE:public md{
public:
	mdNVE(int n=1):md(1){};
        void	runNVE();
	void RescaleVelocity(double KE,double newKE);
	double Conserved(double KE, double PE);
	double MAXVALF();
private:
	                      //>> INTERNAL VARIABLES <<//
  int nfrozen;               // Number of frozen atoms
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
*/

void mdNVE::runNVE(int step1){
//-------------------------------------------------------------------------------
//
// REFERENCES:
//   Martyna et. al., "Explicit reversible integrators for extended system
//      dynamics," Mol. Phys. 87 (1996), 1117.
//   Allen and Tildesley, "Computer Simulation of Liquids" (predictor-corrector)
//   Frenkel and Smit, "Understanding Molecular Simulation".
//
//------------------------------------------------------------------------------
// REVISION LOG:
//     04-17-2013: revised from NVT code.
//
//------------------------------------------------------------------------------



 // Optimizer // A subroutine that is called to optimize the electron density relative to the ion positions

  
   // rho              // Electron density, realspace
   // energy                    // Energy from electronic optimization
  //Vector3<double>  forces;           // Total forces.  
                     // First is ion number, second direction (1,2,3 for x,y,z)

  

                        //>> INITIALIZATION <<//
  TITLE("mdNVE","runNVE");
   timer::tick("mdNVE","runNVE",'C');
  step=step1+step_rst;
  if(step==1)cout<<" (NVE) Start."<<endl;
  else {
      cout<<" (NVE) step: "<<step<<endl;
      }
  //watch2 = TimerStart()
  //watch = TimerStart()
  //numIon = sizeof(cell%ionTable);
  
  cout<<" (NVE) Fronzen Atom Number : "<<nfrozen<<endl;


  // Calculate the number of total MD steps.
 // nStep = int(timeTot/dt);
  //f.zero_out();//forces = 0;

  // Allocate arrays for velocities, cartesian coordinates.
  //velocity=new Vector3<double>[numIon];
  //cartNoWrap=new Vector3<double>[numIon];
  

  // Set up extra output to ion optimizer / MD header
  cout<<"Time interval       : "<< dt*fundamentalTime*1e15<< " (fs)"<<endl;
  cout<<"Target temperature  : "<< temperature/K_BOLTZMAN_AU<< " (K)"<<endl;
  //cout<<"Number of NVE steps : "<< nStep<<endl;

  //if (rankGlobal==0)  GeometryMinimizerReportHeader(parameterInfo);
  if(step==1){
      for(k=0;k<numIon;k++){
        if(ionmbl[k].x==0)vel[k].x=0;
        if(ionmbl[k].y==0)vel[k].y=0;
        if(ionmbl[k].z==0)vel[k].z=0;
      }
      RemoveMovementOfCenterOfMass();
      scalevel();
  }
  moveatoms(step);
  twiceKE=GetAtomKE();
  twiceKE = twiceKE * 2;

  tempNow = twiceKE/(3*double(numIon-nfrozen))/K_BOLTZMAN_AU;
  cout<<" start temperature = "<< tempNow/K_BOLTZMAN_AU<< " (k)"<<endl;





  // Set up forces 
  callforce();
  for(k=0;k<numIon;k++){
        if(ionmbl[k].x==0)force[k].x=0;
        if(ionmbl[k].y==0)force[k].y=0;
        if(ionmbl[k].z==0)force[k].z=0;
  }

  //cout<<"begin maxForce"<<endl;
  maxForce = MAXVALF();
  cout<<"maxForce: "<<sqrt(maxForce)<<endl; 
  mdenergy=en.etot/2;
  conservedE = Conserved(twiceKE/2, mdenergy);
 // oldEtot=mdenergy;

  // Print intial output
  //GeometryMinimizerReportSteps(startStep, TimerStop(watch), forces, maxForce, &
  //                                  0, .TRUE., hamiltonian, .TRUE., &
  //                                  temperature=temperature);


  cout<< "NVE_STEP "<<" "<<"SystemEnergy"<<" "<< "Conserved"<<" "<< "DeltaE"<<" "<< "Temperature"<<endl;

  // Output the message to the screen.
  cout<<step<<" "<< mdenergy<<" "<< conservedE<<" "<< mdenergy-oldEtot<<" "<< twiceKE/(3*double(numIon-nfrozen))/K_BOLTZMAN_AU<<endl;
 // oldEtot=mdenergy;
                        //>> FUNCTION BODY <<//

  //xiaohui move this line 2015-09-15
  //cout<<"ready part is ready!"<<" "<<endl;

  // Loop for MD step 
 // for(step = startStep+1;step< nStep+1;step++){

    cout<<"------------------------------------------------------------------------------"<<endl;
    cout<< "MD(NVE) STEP "<< step<<endl;
    cout<<"------------------------------------------------------------------------------"<<endl;
    if (!MY_RANK){
         out<<"------------------------------------------------------------------------------"<<endl;
         out<< "MD(NVE) STEP "<< step<<endl;
         out<<"------------------------------------------------------------------------------"<<endl;
         ofs_running<<"------------------------------------------------------------------------------"<<endl;
         ofs_running<< "MD(NVE) STEP "<< step<<endl;
         ofs_running<<"------------------------------------------------------------------------------"<<endl;
    }
   // callforce();
  /*  cout<<"force to test:"<<endl;
    cout<<force[0].x<<" "<<force[0].y<<" "<<force[0].z<<endl;
    cout<<force[1].x<<" "<<force[1].y<<" "<<force[1].z<<endl;

    cout<<"coord to test:"<<endl;
    cout<<cartNoWrap[0].x<<" "<<cartNoWrap[0].y<<" "<<cartNoWrap[0].z<<endl;
    cout<<cartNoWrap[1].x<<" "<<cartNoWrap[1].y<<" "<<cartNoWrap[1].z<<endl;
    // Calculate the Mean-Square-Displacement.
    cout<<"vel to test"<<endl;
    cout<<vel[0].x<<" "<<vel[0].y<<" "<<vel[0].z<<endl;
    cout<<vel[1].x<<" "<<vel[1].y<<" "<<vel[1].z<<endl;*/
    //MonitorMeanSquareDisplacement(step, msd,diffuCoeff);
  
    // (1) 1st step of Verlet-Velocity
    // New velocity are obtained
	for(k=0;k<numIon;k++){       
           mass = allmass[k];
     //2014/6/28
           vel[k] = vel[k] + force[k]/mass*dt/2.0;
       	} 
    // (2) 2nd step of Verlet-Velocity 
    // Update the Non-Wrapped cartesion coordinate
       twiceKE = GetAtomKE();
       twiceKE = 2 * twiceKE;
       if(step!=1||rstMD==1)for(k=0;k<numIon;k++){       
           mass = allmass[k];
     //2014/6/28
           vel[k] = vel[k] + force[k]/mass*dt/2.0;
       	} 
        RemoveMovementOfCenterOfMass();  
       for(i=0;i<numIon;i++){
           cartNoWrap[i] = vel[i]*dt + cartNoWrap[i];
		
	}
    // Calculate the maximal velocities.
    // The step in fractional coordinates
       maxStep = 0;
       for( i = 0;i< numIon;i++){
		if((pow(vel[i].x,2.0)+pow(vel[i].y,2.0)+pow(vel[i].z,2.0))>maxStep)
      maxStep = pow(vel[i].x,2.0)+pow(vel[i].y,2.0)+pow(vel[i].z,2.0);
      fracStep = ionlatvec.Inverse()* vel[i]*dt/ucell.lat0;
      taudirac[i] = taudirac[i] + fracStep;
    }
   // cout<<"dx: "<<fracStep.x<<" "<<fracStep.y<<" "<<fracStep.z<<endl;
    moveatoms(step);
    string t("md_pos_");
    t=intTurnTostring(step,t);
	connection2();
    printpos(t,step);
    maxStep = sqrt(maxStep)*dt;
/*    if(step1){
         PDF(step1);
         printRDF(step1);
    }*/

    // (3.1) OFDFT optimizer: in order to calculate the force
    // Refresh the ionpositions, including the Ewald term
    // and the local pseudopotentials term.
   // RefreshIonPositions(grids); // This must be called every time the 
                                          // ion positions are changed

    // (3.2)
    // Get the new charge density in the new ion positions.

    // (3.3) Calculate the forces. 
    //CalculateForces(rho, forces);
	//xiaohui add 2014-06-13
	//self_consistent(step);
/*	Force_LCAO FL;
	FL.allocate (); 
	FL.start_force();

	int ion;
	for(ion=0;ion<numIon;ion++){
			force[ion].x=FL.fcs(ion,1)/2;
			force[ion].y=FL.fcs(ion,2)/2;
			force[ion].z=FL.fcs(ion,3)/2;
	}

    //callforce();
	for(i=0;i<numIon;i++){
        if(ionmbl[k].x==0)force[k].x=0;
        if(ionmbl[k].y==0)force[k].y=0;
        if(ionmbl[k].z==0)force[k].z=0;
	}
    maxForce = MAXVALF();*/
	// (3.4) 3rd step to Update velocity
  //2014/6/28
  /*  for( ii=0;ii<numIon;ii++){        
      mass =  allmass[ii];
      vel[ii] = vel[ii] + force[ii]/mass*dt/2.0;	  
    }*/

    // calculate the kinetic energy of ions

    //------------------------------------------------------
    // this part should be not necessary in NVE
    // mohan 2013-05-09
    // calculate the conserved quantity during MD 
    //hamiltonian = Conserved(twiceKE/2._DP, energy(1))
    //newKE = 0.5*twiceKE-(hamiltonian-conservedE) 
    //newKE = 2.0*newKE;

    // updates velocity in order to adjust temperature
    //CALL RescaleVelocity(twiceKE, newKE, numIon, nfrozen)
    //------------------------------------------------------

    // calculate the conserved quantity during MD 
    hamiltonian = Conserved(twiceKE/2, mdenergy);

// kinetic, external, coulombic, exchange-correlation,
// ion-ion, Thomas-Fermi, von Weiszacker and
// third term Wang-Teter, WGC,



    // For monitoring MD or restarting MD
  //   OutputMDGeom(step);
     cout<< setprecision (9)<<hamiltonian<<" "<< setprecision (9)<<twiceKE/2<<endl;
     if (doMSD ==1){ 
          OutputIonVelocityAndStats(step);
     }
     MonitorMeanSquareDisplacement(step);
    // output the density
    //OutputDen(step, grids(0)%rho);




    // Print output
    //GeometryMinimizerReportSteps(step, TimerStop(watch), &
    //                                  forces, maxForce, maxStep, &
    //                                  .TRUE., hamiltonian, .TRUE., &
     //                                 temperature=twiceKE/(3*double(numIon-nfrozen)))

   // watch = TimerStart();
    
    // Output the message to the screen.
    if (!MY_RANK){ 
         out<<step<<" "<< mdenergy<<" "<< hamiltonian<<" "<< mdenergy-oldEtot<<" "<< twiceKE/(3*double(numIon-nfrozen))/K_BOLTZMAN_AU<<endl;
         out<<step<<" "<< mdenergy<<" "<< conservedE<<" "<< mdenergy-oldEtot<<" "<< twiceKE/(3*double(numIon-nfrozen))/K_BOLTZMAN_AU<<endl;
         ofs_running<<step<<" "<< mdenergy<<" "<< hamiltonian<<" "<< mdenergy-oldEtot<<" "<< twiceKE/(3*double(numIon-nfrozen))/K_BOLTZMAN_AU<<endl;
         ofs_running<<step<<" "<< mdenergy<<" "<< conservedE<<" "<< mdenergy-oldEtot<<" "<< twiceKE/(3*double(numIon-nfrozen))/K_BOLTZMAN_AU<<endl;
    }
    oldEtot=mdenergy;
    mstout(step);

//	}
  //----------------------------------------
  // End of NVE MD loop
  //----------------------------------------

  //GeometryMinimizerReportFooter(TimerStop(watch2));

  //cout<<"(NVE):this step finished."<<endl;
/*
  if(force!=NULL)delete(force);
  if(ionmbl!=NULL)delete(ionmbl);
  if(cartNoWrap!=NULL)delete(cartNoWrap);
  if(vel!=NULL)delete(vel);
 // if(msdcoords0!=NULL) delete(msdcoords0);
 // if(msdcoordst!=NULL) delete(msdcoordst);
  if(taudirac!=NULL)delete(taudirac);
  if(allmass!=NULL)delete(allmass);*/
     cout<<"(NVE): this step finished."<<endl;
     if (!MY_RANK){
         out<<"(NVE): this step finished."<<endl;
         ofs_running<<"(NVE): this step finished."<<endl;
     }
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

     timer::tick("mdNVE","runNVE",'C');
     return;
}


void mdNVE::RescaleVelocity(double KE,double newKE){
//---------------------------------------------------------------------------
// DESCRIPTION:
//   This function propagates the Nose-Hoover thermostat extended-system
//   variables
//----------------------------------------------------------------------------
                        //>> INTERNAL VARIABLES <<//

  double scaling,tempNow, tempTarget;


                        //>> INITIALIZATION <<//

  tempNow = KE/(3*double(numIon-nfrozen))/K_BOLTZMAN_AU;
  tempTarget = newKE/(3*double(numIon-nfrozen))/K_BOLTZMAN_AU;

  out<< "Temperature now     : "<< tempNow<< " (K)"<<endl;
  out<< "New temperature  : "<< temperature/K_BOLTZMAN_AU<< " (K)"<<endl;

  scaling = sqrt(tempTarget/tempNow);
                        //>> FUNCTION BODY <<//

  out<<" scaling factor of velocity : "<< scaling<<endl;
 for(int i=0;i<numIon;i++) vel[i] = vel[i]*scaling;
  twiceKE = twiceKE*scaling*scaling;
  return;
}

double mdNVE::Conserved(double KE, double PE){
//---------------------------------------------------------------------------
// DESCRIPTION:
//   This function calculates the conserved quantity for the NVE system. 
//----------------------------------------------------------------------------

   // KE   // Kinetic energy of particles (due to their velocity)
   // PE   // Potential energy of particles (energy calculated using OFDFT)

  double Conserved; // The conserved quantity

                        //>> INITIALIZATION <<//
                        //>> FUNCTION BODY <<//

   Conserved = KE + PE ;
   
   if (!MY_RANK){                 
       out<< "NVE Conservation     : "<< Conserved<<" (Hartree)"<<endl;
       out<< "NVE Temperature      : "<< twiceKE/(3*double(numIon-nfrozen))/K_BOLTZMAN_AU<<" (K)"<<endl;
       out<< "NVE Kinetic energy   : "<< KE<<" (Hartree)"<<endl;
       out<< "NVE Potential energy : "<< PE<<" (Hartree)"<<endl;
       ofs_running<< "NVE Conservation     : "<< Conserved<<" (Hartree)"<<endl;
       ofs_running<< "NVE Temperature      : "<< twiceKE/(3*double(numIon-nfrozen))/K_BOLTZMAN_AU<<" (K)"<<endl;
       ofs_running<< "NVE Kinetic energy   : "<< KE<<" (Hartree)"<<endl;
       ofs_running<< "NVE Potential energy : "<< PE<<" (Hartree)"<<endl;
   }
   return Conserved;
}
double mdNVE::MAXVALF(){
	//cout<<"enter in MAXVALF"<<endl;
	double max,force0;
	int i,j;
	max=0;
	for(i=0;i<numIon;i++){
		force0=0;
		force0+=pow(force[i].x,2)+pow(force[i].y,2)+pow(force[i].z,2);
		if(max<force0)max=force0;
	}
	return max;
}

