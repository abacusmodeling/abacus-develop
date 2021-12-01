#include "MD_basic.h"
#include "../input.h"

//define in MD_basic.h
//class MD_basic

MD_basic::MD_basic(MD_parameters& MD_para_in, UnitCell_pseudo &unit_in):
    mdp(MD_para_in),
    ucell(unit_in)
{	
	mdp.dt=mdp.dt/ModuleBase::AU_to_FS;
	temperature_=mdp.tfirst/ModuleBase::Hartree_to_K;


    // CMD parameters
    mdp.rcut_lj *= ModuleBase::ANGSTROM_AU;
	mdp.epsilon_lj /= ModuleBase::Hartree_to_eV;
	mdp.sigma_lj *= ModuleBase::ANGSTROM_AU;


	//other parameter default
	outputstressperiod_ = 1 ;// The period to output stress
	step_rst_=0;
    step_=0;

    vel=new ModuleBase::Vector3<double>[ucell.nat]; 
	cart_change=new ModuleBase::Vector3<double>[ucell.nat];
	tauDirectChange=new ModuleBase::Vector3<double>[ucell.nat];
	allmass=new double[ucell.nat];
	ionmbl=new ModuleBase::Vector3<int>[ucell.nat];

    frozen_freedom_ = MD_func::getMassMbl(ucell, mdp, allmass, ionmbl);
    MD_func::InitVel(ucell, temperature_, allmass, frozen_freedom_, ionmbl, vel);
    // if (ucell.set_vel)    //  Yuanbo Li 2021-08-01
    // {
    //     int iat=0;    //initialize velocity of atoms from STRU  liuyu 2021-07-14
    //     for(int it=0; it<ucell.ntype; ++it)
    //     {
    //         for(int ia=0; ia<ucell.atoms[it].na; ++ia)
    //         {
    //             vel[iat] = ucell.atoms[it].vel[ia];
    //             ++iat;
    //         }
    //     }
    //     assert(iat==ucell.nat);
    // }
    // else
    // {
    //     MD_func::InitVelocity(ucell.nat, temperature_, fundamentalTime, allmass, vel);
    // }

    //MD starting setup
    if(mdp.rstMD)
    {
        if(!MD_func::RestartMD(ucell.nat, vel, step_rst_))
        {
			std::cout<<"error in restart MD!"<<std::endl;
			return;
		}
    }
/*	if(!mdp.rstMD){
		frozen_freedom_ = MD_func::getMassMbl(ucell, allmass, ionmbl);
		//MD_func::gettauDirectChange(ucell, tauDirectChange);
		//A fresh new MD: Do not restart MD
		MD_func::InitVelocity(ucell.nat, temperature_, fundamentalTime, allmass, vel) ;
		// Initialize thermostat, and barostat
		
	}
	else{
		frozen_freedom_ = MD_func::getMassMbl(ucell, allmass, ionmbl);
		//MD_func::gettauDirectChange(ucell, tauDirectChange);
		MD_func::InitVelocity(ucell.nat, temperature_, fundamentalTime, allmass, vel) ;
		if(!MD_func::RestartMD(ucell.nat, vel, step_rst_)){
			std::cout<<"error in restart MD!"<<std::endl;
			return;
		}
	}*/
	if(mdp.NVT_tau<1e-10)
	{
		if(mdp.NVT_control == 1) mdp.NVT_tau = 20 * mdp.dt; //NHC
		else if(mdp.NVT_control == 2) mdp.NVT_tau = 20 * mdp.dt; //LGV
		else if(mdp.NVT_control == 3) mdp.NVT_tau = 20 * mdp.dt; //ADS
		else ModuleBase::WARNING_QUIT("md:InitMD", "please choose available reservoir!!!");
	}
	else
	{
		mdp.NVT_tau /= ModuleBase::AU_to_FS;
	}
	//tau = 1.0/(NVT_tau*fundamentalTime*1e15);
	// Q=KbT*tau**2
        
    mdp.Qmass=mdp.Qmass/ModuleBase::AU_to_MASS;
	if(mdp.Qmass<1e-10)
    {
        mdp.Qmass = pow(mdp.NVT_tau,2)*temperature_;///beta;
    }
	
	//NHC setup
	if(mdp.mdtype==1 && mdp.NVT_control == 1)
	{
		mdt.init_NHC(
        mdp.MNHC, 
        mdp.Qmass, 
        mdp.NVT_tau, 
        mdp.dt,
        mdp.NVT_control, 
        GlobalV::ofs_running, 
        ucell.nat,
        temperature_,
        vel,
        allmass
        );
        if(mdp.rstMD) mdt.NHC_restart();
	}
}

MD_basic::~MD_basic()
{
    //delete []force;
    delete []ionmbl;
    delete []cart_change;
    delete []vel;
    delete []tauDirectChange;
    delete []allmass;
}

void MD_basic::runNVT(int step1, double potential, ModuleBase::Vector3<double> *force, const ModuleBase::matrix &stress)
{
//------------------------------------------------------------------------------
// DESCRIPTION:
// Molecular dynamics calculation with fixed Volume and slight fluctuated temperature
// Using thermostat : 1, Nose-Hoover Chains; 2, Langevin; 3, Anderson
// Normal Nose-Hoover thermostat method is retained for test.
//------------------------------------------------------------------------------

	ModuleBase::TITLE("MD_basic","runnvt");
	ModuleBase::timer::tick("MD_basic","runnvt");
	step_=step1+step_rst_;
	//the real MD step
	
	//change target temperature
//	if (step_!=1)MD_func::ReadNewTemp( step_ );
	
	std::cout << " ------------------------------------------------------------" << std::endl;
	std::cout << " Target temperature  : " << temperature_*ModuleBase::Hartree_to_K << " (K)"<< std::endl;
	
	if(step_==1||step_%mdp.fixTemperature==1)
    {
		for(int k=0;k<ucell.nat;k++)
        {
			if(ionmbl[k].x==0)vel[k].x=0;
			if(ionmbl[k].y==0)vel[k].y=0;
			if(ionmbl[k].z==0)vel[k].z=0;
		}
		MD_func::scalevel(ucell.nat, frozen_freedom_, temperature_, vel, allmass);
	}

	// get the kinetic energy
	double twiceKE = MD_func::GetAtomKE(ucell.nat, vel, allmass);
	twiceKE = twiceKE * 2;

	//print total stress + stress_MD
	if(GlobalV::STRESS)
	{
		outStressMD(stress, twiceKE);
	}

	for(int k=0;k<ucell.nat;k++){
		if(ionmbl[k].x==0)force[k].x=0;
		if(ionmbl[k].y==0)force[k].y=0;
		if(ionmbl[k].z==0)force[k].z=0;
	}
	double maxForce = MD_func::MAXVALF(ucell.nat, force);

	energy_=potential;

    double hamiltonian;
	//----------------------------------------------
	// big loop
	//-----------------------------------------------
	std::cout<<" "<<std::left<<std::setw(12)<<"MD_STEP"<<std::left<<std::setw(12)<< "SystemE"<<std::left<<std::setw(12)<< "Conserved"<<std::left<<std::setw(12)<< "DeltaE"<<std::left<<std::setw(12)<< "Temperature"<<std::endl;
	std::cout<<" "<<std::left<<std::setw(12)<<step_<<std::left<<std::setw(12)<< energy_<<std::left<<std::setw(12)<< hamiltonian<<std::left<<std::setw(12)<< energy_-oldEtot_<<std::left<<std::setw(12)<<twiceKE/(double(3*ucell.nat-frozen_freedom_))*ModuleBase::Hartree_to_K<<std::endl;
	std::cout << " ------------------------------------------------------------" << std::endl;

	oldEtot_=energy_;

	if (!GlobalV::MY_RANK)
    {
        GlobalV::ofs_running<<"\n --------------------------------------------------"<<std::endl;
        GlobalV::ofs_running<<" Molecular Dynamics (NVT) STEP "<< unsigned(step_) <<std::endl;
        GlobalV::ofs_running<< " --------------------------------------------------"<<std::endl;
	}
	
	// Calculate the Mean-Square-Displacement.
	if(step_==1&&mdp.rstMD==0)
    { 
        //Note: side scheme position before
        //Now turn to middle scheme.
	
	    // New vel obtained (Verlet-like step)
        if(mdp.NVT_control == 1)
        {
            hamiltonian = mdt.NHChamiltonian(
                twiceKE/2, 
                energy_,
                temperature_,
                frozen_freedom_);
        }
        else hamiltonian = MD_func::Conserved(twiceKE/2, energy_, 3*ucell.nat-frozen_freedom_);
	    this->update_half_velocity(force);
	}
	else 
    {
		this->update_half_velocity(force);
		
		twiceKE=MD_func::GetAtomKE(ucell.nat, vel, allmass);
		twiceKE = 2 * twiceKE;
		if(mdp.NVT_control==1)
        {
            hamiltonian = mdt.NHChamiltonian(
                twiceKE/2, 
                energy_,
                temperature_,
                frozen_freedom_);
        } 
		else hamiltonian = MD_func::Conserved(twiceKE/2, energy_, ucell.nat-frozen_freedom_);
        //Note: side scheme position before
        //Now turn to middle scheme. 
		this->update_half_velocity(force);
	}

	// Update the Non-Wrapped cartesion coordinates
    if(mdp.mdtype==2) MD_func::scalevel(ucell.nat, frozen_freedom_, temperature_, vel, allmass);//choose velocity scaling method

    update_half_direct(1);
	
	mdt.Integrator(mdp.NVT_control, temperature_, vel, allmass);//thermostat interact with velocity
    twiceKE=MD_func::GetAtomKE(ucell.nat, vel, allmass);
    twiceKE = 2 * twiceKE;
	
	update_half_direct(0);

	// Calculate the maximal velocities.
	// The step in fractional coordinates
	double maxStep = 0;
    for(int  ii=0;ii<ucell.nat;ii++)
    { 
		if((pow(vel[ii].x,2)+pow(vel[ii].y,2)+pow(vel[ii].z,2))>maxStep)
        {
		    maxStep = pow(vel[ii].x,2)+pow(vel[ii].y,2)+pow(vel[ii].z,2);
        }
    }
	getTaudUpdate();

	//save the atom position change to DFT module
	save_output_position();
	maxStep = sqrt(maxStep)*mdp.dt;

    if (!GlobalV::MY_RANK)
    {
        GlobalV::ofs_running << " maxForce             : " << std::setw(10) << maxForce << std::endl;
        GlobalV::ofs_running << " maxStep              : " << std::setw(10) << maxStep << std::endl;
        GlobalV::ofs_running << " " <<std::left<<std::setw(20)<<"MD_STEP"<<std::left<<std::setw(20)<< "SystemE"<<std::left<<std::setw(20)<< "Conserved"<<std::left<<std::setw(20)<< "DeltaE"<<std::left<<std::setw(20)<< "Temperature"<<std::endl;
	    GlobalV::ofs_running << " " <<std::left<<std::setw(20)<<step_<<std::left<<std::setw(20)<< energy_<<std::left<<std::setw(20)<< hamiltonian<<std::left<<std::setw(20)<< energy_-oldEtot_<<std::left<<std::setw(20)<<twiceKE/(double(3*ucell.nat-frozen_freedom_))/ModuleBase::K_BOLTZMAN_AU<<std::endl;
    }
    oldEtot_=energy_;
    //output basic restart info
    MD_func::mdRestartOut(step_, mdp.recordFreq, ucell.nat, vel);
    //output NHC thermo restart info
    if(mdp.NVT_control==1) mdt.NHC_info_out(step_, mdp.recordFreq, ucell.nat);

//2015-09-25, xiaohui
#ifdef __MPI
    if(step_%100==0)
    {
        MPI_Bcast(vel,ucell.nat*3,MPI_DOUBLE,0,MPI_COMM_WORLD);
    //    MPI_Bcast(&xLogS,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    //    MPI_Bcast(&vLogS,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
#endif
    ModuleBase::timer::tick("MD_basic","runnvt");
    return;
}

void MD_basic::runNVE(int step1, double potential, ModuleBase::Vector3<double> *force, const ModuleBase::matrix &stress)
{
//-------------------------------------------------------------------------------
// Adiabatic ensemble 
// Molecular dynamics calculation with Verlet algorithm
//-------------------------------------------------------------------------------

    ModuleBase::TITLE("MD_basic","runNVE");
    ModuleBase::timer::tick("MD_basic","runNVE");
    step_=step1+step_rst_;

  // Set up extra output to ion optimizer / MD header
  //std::cout<<"Time interval       : "<< mdp.dt*fundamentalTime*1e15<< " (fs)"<<std::endl;
  //std::cout<<"Target temperature  : "<< temperature_/ModuleBase::K_BOLTZMAN_AU<< " (K)"<<std::endl;

    if(step_==1){
        for(int k=0;k<ucell.nat;k++)
        {
            if(ionmbl[k].x==0)vel[k].x=0;
            if(ionmbl[k].y==0)vel[k].y=0;
            if(ionmbl[k].z==0)vel[k].z=0;
        }
	if (ucell.set_vel==false)   // Yuanbo Li 2021/8/20
	{
	MD_func::scalevel(ucell.nat, frozen_freedom_, temperature_, vel, allmass);
	}
    }
    double twiceKE=MD_func::GetAtomKE(ucell.nat, vel, allmass);
    twiceKE = twiceKE * 2;

    double tempNow = twiceKE/(double(3*ucell.nat-frozen_freedom_))/ModuleBase::K_BOLTZMAN_AU;
    //std::cout<<" start temperature = "<< tempNow/ModuleBase::K_BOLTZMAN_AU<< " (k)"<<std::endl;

    if(GlobalV::STRESS)
	{
		outStressMD(stress, twiceKE);
	}
  
    for(int k=0;k<ucell.nat;k++)
    {
        if(ionmbl[k].x==0)force[k].x=0;
        if(ionmbl[k].y==0)force[k].y=0;
        if(ionmbl[k].z==0)force[k].z=0;
    }

    //std::cout<<"begin maxForce"<<std::endl;
    double maxForce = MD_func::MAXVALF(ucell.nat, force);
    //std::cout<<"maxForce: "<<sqrt(maxForce)<<std::endl; 

    energy_=potential;

    //double conservedE = MD_func::Conserved(twiceKE/2, energy_, 3*ucell.nat-frozen_freedom_);

    // Output the message to the screen.
    // std::cout << " ------------------------------------------------------------" << std::endl;
    // std::cout<<" "<<std::left<<std::setw(12)<<"NVE_STEP"<<std::left<<std::setw(12)<< "SystemE"<<std::left<<std::setw(12)<< "Conserved"<<std::left<<std::setw(12)<< "DeltaE"<<std::left<<std::setw(12)<< "Temperature"<<std::endl;
	// std::cout<<" "<<std::left<<std::setw(12)<<step_<<std::left<<std::setw(12)<< energy_<<std::left<<std::setw(12)<< conservedE <<std::left<<std::setw(12)<< energy_-oldEtot_<<std::left<<std::setw(12)<<twiceKE/(double(3*ucell.nat-frozen_freedom_))/ModuleBase::K_BOLTZMAN_AU<<std::endl;
	// std::cout << " ------------------------------------------------------------" << std::endl;
  
    // (1) 1st step of Verlet-Velocity
    // New velocity are obtained
	this->update_half_velocity(force);
    // (2) 2nd step of Verlet-Velocity 
    // Update the Non-Wrapped cartesion coordinate
    twiceKE = MD_func::GetAtomKE(ucell.nat, vel, allmass);
    twiceKE = 2 * twiceKE;
    if(step_!=1||mdp.rstMD==1)this->update_half_velocity(force);
    update_half_direct(1);
    update_half_direct(0);
    // Calculate the maximal velocities.
    // The step in fractional coordinates
    double maxStep = 0;
    for(int  i = 0;i< ucell.nat;i++)
    {
        if((pow(vel[i].x,2.0)+pow(vel[i].y,2.0)+pow(vel[i].z,2.0))>maxStep)
        {
            maxStep = pow(vel[i].x,2.0)+pow(vel[i].y,2.0)+pow(vel[i].z,2.0);
        }
    }
    getTaudUpdate();
    save_output_position();
    maxStep = sqrt(maxStep)*mdp.dt;


    // calculate the conserved quantity during MD 
    double hamiltonian = MD_func::Conserved(twiceKE/2, energy_, 3*ucell.nat-frozen_freedom_);

    //std::cout<< std::setprecision (9)<<hamiltonian<<" "<< std::setprecision (9)<<twiceKE/2<<std::endl;

    
    // Output the message to the screen.
    if (!GlobalV::MY_RANK)
    { 
        GlobalV::ofs_running << " ------------------------------------------------------------" << std::endl;
        GlobalV::ofs_running << " "<<std::left<<std::setw(12)<<"NVE_STEP"<<std::left<<std::setw(12)<< "SystemE"<<std::left<<std::setw(12)<< "Conserved"<<std::left<<std::setw(12)<< "DeltaE"<<std::left<<std::setw(12)<< "Temperature"<<std::endl;
	    GlobalV::ofs_running << " "<<std::left<<std::setw(12)<<step_<<std::left<<std::setw(12)<< energy_<<std::left<<std::setw(12)<< hamiltonian <<std::left<<std::setw(12)<< energy_-oldEtot_<<std::left<<std::setw(12)<<twiceKE/double((3*ucell.nat-frozen_freedom_))/ModuleBase::K_BOLTZMAN_AU<<std::endl;
	    GlobalV::ofs_running << " ------------------------------------------------------------" << std::endl;
    }
    std::cout << " ------------------------------------------------------------" << std::endl;
    std::cout << " "<<std::left<<std::setw(12)<<"NVE_STEP"<<std::left<<std::setw(12)<< "SystemE"<<std::left<<std::setw(12)<< "Conserved"<<std::left<<std::setw(12)<< "DeltaE"<<std::left<<std::setw(12)<< "Temperature"<<std::endl;
	std::cout << " "<<std::left<<std::setw(12)<<step_<<std::left<<std::setw(12)<< energy_<<std::left<<std::setw(12)<< hamiltonian <<std::left<<std::setw(12)<< energy_-oldEtot_<<std::left<<std::setw(12)<<twiceKE/(double(3*ucell.nat-frozen_freedom_))/ModuleBase::K_BOLTZMAN_AU<<std::endl;
	std::cout << " ------------------------------------------------------------" << std::endl;

    oldEtot_=energy_;
    MD_func::mdRestartOut(step_, mdp.recordFreq, ucell.nat, vel);

    //std::cout<<"(NVE): this step finished."<<std::endl;
    // if (!GlobalV::MY_RANK)
    // {
    //     GlobalV::ofs_running<<"(NVE): this step finished."<<std::endl;
    // }
//2015-09-25, xiaohui
#ifdef __MPI
    if(step_%100==0)
    {
        MPI_Bcast(vel,ucell.nat*3,MPI_DOUBLE,0,MPI_COMM_WORLD);
 //       MPI_Bcast(&xLogS,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
 //       MPI_Bcast(&vLogS,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
#endif

    ModuleBase::timer::tick("MD_basic","runNVE");
    return;
}

bool MD_basic::runFIRE(int step1, double potential, ModuleBase::Vector3<double> *force, const ModuleBase::matrix &stress)
{
//-------------------------------------------------------------------------------
// REFERENCES:
//   
    ModuleBase::TITLE("MD_basic","runFIRE");
    ModuleBase::timer::tick("MD_basic","runFIRE");
    step_ = step1;
    // if(step_==1)
    // {
    //     std::cout<<" FIRE Start."<<std::endl;
    // }
    // else 
    // {
    //     std::cout<<" FIRE step: "<<step_<<std::endl;
    // }
  
    //std::cout<<" FIRE Fronzen Atom Number : "<<double(frozen_freedom_)/3.0<<std::endl;
    

    // Set up extra output to ion optimizer / MD header
    //std::cout<<"Time interval       : "<< mdp.dt*fundamentalTime*1e15<< " (fs)"<<std::endl;

    if(step_==1)
    {
        for(int k=0;k<ucell.nat;k++)
        {
            if(ionmbl[k].x==0)vel[k].x=0;
            if(ionmbl[k].y==0)vel[k].y=0;
            if(ionmbl[k].z==0)vel[k].z=0;
        }
    //      if(frozen_freedom_==0)RemoveMovementOfCenterOfMass();
        MD_func::scalevel(ucell.nat, frozen_freedom_, temperature_, vel, allmass);
    }
    double twiceKE=MD_func::GetAtomKE(ucell.nat, vel, allmass);
    twiceKE = twiceKE * 2;

    double tempNow = twiceKE/(double(3*ucell.nat-frozen_freedom_))/ModuleBase::K_BOLTZMAN_AU;
    //std::cout<<" start temperature = "<< tempNow/ModuleBase::K_BOLTZMAN_AU<< " (k)"<<std::endl;

    for(int k=0;k<ucell.nat;k++)
    {
        if(ionmbl[k].x==0)force[k].x=0;
        if(ionmbl[k].y==0)force[k].y=0;
        if(ionmbl[k].z==0)force[k].z=0;
    }

    //std::cout<<"begin maxForce"<<std::endl;
    double maxForce = MD_func::MAXVALF(ucell.nat, force);
    //std::cout<<"maxForce: "<<sqrt(maxForce)<<std::endl; 
    
    energy_=potential;

    // double conservedE = MD_func::Conserved(twiceKE/2, energy_, 3*ucell.nat-frozen_freedom_);

    // std::cout << " ------------------------------------------------------------" << std::endl;
    // std::cout << " "<<std::left<<std::setw(12)<<"FIRE_STEP"<<std::left<<std::setw(12)<< "SystemE"<<std::left<<std::setw(12)<< "Conserved"<<std::left<<std::setw(12)<< "DeltaE"<<std::left<<std::setw(12)<< "Temperature"<<std::endl;
	// std::cout << " "<<std::left<<std::setw(12)<<step_<<std::left<<std::setw(12)<< energy_<<std::left<<std::setw(12)<< conservedE <<std::left<<std::setw(12)<< energy_-oldEtot_<<std::left<<std::setw(12)<<twiceKE/(double(3*ucell.nat-frozen_freedom_))/ModuleBase::K_BOLTZMAN_AU<<std::endl;
	// std::cout << " FIRE Fronzen Atom Number : "<<double(frozen_freedom_)/3.0<<std::endl;
    // std::cout << " ------------------------------------------------------------" << std::endl;

    // std::cout<<"------------------------------------------------------------------------------"<<std::endl;
    // std::cout<< "FIRE STEP "<< step_<<std::endl;
    // std::cout<<"------------------------------------------------------------------------------"<<std::endl;
    // if (!GlobalV::MY_RANK){
    //      GlobalV::ofs_running<<"------------------------------------------------------------------------------"<<std::endl;
    //      GlobalV::ofs_running<< "FIRE STEP "<< step_<<std::endl;
    //      GlobalV::ofs_running<<"------------------------------------------------------------------------------"<<std::endl;
    // }
  
    // (1) 1st step of Verlet-Velocity
    // New velocity are obtained
	this->update_half_velocity(force);

    // (2) 2nd step of Verlet-Velocity 
    // Update the Non-Wrapped cartesion coordinate
    twiceKE = MD_func::GetAtomKE(ucell.nat, vel, allmass);
    twiceKE = 2 * twiceKE;
    if(step_!=1)this->update_half_velocity(force);

	double largest_grad_FIRE = 0.0;
	for(int i=0;i<ucell.nat;i++)
    {
        if(largest_grad_FIRE < abs(force[i].x))
        {
            largest_grad_FIRE= abs(force[i].x);
        } 
        if(largest_grad_FIRE < abs(force[i].y))
        {
            largest_grad_FIRE= abs(force[i].y);
        } 
        if(largest_grad_FIRE < abs(force[i].z))
        {
            largest_grad_FIRE= abs(force[i].z);
        } 
	}
	
	GlobalV::ofs_running << " LARGEST GRAD (eV/A)  : " << largest_grad_FIRE * ModuleBase::Ry_to_eV / 0.529177 << std::endl;

	if(largest_grad_FIRE*ModuleBase::Ry_to_eV/0.529177 < 0.01)
    {
	//"convergency reach"
		//std::cout <<"CONVERGENCY REACH of FIRE in the "<<step <<" steps " <<std::endl;
		ModuleBase::timer::tick("MD_basic","runFIRE");
		return true;
	}

	fire.check_FIRE(ucell.nat, force, mdp.dt, vel);
    
    update_half_direct(1);
    update_half_direct(0);
    // Calculate the maximal velocities.
    // The step in fractional coordinates
    double maxStep = 0;
    for(int  i = 0;i< ucell.nat;i++)
    {
        if((pow(vel[i].x,2.0)+pow(vel[i].y,2.0)+pow(vel[i].z,2.0))>maxStep)
        {
            maxStep = pow(vel[i].x,2.0)+pow(vel[i].y,2.0)+pow(vel[i].z,2.0);
        }
    }
    getTaudUpdate();

    save_output_position();
    maxStep = sqrt(maxStep)*mdp.dt;

    double hamiltonian = MD_func::Conserved(twiceKE/2, energy_, 3 * ucell.nat - frozen_freedom_);

    // Output the message to the screen.
    if (!GlobalV::MY_RANK)
    { 
        GlobalV::ofs_running << " ------------------------------------------------------------" << std::endl;
        GlobalV::ofs_running << " "<<std::left<<std::setw(12)<<"NVE_STEP"<<std::left<<std::setw(12)<< "SystemE"<<std::left<<std::setw(12)<< "Conserved"<<std::left<<std::setw(12)<< "DeltaE"<<std::left<<std::setw(12)<< "Temperature"<<std::endl;
	    GlobalV::ofs_running << " "<<std::left<<std::setw(12)<<step_<<std::left<<std::setw(12)<< energy_<<std::left<<std::setw(12)<< hamiltonian <<std::left<<std::setw(12)<< energy_-oldEtot_<<std::left<<std::setw(12)<<twiceKE/(double(3*ucell.nat-frozen_freedom_))/ModuleBase::K_BOLTZMAN_AU<<std::endl;
	    GlobalV::ofs_running << " ------------------------------------------------------------" << std::endl;
    }
    std::cout << " ------------------------------------------------------------" << std::endl;
    std::cout << " "<<std::left<<std::setw(12)<<"NVE_STEP"<<std::left<<std::setw(12)<< "SystemE"<<std::left<<std::setw(12)<< "Conserved"<<std::left<<std::setw(12)<< "DeltaE"<<std::left<<std::setw(12)<< "Temperature"<<std::endl;
	std::cout << " "<<std::left<<std::setw(12)<<step_<<std::left<<std::setw(12)<< energy_<<std::left<<std::setw(12)<< hamiltonian <<std::left<<std::setw(12)<< energy_-oldEtot_<<std::left<<std::setw(12)<<twiceKE/(double(3*ucell.nat-frozen_freedom_))/ModuleBase::K_BOLTZMAN_AU<<std::endl;
	std::cout << " FIRE Fronzen Atom Number : "<<double(frozen_freedom_)/3.0<<std::endl;
    std::cout << " ------------------------------------------------------------" << std::endl;

    // Output the message to the screen.
    oldEtot_=energy_;

    // std::cout<<"(NVE): this step finished."<<std::endl;
    // if (!GlobalV::MY_RANK){
    //     GlobalV::ofs_running<<"(NVE): this step finished."<<std::endl;
    // }
//2015-09-25, xiaohui
#ifdef __MPI
    if(step1%100==0){
        MPI_Bcast(vel,ucell.nat*3,MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
#endif

    ModuleBase::timer::tick("MD_basic","runFIRE");
    return false;
}

//update velocities of ions for half MD step
void MD_basic::update_half_velocity(ModuleBase::Vector3<double> *force)
{
    for(int  ii=0;ii<ucell.nat;ii++){ 
        vel[ii] = vel[ii] + force[ii]/allmass[ii]*mdp.dt/2.0;
    }
}

//update cartesian coordinate changes for half MD step, would turn to direct coordinate changes later
void MD_basic::update_half_direct(const bool is_restart)
{
    if(is_restart)
    {//first half
        for(int  ii=0;ii<ucell.nat;ii++)
        { 
            cart_change[ii] = vel[ii]*mdp.dt/2.0/ucell.lat0;
        }
    }
    else
    {//second half
        for(int  ii=0;ii<ucell.nat;ii++)
        { 
            cart_change[ii] += vel[ii]*mdp.dt/2.0/ucell.lat0;
        }
    }
}

//output Structure information for each MD step
void MD_basic::save_output_position()
{
    std::string posOutName("md_pos_");
	posOutName=MD_func::intTurnTostring(step_,posOutName);
	ucell.update_pos_taud(tauDirectChange);
	MD_func::printpos(posOutName, step_, mdp.recordFreq, ucell);
}

int MD_basic::getRealStep()
{
    return this->step_;
}

//output pressure of total MD system, P = tr(stress) + P_kin
void MD_basic::outStressMD(const ModuleBase::matrix& stress, const double& twiceKE)
{
    GlobalV::ofs_running<<"\noutput Pressure for check!"<<std::endl;
    double press = 0.0;
    for(int i=0;i<3;i++)
    {
        press += stress(i,i)/3;
    }
    double virial = press;
    press += twiceKE/3/ucell.omega; //output virtual press = 2/3 *Ek/V + sum(sigma[i][i])/3
    const double unit_transform = ModuleBase::HARTREE_SI / pow(ModuleBase::BOHR_RADIUS_SI,3) * 1.0e-8 ;
    GlobalV::ofs_running<<"Virtual Pressure is "<<press*unit_transform<<" Kbar "<<std::endl;
    GlobalV::ofs_running<<"Virial Term is "<<virial*unit_transform<<" Kbar "<<std::endl;
    GlobalV::ofs_running<<"Kenetic Term is "<<(press-virial)*unit_transform<<" Kbar "<<std::endl;
}

//turn cartesian coordinate changes to direct changes
void MD_basic::getTaudUpdate()
{
    ModuleBase::Vector3<double> fracStep;
	for(int  ii=0;ii<ucell.nat;ii++)
    { 
		ModuleBase::Mathzone::Cartesian_to_Direct(cart_change[ii].x,cart_change[ii].y,cart_change[ii].z,
					ucell.latvec.e11,ucell.latvec.e12,ucell.latvec.e13,
					ucell.latvec.e21,ucell.latvec.e22,ucell.latvec.e23,
					ucell.latvec.e31,ucell.latvec.e32,ucell.latvec.e33,
					fracStep.x,fracStep.y,fracStep.z);

		tauDirectChange[ii] = fracStep;
	}
}
