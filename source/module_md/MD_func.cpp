#include "MD_func.h"

bool MD_func::RestartMD(const int& numIon, ModuleBase::Vector3<double>* vel, int& step_rst)
{
	int error(0);
	double *vell=new double[numIon*3];
	if (!GlobalV::MY_RANK)
	{
		std::stringstream ssc;
		ssc << GlobalV::global_readin_dir << "Restart_md.dat";
		std::ifstream file(ssc.str().c_str());

		if(!file)
		{
			std::cout<<"please ensure whether 'Restart_md.dat' exists!"<<std::endl;
			error = 1;
		}
		if (!error)
		{
			file.ignore(11, '\n');

        //----------------------------------------------------------
        // main parameters
        //----------------------------------------------------------
			file.ignore(13, '\n');//read and assert number of atoms
			int num;
			file>>num;
			if(num != numIon)
			{
				std::cout<<"please ensure whether 'Restart_md.dat' right!"<<std::endl;
				error = 1;
			}
		}
		if (!error)
		{
			file.get();
			file.ignore(23, '\n');//read velocities
			for(int i = 0;i<numIon*3;i++)
			{
				file>>vell[i];
			}
			file.get();
			file.ignore(6, '\n');//read start step of MD
			file>>step_rst;
			file.close();
		}
	}
#ifdef __MPI
	MPI_Bcast(&error,1,MPI_INT,0,MPI_COMM_WORLD);
#endif
	if(error)
	{
		delete[] vell;
		exit(0);
	}
#ifdef __MPI
	MPI_Bcast(&step_rst,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(vell,numIon*3,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif
	for(int i=0;i<numIon;i++)
	{
		vel[i].x=vell[i*3];
		vel[i].y=vell[i*3+1];
		vel[i].z=vell[i*3+2];
	}

	delete []vell;
	return true;
}

void MD_func::mdRestartOut(const int& step, const int& recordFreq, const int& numIon, ModuleBase::Vector3<double>* vel)
{
//this function used for outputting the information of restart file
	bool pass;
	pass = 0;

	if (recordFreq==1||step==1||( recordFreq > 1&& step%recordFreq==0 ) )
		pass =1;
	if (!pass) return;

	if(!GlobalV::MY_RANK){
		std::stringstream ssc;
		ssc << GlobalV::global_out_dir << "Restart_md.dat";
		std::ofstream file(ssc.str().c_str());
		file<<"MD_RESTART"<<std::endl;
		file<<"ATOM_NUMBERS: "<<numIon<<std::endl;
		file<<"ION_VELOCITIES_(a.u.): "<<std::endl;
		for(int i=0;i<numIon;i++){
			file<<std::setprecision (12)<<vel[i].x<<" "<<std::setprecision (12)<<vel[i].y<<" "<<std::setprecision (12)<<vel[i].z<<std::endl;
		}
		file<<"step: "<<step<<std::endl;                
		file.close();
	}
	return;
}

double MD_func::GetAtomKE(const int& numIon, const ModuleBase::Vector3<double>* vel, const double * allmass){
//---------------------------------------------------------------------------
// DESCRIPTION:
//   This function calculates the classical kinetic energy of a group of atoms.
//----------------------------------------------------------------------------

	double ke = 0.0;

	// kinetic energy = \sum_{i} 1/2*m_{i}*(vx^2+vy^2+vz^2)
	for(int ion=0;ion<numIon;ion++){
		ke +=   0.5 * allmass[ion] * (vel[ion].x*vel[ion].x+vel[ion].y*vel[ion].y+vel[ion].z*vel[ion].z);
	}
	//std::cout<<"in GetAtomKE KE="<< ke<<std::endl;
	return ke;
}

// Read Velocity from STRU liuyu 2021-09-24
void MD_func::ReadVel(
	const UnitCell_pseudo &unit_in, 
	ModuleBase::Vector3<double>* vel)
{
	int iat=0;
    for(int it=0; it<unit_in.ntype; ++it)
    {
        for(int ia=0; ia<unit_in.atoms[it].na; ++ia)
        {
            vel[iat] = unit_in.atoms[it].vel[ia];
			if(unit_in.atoms[it].mbl[ia].x==0) vel[iat].x = 0;
			if(unit_in.atoms[it].mbl[ia].y==0) vel[iat].y = 0;
			if(unit_in.atoms[it].mbl[ia].z==0) vel[iat].z = 0;
            ++iat;
        }
    }
    assert(iat==unit_in.nat);
}

// Initial velocity randomly
void MD_func::RandomVel(
	const int& numIon, 
	const double& temperature, 
	const double* allmass,
	const int& frozen_freedom,
	const ModuleBase::Vector3<int>* ionmbl,
	ModuleBase::Vector3<double>* vel)
{
	if(!GlobalV::MY_RANK)
	{
		srand(time(0));
		ModuleBase::Vector3<double> average;
		ModuleBase::Vector3<double> mass;
		average.set(0,0,0);
		mass.set(0,0,0);
		for(int i=0; i<numIon; i++)
		{
			if(ionmbl[i].x==0)
			{
				vel[i].x = 0;
			}
			else
			{
				vel[i].x = rand()/double(RAND_MAX)-0.5;
				mass.x += allmass[i];
			}
			if(ionmbl[i].y==0)
			{
				vel[i].y = 0;
			}
			else
			{
				vel[i].y = rand()/double(RAND_MAX)-0.5;
				mass.y += allmass[i];
			}
			if(ionmbl[i].z==0)
			{
				vel[i].z = 0;
			}
			else
			{
				vel[i].z = rand()/double(RAND_MAX)-0.5;
				mass.z += allmass[i];
			}
			average += allmass[i]*vel[i];
		}

		average.x = average.x / mass.x;
		average.y = average.y / mass.y;
		average.z = average.z / mass.z;

		for(int i=0; i<numIon; i++)
    	{
			if(ionmbl[i].x) vel[i].x -= average.x;
			if(ionmbl[i].y) vel[i].y -= average.y;
			if(ionmbl[i].z) vel[i].z -= average.z;
		}
	
		double factor = 0.5*(3*numIon-frozen_freedom)*temperature/GetAtomKE(numIon, vel, allmass);
		for(int i=0; i<numIon; i++)
    	{
        	vel[i] = vel[i]*sqrt(factor);
    	}
	}

#ifdef __MPI
	MPI_Bcast(vel,numIon*3,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif
	return;
}

void MD_func::InitVel(
	const UnitCell_pseudo &unit_in, 
	const double& temperature, 
	double* allmass,
	int& frozen_freedom,
	ModuleBase::Vector3<int>* ionmbl,
	ModuleBase::Vector3<double>* vel)
{
	//frozen_freedom = getMassMbl(unit_in, allmass, ionmbl);

	if(unit_in.set_vel)
    {
        ReadVel(unit_in, vel);
    }
    else
    {
        RandomVel(unit_in.nat, temperature, allmass, frozen_freedom, ionmbl, vel);
    }
}

void MD_func::InitVelocity(
	const int& numIon, 
	const double& temperature, 
	const double& fundamentalTime, 
	const double* allmass,
	ModuleBase::Vector3<double>* vel)
{
	if(!GlobalV::MY_RANK){ //xiaohui add 2015-09-25
		srand(time(0));
		int iy=rand()%21,im=rand()%11,id=rand()%21,ih=rand()%31,in=rand()%61,is=rand()%41;
		int jy=rand()%21,jm=rand()%11,jd=rand()%21,jh=rand()%31,jn=rand()%61,js=rand()%41;
		//    int iy=13,im=10,id=1,ih=21,in=54,is=23;
		//  int jy=14,jm=6,jd=19,jh=4,jn=29,js=32;
		long int i,j;
		long int i0,i1,i2,j0,j1,j2;
		long int a=16807,m=2147483647,q=127773,r=2836;
		long int M;  
		double x,y,u,v;
		i=iy+70*(im+12*(id+31*(ih+23*(in+59*is))));				//initial the seeds
		j=jy+70*(jm+12*(jd+31*(jh+23*(jn+59*js))));	
		i0=a*(i%q)-r*((i-i%q)/q);
		j0=a*(j%q)-r*((j-j%q)/q);
		if(i0>=0)	i1=i0;  	else i1=i0+m;	
		i0=a*(j%q)-r*((j-j%q)/q);
		if(j0>=0)	j1=j0;  	else j1=j0+m;	
		for(M=0;M<3*numIon;){	
			i0=a*(i1%q)-r*((i1-i1%q)/q);                        
			if(i0>=0)	i2=i0;  
			else i2=i0+m;									
			u=((double)i1)/((double)m);		       //first ramdon number        
			i1=i2;
			j0=a*(j1%q)-r*((j1-j1%q)/q);
			if(j0>=0)	j2=j0;  
			else j2=j0+m;										
			v=((double)j1)/((double)m);		        //second ramdon number
			j1=j2;
			x=tan(ModuleBase::PI*(u-0.5));
			y=1.6*v/(ModuleBase::PI*(x*x+1.0));
			if(y<=(1.0/sqrt(2.0*ModuleBase::PI)*exp(-0.5*x*x))){
				if(M<numIon) vel[M].x=x;
				else if(M<2*numIon) vel[M-numIon].y=x;
				else vel[M-2*numIon].z=x;
					// std::cout<<"this vel is:"<<x<<std::endl;
					M++;
			}
		}
      
		int ion = 0;
		for(ion=0;ion<numIon;ion++){
			vel[ion]/=sqrt(allmass[ion]*9.01938e-31);				 
		}
		for(ion=0;ion<numIon;ion++){
			vel[ion]*=sqrt(temperature/ModuleBase::K_BOLTZMAN_AU * ModuleBase::K_BOLTZMAN_SI);
			vel[ion]*= fundamentalTime/ModuleBase::BOHR_RADIUS_SI;
		//rescale the velocities to a.u.
		}
	} //xiaohui add 2 lines 2015-09-25, bcast vel
//2015-09-25, xiaohui
#ifdef __MPI
	MPI_Bcast(vel,numIon*3,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif
	return;
}

/*void MD_func::ReadNewTemp(int step )
{
//-----------------------------------------------------------------
//  If fixTemperature == 0, then this subroutine will be skipped
//  otherwise, we are going to read a new temperature from the disk.
//  You only need to create a file named ChangeTemp.dat, and then 
//  change the temperature in it. 
//-----------------------------------------------------------------

	double intemp;
	int sstep=0;

	if ( fixTemperature > 0 )
	{

	// Change the temperature every 'fixTemperature' steps.
		if( (step!=1)&&(step%fixTemperature == 1) )
		{
			
			// Read in new temperature from file.
			std::ifstream file;
			file.open("ChangeTemp.dat");
			if (!file){
				std::cout<<"ERROR IN OPENING ChangeTemp.dat, CODE STOP!"<<std::endl;
				exit(0);
			}
			while((sstep+fixTemperature)<step){
				file>> intemp;
				sstep+=fixTemperature;
			}
			file.close();

			// Renew information.
			intemp =  intemp * ModuleBase::K_BOLTZMAN_AU;
			if ( fabs(intemp-temperature) >1e-6 ) {
				std::cout <<"(ReadNewTemp): Read in new temp:"<< intemp/ModuleBase::K_BOLTZMAN_AU 
					<<" previous temp:"<< temperature/ModuleBase::K_BOLTZMAN_AU<<std::endl;
				temperature = intemp;
			}
			else{
				std::cout<<"(ReadNewTemp): new temp:"<< intemp/ModuleBase::K_BOLTZMAN_AU
					<<" previous temp:"<<temperature/ModuleBase::K_BOLTZMAN_AU
					<< ". No change of temp."<<std::endl;
			}
		}
	}

	return;
}*/

//int to std::string and add to output path
std::string MD_func::intTurnTostring(long int iter, std::string path)
{
	long int i[10],k=0;
	if(iter>9999999999) return "error!";
	for(int j=9; j>-1; j--)
	{
		if(iter==0) continue;
		if(iter>pow(10,j)-1)
		{
			i[k] = iter%10;
			iter /= 10;
			k++;
		}
	}
	for(int j=k-1;j>-1;j--){
		path+=(i[j]+48);
	}
	return path;
}

int MD_func::getMassMbl(const UnitCell_pseudo &unit_in, 
			const MD_parameters &mdp,
			double* allmass, 
			ModuleBase::Vector3<int>* ionmbl)
{
//some prepared information
//mass and degree of freedom
	int ion=0;
	int frozen_freedom=0;
	for(int it=0;it<unit_in.ntype;it++){
		for(int i=0;i<unit_in.atoms[it].na;i++){
			allmass[ion]=unit_in.atoms[it].mass/ModuleBase::AU_to_MASS;
			ionmbl[ion]=unit_in.atoms[it].mbl[i];
			if (ionmbl[ion].x==0) frozen_freedom++;
			if (ionmbl[ion].y==0) frozen_freedom++;
			if (ionmbl[ion].z==0) frozen_freedom++;

			ion++;
		}
	}

	// the center of mass is fixed except for NVT Anderson
	if(!(mdp.mdtype==1 && mdp.NVT_control==3))
	{
		frozen_freedom += 3;
	}

	return frozen_freedom;
}

void MD_func::printpos(const std::string& file, const int& iter, const int& recordFreq, const UnitCell_pseudo& unit_in)
{
//intend to output the positions of atoms to ordered file
	bool pass;
	pass = 0;

	if (recordFreq==1||iter==1||( recordFreq > 1&& iter%recordFreq==0 ) )
		pass =1;
	if (!pass) return;

	std::string file1=file+".xyz";
	std::string file2=file+".cif";

	//xiaohui add 'GlobalV::OUT_LEVEL', 2015-09-16
	if(GlobalV::OUT_LEVEL == "i"||GlobalV::OUT_LEVEL == "ie") unit_in.print_tau();
	if(GlobalV::OUT_LEVEL == "i"||GlobalV::OUT_LEVEL == "ie") unit_in.print_cell_xyz(file1);
	unit_in.print_cell_cif(file2);
	std::stringstream ss;

	ss << GlobalV::global_out_dir << "STRU_MD";

	//zhengdy modify 2015-05-06, outputfile "STRU_Restart"
#ifdef __LCAO
	unit_in.print_stru_file(GlobalC::ORB, ss.str(),2);
#else
	unit_in.print_stru_file(ss.str(),2);
#endif

	return;
}

//rescale velocities to target temperature.
void MD_func::scalevel(
	const int& numIon,
	const int& nfrozen,
	const double& temperature,
	ModuleBase::Vector3<double>* vel,
	const double* allmass
)
{
	double ke=GetAtomKE(numIon, vel, allmass);
	if(ke>1e-9)
	{
		for(int i=0;i<numIon;i++)
		{
			vel[i]*=sqrt((3*numIon-nfrozen)*temperature/ke/2);
		}
	}
	return;
}

double MD_func::Conserved(const double KE, const double PE, const int nfreedom){
//---------------------------------------------------------------------------
//   This function calculates the conserved quantity for the NVE system. 
//----------------------------------------------------------------------------

   	// KE   // Kinetic energy of particles (due to their velocity)
   	// PE   // Potential energy (DFT total energy)
	// number //number of atoms which have full freedoms

  	double Conserved; // The conserved quantity

   	Conserved = KE + PE ;
   
	if (!GlobalV::MY_RANK)
	{               
		GlobalV::ofs_running<< "--------------------------------------------------"<<std::endl;
        GlobalV::ofs_running<< "            SUMMARY OF NVE CALCULATION            "<<std::endl;
        GlobalV::ofs_running<<" --------------------------------------------------"<<std::endl;  
		GlobalV::ofs_running<< "NVE Conservation     : "<< Conserved<<" (Hartree)"<<std::endl;
		GlobalV::ofs_running<< "NVE Temperature      : "<< 2*KE/(nfreedom)*ModuleBase::Hartree_to_K<<" (K)"<<std::endl;
		GlobalV::ofs_running<< "NVE Kinetic energy   : "<< KE<<" (Hartree)"<<std::endl;
		GlobalV::ofs_running<< "NVE Potential energy : "<< PE<<" (Hartree)"<<std::endl;
	}
   	return Conserved;
}

double MD_func::MAXVALF(const int numIon, const ModuleBase::Vector3<double>* force){
	//std::cout<<"enter in MAXVALF"<<std::endl;
	double max=0;
	for(int i=0;i<numIon;i++){
		double force0 = pow(force[i].x,2)+pow(force[i].y,2)+pow(force[i].z,2);
		if(max<force0) max = force0;
	}
	return max;
}
