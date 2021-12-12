#include "MD_func.h"
#include "cmd_neighbor.h"
#include "LJ_potential.h"
#include "DP_potential.h"
#include "../input.h"
#include "../module_neighbor/sltk_atom_arrange.h"
#include "../module_neighbor/sltk_grid_driver.h"
#include "../module_base/global_variable.h"

#ifndef __CMD
#include "../src_pw/run_md_pw.h"
#endif

#ifdef __LCAO
#include "../src_lcao/run_md_lcao.h"
#endif


double MD_func::gaussrand()
{
    static double V1, V2, S;
    static int phase = 0;
    double X;
     
    if ( phase == 0 ) 
    {
        do 
        {
            double U1 = rand()/double(RAND_MAX);
            double U2 = rand()/double(RAND_MAX);
             
            V1 = 2.0 * U1 - 1.0;
            V2 = 2.0 * U2 - 1.0;
            S = V1 * V1 + V2 * V2;
        } while(S >= 1 || S == 0);
         
        X = V1 * sqrt(-2.0 * log(S) / S);
    } 
    else
    {
        X = V2 * sqrt(-2.0 * log(S) / S);
    }

    phase = 1 - phase;
 
    return X;
}

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

double MD_func::GetAtomKE(
		const int &numIon,
		const ModuleBase::Vector3<double> *vel, 
		const double *allmass)
{
	double ke = 0;

	for(int ion=0; ion<numIon; ++ion)
	{
		ke += 0.5 * allmass[ion] * vel[ion].norm2();
	}

	return ke;
}

void MD_func::kinetic_stress(
		const UnitCell_pseudo &unit_in,
		const ModuleBase::Vector3<double> *vel, 
		const double *allmass, 
		double &kinetic,
		ModuleBase::matrix &stress)
{
//---------------------------------------------------------------------------
// DESCRIPTION:
//   This function calculates the classical kinetic energy of atoms
//   and its contribution to stress.
//----------------------------------------------------------------------------

	kinetic = MD_func::GetAtomKE(unit_in.nat, vel, allmass);

	ModuleBase::matrix temp;
	temp.create(3,3);    // initialize

	for(int ion=0; ion<unit_in.nat; ++ion)
	{
		for(int i=0; i<3; ++i)
		{
			for(int j=i; j<3; ++j)
			{
				temp(i, j) += allmass[ion] * vel[ion][i] * vel[ion][j];
			}
		}
	}

	for(int i=0; i<3; ++i)
	{
		for(int j=0; j<3; ++j)
		{
			if(j<i) 
			{
				stress(i, j) = stress(j, i);
			}
			else
			{
				stress(i, j) = temp(i, j)/unit_in.omega;
			}
		}
	}
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
	const ModuleBase::Vector3<int> frozen,
	const ModuleBase::Vector3<int>* ionmbl,
	ModuleBase::Vector3<double>* vel)
{
	if(!GlobalV::MY_RANK)
	{
		ModuleBase::Vector3<double> average;
		ModuleBase::Vector3<double> mass;
		average.set(0,0,0);
		mass.set(0,0,0);
		for(int i=0; i<numIon; i++)
		{
			for(int k=0; k<3; ++k)
			{
				if(ionmbl[i][k]==0)
				{
					vel[i][k] = 0;
				}
				else
				{
					vel[i][k] = rand()/double(RAND_MAX)-0.5;
					mass[k] += allmass[i];
				}
			}
			average += allmass[i]*vel[i];
		}

		for(int i=0; i<numIon; i++)
    	{
			for(int k=0; k<3; ++k)
			{
				if(ionmbl[i][k] && frozen[k]==0)
				{
					vel[i][k] -= average[k] / mass[k];
				}
			}
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
	ModuleBase::Vector3<int> frozen;
	getMassMbl(unit_in, allmass, frozen, ionmbl);
	frozen_freedom = frozen.x + frozen.y + frozen.z;
	if(frozen.x == 0) ++frozen_freedom;
	if(frozen.y == 0) ++frozen_freedom;
	if(frozen.z == 0) ++frozen_freedom;

	if(unit_in.set_vel)
    {
        ReadVel(unit_in, vel);
    }
    else
    {
        RandomVel(unit_in.nat, temperature, allmass, frozen_freedom, frozen, ionmbl, vel);
    }
}

void MD_func::InitPos(
	const UnitCell_pseudo &unit_in, 
	ModuleBase::Vector3<double>* pos)
{
	int ion=0;
	for(int it=0;it<unit_in.ntype;it++)
	{
		for(int i=0;i<unit_in.atoms[it].na;i++)
		{
			pos[ion] = unit_in.atoms[it].tau[i]*unit_in.lat0;
			ion++;
		}
	}
}

//calculate potential, force and virial
void MD_func::force_virial(const int &istep,
		const MD_parameters &mdp,
		const UnitCell_pseudo &unit_in,
		double &potential,
		ModuleBase::Vector3<double> *force,
		ModuleBase::matrix &stress)
{
	ModuleBase::TITLE("MD_func", "force_stress");
    ModuleBase::timer::tick("MD_func", "force_stress");

	if(mdp.md_potential == "LJ")
	{
		bool which_method = unit_in.judge_big_cell();
		if(which_method)
		{
			CMD_neighbor cmd_neigh;
			cmd_neigh.neighbor(unit_in);

			potential = LJ_potential::Lennard_Jones(
								unit_in,
								cmd_neigh,
								force,
								stress);
		}
		else
		{
			Grid_Driver grid_neigh(GlobalV::test_deconstructor, GlobalV::test_grid_driver, GlobalV::test_grid);
			atom_arrange::search(
					GlobalV::SEARCH_PBC,
					GlobalV::ofs_running,
					grid_neigh,
					unit_in, 
					GlobalV::SEARCH_RADIUS, 
					GlobalV::test_atom_input,
					INPUT.test_just_neighbor);

			potential = LJ_potential::Lennard_Jones(
								unit_in,
								grid_neigh,
								force,
								stress);
		}
	}
	else if(mdp.md_potential == "DP")
	{
		DP_potential::DP_pot(unit_in, potential, force, stress);
	}
#ifndef __CMD
	else if(mdp.md_potential == "FP")
	{
		if(GlobalV::BASIS_TYPE=="pw" || GlobalV::BASIS_TYPE=="lcao_in_pw")
		{
			Run_MD_PW md_pw;
			md_pw.md_force_virial(istep, unit_in.nat, potential, force, stress);
		}
#ifdef __LCAO
		else if(GlobalV::BASIS_TYPE=="lcao")
		{
			//Run_lcao::lcao_line();
		}
#endif
	}
#endif
	else
	{
		ModuleBase::WARNING_QUIT("md_force_stress", "Unsupported MD potential !");
	}

	ModuleBase::timer::tick("MD_func", "md_force_stress");
}

void MD_func::outStress(const ModuleBase::matrix &virial, const ModuleBase::matrix &stress)
{
	GlobalV::ofs_running<<"\noutput Pressure for check!"<<std::endl;
    double stress_scalar = 0.0, virial_scalar = 0.0;
    for(int i=0;i<3;i++)
    {
        stress_scalar += stress(i,i)/3;
		virial_scalar += virial(i,i)/3;
    }
    const double unit_transform = ModuleBase::HARTREE_SI / pow(ModuleBase::BOHR_RADIUS_SI,3) * 1.0e-8;
    GlobalV::ofs_running<<"Virtual Pressure is "<<stress_scalar*unit_transform<<" Kbar "<<std::endl;
    GlobalV::ofs_running<<"Virial Term is "<<virial_scalar*unit_transform<<" Kbar "<<std::endl;
    GlobalV::ofs_running<<"Kenetic Term is "<<(stress_scalar-virial_scalar)*unit_transform<<" Kbar "<<std::endl;

	GlobalV::ofs_running << std::setprecision(6) << std::setiosflags(ios::showpos) << std::setiosflags(ios::fixed) << std::endl;
	ModuleBase::GlobalFunc::NEW_PART("MD STRESS (KBAR)");
	for (int i=0; i<3; i++)
	{
		GlobalV::ofs_running << " " << std::setw(15) << stress(i,0)*unit_transform 
			<< std::setw(15)<< stress(i,1)*unit_transform 
			<< std::setw(15) << stress(i,2)*unit_transform << std::endl;

	}
	GlobalV::ofs_running << std::setiosflags(ios::left);
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

void MD_func::MDdump(const int &step, 
		const int &natom,
		const ModuleBase::matrix &virial, 
		const ModuleBase::Vector3<double> *force)
{
	if(GlobalV::MY_RANK) return;

	std::stringstream file;
    file << GlobalV::global_out_dir << "MD_dump";
	std::ofstream ofs;
	if(step == 0)
	{
		ofs.open(file.str(), ios::trunc);
	}
	else
	{
		ofs.open(file.str(), ios::app);
	}

	const double unit_virial = ModuleBase::HARTREE_SI / pow(ModuleBase::BOHR_RADIUS_SI,3) * 1.0e-8;
	const double unit_force = ModuleBase::Hartree_to_eV*ModuleBase::ANGSTROM_AU;

	ofs << "MDstep:  " << step << std::endl;

	ofs << "VIRIAL (KBAR)" << std::endl;
	ofs << std::setprecision(12) << std::setiosflags(ios::fixed);
	for(int i=0; i<3; ++i)
	{
		ofs << std::setw(18) << virial(i, 0) * unit_virial 
			<< std::setw(18) << virial(i, 1) * unit_virial 
			<< std::setw(18) << virial(i, 2) * unit_virial << std::endl;
	}

	ofs << "\nFORCE (eV/Angstrom)" << std::endl;
	for(int i=0; i<natom; ++i)
	{
		ofs << std::setw(18) << force[i].x * unit_force 
			<< std::setw(18) << force[i].y * unit_force 
			<< std::setw(18) << force[i].z * unit_force << std::endl;
	}

	ofs << std::endl;
	ofs << std::endl;
	ofs.close();
}

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

void MD_func::getMassMbl(const UnitCell_pseudo &unit_in, 
			double* allmass, 
			ModuleBase::Vector3<int> &frozen,
			ModuleBase::Vector3<int>* ionmbl)
{
//some prepared information
//mass and degree of freedom
	int ion=0;
	frozen.set(0,0,0);
	for(int it=0;it<unit_in.ntype;it++){
		for(int i=0;i<unit_in.atoms[it].na;i++)
		{
			allmass[ion]=unit_in.atoms[it].mass/ModuleBase::AU_to_MASS;
			ionmbl[ion]=unit_in.atoms[it].mbl[i];
			if (ionmbl[ion].x==0) ++frozen.x;
			if (ionmbl[ion].y==0) ++frozen.y;
			if (ionmbl[ion].z==0) ++frozen.z;

			ion++;
		}
	}
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
	unit_in.print_stru_file(GlobalC::ORB, ss.str(), 2, 1);
#else
	unit_in.print_stru_file(ss.str(), 2, 1);
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
