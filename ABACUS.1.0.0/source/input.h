//==========================================================
// Author: Lixin He,mohan
// DATE : 2008-12-24
//==========================================================
#ifndef INPUT_H
#define INPUT_H

#include <string>
using namespace std;

//==========================================================
// CLASS :
// NAME : Input(readin "input_parameters" from "in_data" dir)
//==========================================================
class Input
{
public:
    
	Input();
    ~Input();
    void Init(const string &fn);
    void Print(const string &fn)const;
	void close_log(void)const;

	int npool;

//==========================================================
// MEMBER VARIABLES :
// main parameters
//==========================================================
    string suffix;			// suffix of out put dir
    string atom_file;		//xiaohui modify 2015-02-01
    string pseudo_dir;      // dir of pseudopotential
    
    string pseudo_type;     // the type of pseudopotential, mohan add 2013-05-20, we can accept
			    // UPF format (default) and vwr format. (xiaohui add 2013-06-23)
    string kpoint_file;		//xiaohui modify 2015-02-01
	string wannier_card;	// input card for wannier functions.
    string epm_pseudo_card; //
    string latname;			// lattice name
    string calculation;		// "scf" : self consistent calculation.
    // "nscf" : non-self consistent calculation.
    int ntype;				// number of atom types
    int nbands;				// number of bands
	int nbands_istate;		// number of bands around fermi level for istate calculation.
    int epm_spin_orbital;	// 0: no so effect; 1: plane wave; 2: small box
	bool epm_zeeman;
	double epm_mag_field_x;
	double epm_mag_field_y;
	double epm_mag_field_z;

    bool berry_phase;		// berry phase calculation
    int efield;				// add electrical field
	int edir;
	double emaxpos;
	double eopreg;
	double eamp;
	int bfield;				// add magnetic field
	double bfield_teslax;	//
	double bfield_teslay;	//
	double bfield_teslaz;	//
	double bfield_gauge_x;  // add by sunzhuyuan,
	double bfield_gauge_y;
	double bfield_gauge_z;
	bool opt_epsilon2;		// true : calculate the dielectric functions
	int  opt_nbands;		// number of bands for optical transition matrix
    bool lda_plus_u;		// true : lda plus u calculation
//==========================================================
// electrons / spin
//==========================================================
    string dft_functional;	// input DFT functional.
	int nspin;			// LDA ; LSDA ; non-linear spin
    double nelec;			// total number of electrons
    int lmaxmax;
    double tot_magnetization;
//==========================================================
// new function
//==========================================================
    //int local_basis; xiaohui modify 2013-09-01, for structural adjustment 	// 0:pw 1:lo 2:atomic
    //int linear_scaling; xiaohui modify 2013-09-01	//
	string basis_type; //xiaohui add 2013-09-01, for structural adjustment
	string ks_solver; //xiaohui add 2013-09-01
    bool sparse_matrix;
	int atom_distribution;
    bool symmetry;				//
	bool mlwf_flag;
	int vna;				//control the vna term.
	int grid_speed;			//1:normal 2:fast, mohan add 2012-03-29
    int force;
    bool force_set;
    double force_thr;		// threshold of force in unit (Ry/Bohr)
	double force_thr_ev2;	// invalid force threshold, mohan add 2011-04-17
    double stress_thr;      // Pengfei Li 2017-11-01 //LiuXh update 20180515
	bool stress;			// calculate the stress
	string ion_dynamics;		// methods to move_ion: sd, bfgs, cg...
    double cg_threshold;    // threshold when cg to bfgs, pengfei add 2011-08-15
	double bfgs_w1;			// wolfe condition 1	
	double bfgs_w2;			// wolfe condition 2
	double trust_radius_max;	// trust radius max
	double trust_radius_min;	// trust radius min
	double trust_radius_ini;	// initial move
//	bool gauss_pao_flag
//==========================================================
// planewave
//==========================================================
    //bool gamma_only;          
    bool gamma_only;		// for plane wave.
	bool gamma_only_local;	// for local orbitals.
    double ecutwfc;			// energy cutoff for wavefunctions
    double ecutrho;			// energy cutoff for charge/potential
    int ncx,ncy,ncz;		// three dimension of FFT charge/grid
    int nx,ny,nz;			// three dimension of FFT wavefunc
	int bx,by,bz;			// big mesh ball. mohan add 2011-04-21 
//==========================================================
// technique 
//==========================================================
    //string diago_type; xiaohui modify 2013-09-01		// "cg","davidson","fs"
	int diago_proc;			// the number of procs used to diag. mohan add 2012-01-13
    int diago_cg_maxiter;
	int diago_cg_prec;		// mohan add 2012-03-31
    int diago_david_ndim;
    double ethr;			// used in cg method
    double fs_ref_energy;	// folded-spectrum reference energy
	int nb2d;				// matrix 2d division.
	int nurse;				// used for debug.
	bool colour;			// used for fun.
	int t_in_h;			// calculate the T or not.
	int vl_in_h;			// calculate the T or not.
	int vnl_in_h;			// calculate the T or not.
	int zeeman_in_h;		// calculate the zeeman energy or not.
	int test_force;			// test the force.
	int test_stress;		// test the stress.
//==========================================================
// iteration
//==========================================================
    double dr2;				// \sum |rhog_out - rhog_in |^2
    int niter;				// number of max elec iter
    int nstep;				// number of max ionic iter
	int out_stru;			// outut stru file each ion step
	string out_level;		// control the output information.
//==========================================================
// occupation
//==========================================================
    string occupations;		// "fixed","smearing","tetrahedra","from_input"
    string smearing;		// "gaussian",
    // "mp","methfessel-paxton"
    // "mv","marzari-vanderbilt","cold"
    // "fd","fermi-dirac"
    double degauss;			//
//==========================================================
// charge mixing
//==========================================================
    string mixing_mode;		// "plain","broden",...
    double mixing_beta;		// 0 : no_mixing
    int mixing_ndim;		// used in Broden method
	double mixing_gg0;      // used in kerker method. mohan add 2014-09-27
//==========================================================
// potential / charge / wavefunction / energy
//==========================================================
    string restart_mode;	//
    string start_wfc;		// "file","atomic","random"
    string start_pot;		// "file","atomic"
	string charge_extrap;			// xiaohui modify 2015-02-01
	int mem_saver;			// 1: save psi when nscf calculation.
	int printe;				// mohan add 2011-03-16
    int out_charge;		// output charge density.
	int out_dm; // output density matrix.
	int out_potential;		// yes or no
    bool out_wf;			// yes or no
	int out_dos;			// dos calculation. mohan add 20090909
        int out_band;                   // band calculation pengfei 2014-10-13
	int out_hs;			// output H matrix and S matrix in local basis.
	bool out_lowf;			// output the wave functions in local basis.
	bool out_alllog; 		// output all logs.
	double dos_emin_ev;
	double dos_emax_ev;
	double dos_edelta_ev;
        double b_coef;           //  pengfei 2014-10-13
//==========================================================
// parameters about LCAO calculation. 
// mohan add 2009-11-11 
//==========================================================
	double lcao_ecut;		// ecut of two center integral
	double lcao_dk;			// delta k used in two center integral
	double lcao_dr;			// dr used in two center integral
	double lcao_rmax;		// rmax(a.u.) to make table.
    double search_radius;	//11.1
	bool search_pbc; // 11.2

	// Divide&Conqure
	int dc_nx;
	int dc_ny;
	int dc_nz;

	// selinv method parameter (cooperate with LinLin)
	int selinv_npole;
	double selinv_temp;
	double selinv_gap;
	double selinv_deltae;
	double selinv_mu;
	double selinv_threshold;
	int selinv_niter;

/*
	// md and related parameters.
	double md_dt; //dt: time step for MD.
	int    md_restart; // mohan add: restart from file.
	double md_tolv; // tolv: tolerence for velocity scaling.
	string md_thermostat; //ionic temperature: various md_thermostat  
	double md_temp0; // tempw: start temperature. 
	int md_tstep; // nraise: the temperature will reduce every md_tstep.
	double md_delt;	// delta_t: the reduce amount of temperature.
*/

//md and related parameters(added by zheng da ye)
        int md_mdtype;                   //choose ensemble
        double md_tauthermo;
        double md_taubaro;
        double md_dt;                    //time step
        int md_nresn;                     //parameter during integrater
        int md_nyosh;                      //parameter during integrater
        double md_qmass;                   //mass of thermostat
        double md_tfirst;                    //temperature begin
        double md_tlast;                    //temperature end
        int md_dumpmdfred;                  //The period to dump MD information for monitoring and restarting MD
        string md_mdoutpath;                //output path for md
	bool md_domsd;                   //whether compute <r(t)-r(0)>
        bool md_domsdatom;                //whether compute msd for each atom
        int md_rstmd;                    //whether restart;
        int md_outputstressperiod;      //period to output stress
        int md_fixtemperature;          //period to change temperature
        double md_ediff;             //parameter for constraining total energy change
	double md_ediffg;             //parameter for constraining max force change
        int md_msdstartTime;            //choose which step that msd be calculated 
//end of zhengdaye's add
//==========================================================
// vdwD2
// Peize Lin add 2014-03-31
//==========================================================
	bool vdwD2;					//calculate vdwD2 or not
	double vdwD2_scaling;
	double vdwD2_d;				//damping parameter
	string vdwD2_C6_file;
	string vdwD2_C6_unit;		//"Bohr" or "Angstrom"
	string vdwD2_R0_file;
	string vdwD2_R0_unit;		//"Bohr" or "Angstrom"
	string vdwD2_model;			//"period" or "radius"
	int vdwD2_period[3];
	double vdwD2_radius;
	string vdwD2_radius_unit;	//"Bohr" or "Angstrom"



//==========================================================
// spectrum
// pengfei Li add 2016-11-23
//==========================================================    
        //bool     epsilon;               // calculate epsilon or not
		string   spectral_type;          // the type of the calculated spectrum
		int      spectral_method;        // 0: tddft(linear response)
		int      eels_method;            // 0: hilbert_transform method; 1: standard method
		int      absorption_method;      // 0: vasp's method  1: pwscf's method
		//int		 epsilon_choice;         // 0: hilbert_transform method; 1: standard method
		
		string   kernel_type;           // the kernel type: rpa, tdlda ...
 
        string system;                 // bulk or surface
        double  eta;                   // unit(Ry)
        double  domega;                // unit(Ry)
        int     nomega;
        int     ecut_chi;                   // the dimension of G 
        //int     oband;                 // the number of "occupied" bands  
		double  q_start[3];            // the position of the first q point in direct coordinate
		double  q_direct[3];             // the q direction
        //int     start_q;               // the serial number of the start qpoint
        //int     interval_q;            // the interval of the qpoints
        int     nq;                    // the total number of qpoints for calculation
        bool     out_epsilon;           // output epsilon or not
        bool     out_chi;               // output chi or not
        bool     out_chi0;              // output chi0 or not
        double  fermi_level;            // the change the fermi level(Ry)
        bool     coulomb_cutoff;         // turn on or off the Coulomb_cutoff 0/1
		
		//bool     epsilon0;              // calculate the macroscopic dielectric constant or not
		//double   intersmear;            // eta
		double   intrasmear;            // Eta
		double   shift;
		bool     metalcalc;             // metal or not
		double   eps_degauss;            // degauss
		
		//int	epsilon0_choice;             // 0: vasp's method  1: pwscf's method

        bool     kmesh_interpolation;          // calculting <i,0|j,R>
        double  qcar[100][3];          // the Cartesian position of q points(unit: 2*PI/lat0) 
        int     lcao_box[3];           // the scale for searching the existence of the overlap <i,0|j,R>
	
	//added by zhengdy-soc
	bool noncolin;
	bool lspinorb;
	bool starting_spin_angle;
	double *angle1;
	double *angle2;
		


	//xiaohui add 2015-09-16
	bool input_error;
	
	private:
//==========================================================
// MEMBER FUNCTIONS :
// NAME : Read()
// NAME : Default(set the default value for variables)
// NAME : Check(check the values)
// NAME : Bcast(only used in paralle case,only read in data
//        from first cpu, and distribute the data to the
//        other processors)
//==========================================================
    bool Read(const string &fn);
    void Default(void);
    void Check(void);
#ifdef __MPI
    void Bcast(void);
#endif

	public:
//==========================================================
// MEMBER FUNCTION :
// NAME : read_value(a module used to read in datas)
//==========================================================
    template <class T>
    static void read_value(ifstream &ifs, T &var)
    {
        ifs >> var;
        ifs.ignore(150, '\n');
        return;
    }

	void strtolower(char *sa, char *sb);
	void readbool(ifstream &ifs, bool &var);
};

extern Input INPUT;
#endif //INPUT
