#ifndef ENERGY_H
#define ENERGY_H
#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "module_elecstate/elecstate.h"
#include "src_lcao/local_orbital_wfc.h"
#include "module_psi/psi.h"
#include "module_elecstate/elecstate.h"

	class LCAO_Hamilt;

class energy
{
	public:
	friend class Electrons_Stochastic;//mohan add 2021-01-31 
	friend class LOOP_elec;
	friend class Local_Orbital_Charge;
	friend class Threshold_Elec;
	friend class Forces<double>;
	friend class Charge;
	friend class Potential;
	friend class Occupy;
	friend class wavefunc;
	friend class Ions;
	friend class Force_Lo;
	friend class Exx_pw;
	friend class ON_Tests;

    energy();
    ~energy();

	//calculate density of states(DOS) and partial density of states(PDOS) and mulliken charge for LCAO base
	void perform_dos(
		const psi::Psi<double> *psid, 
		const psi::Psi<std::complex<double>> *psi, 
		LCAO_Hamilt &uhm,
		const elecstate::ElecState* pelec);
	//calculate density of states(DOS) for PW base
    void perform_dos_pw(const elecstate::ElecState* pelec);

    double etot;    	   // the total energy of the solid
    double ef;             // the fermi energy
    double ef_up;		   // spin up fermi energy
    double ef_dw;		   // spin down fermi energy
	int printe;			   // print energy every 'printe' electron iteration.
	int iter;

	//below defines the bandgap:
	double bandgap;        // bandgap = E_{lumo} - E_{homo}
	double bandgap_up;      
	double bandgap_dw;
	
	public:

    // The variables needed to compute the energies
    double etot_old;
	double etot_harris;	   // Harris-Weinert-Foulkers energy
    double eband;          // the band energy
    double deband;         // correction for variational energy
	double deband_harris;
    double descf;
	double esol_el;		   //correction for implicit solvation energy
	double esol_cav;		   //correction for implicit solvation energy

    double etxcc;          // the nlcc exchange and correlation
	double etxc;
	double vtxc;
	double exx;            // the exact exchange energy.
	double evdw;			// the vdw energy				// Peize Lin add 2021.03.09

    double demet;          // correction for metals

	int out_dos;			// control dos calculation
	int out_band;                    // control band calculation  pengfei 2014-10-13
	int out_proj_band;                    // control projected band calculation  jiyy add 2022-05-11
	const std::string Name_Angular[5][11] = 
	{
    	{"s"},
    	{"px", "py", "pz"},
    	{"d3z^2-r^2", "dxy", "dxz", "dx^2-y^2", "dyz"},
    	{"f5z^2-3r^2", "f5xz^2-xr^2", "f5yz^2-yr^2", "fzx^2-zy^2", "fxyz", "fx^3-3*xy^2", "f3yx^2-y^3"},
    	{"g1", "g2", "g3", "g4", "g5", "g6", "g7", "g8", "g9"}
	};          // name of atomic orbital    jiyy add 2022-05-10

	double dos_emin_ev;
	double dos_emax_ev;
	double dos_edelta_ev;
	double dos_scale;
	double bcoeff;

    double ewald();

    //=========================================================
    // delta_escf = -\int \delta rho(r) V_scf
    // Calculate the difference between Hartree and XC energy
    // at first order in the charge density difference \delta
    // rho
    //=========================================================
	void calculate_etot(void);
	void print_etot(const bool converged, const int &iter, 
	const double &scf_thr, const double &duration, const double &pw_diag_thr=0, const double &avg_iter=0, bool print = true);

	void print_occ(const elecstate::ElecState* pelec);

	void print_format(const std::string &name, const double &value);

	void print_orbital_file();

	void calculate_harris(); //mohan add 2012-06-05

    double delta_e(const elecstate::ElecState* pelec);

    void delta_escf(const elecstate::ElecState* pelec);

	void cal_bandgap(const elecstate::ElecState* pelec);
	void cal_bandgap_updw(const elecstate::ElecState* pelec);

	void cal_converged(elecstate::ElecState* pelec);
	ModuleBase::matrix vnew;
	bool vnew_exist = false;

    
    void set_exx();       // Peize Lin add 2016-12-03
    
};

#endif //energy

