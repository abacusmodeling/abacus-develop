//==========================================================
// AUTHOR : mohan
// DATE : 2008-11-5
//==========================================================
#ifndef EXIMPORT_H
#define EXIMPORT_H
#include <fstream>

using namespace std;

//==========================================================
// CLASS : 
// NAME : eximport(control output,input datas)
//==========================================================
class eximport
{
public:

	eximport();
	~eximport();

//==========================================================
// MEMBER FUNCTION : 
// NAME : write_data(write the information to be compared)
//==========================================================
	void write_data
	(
		const string &fn, // file name of output file
		const string &type// choice : "all" "evc" "band"
						  // "energy" "charge"
	);

//==========================================================
// MEMBER FUNCTION : 
// NAME : read_data(read in datas) 
//==========================================================
	void read_data
	(
		const string &fn // file name of input file
	);

//==========================================================
// MEMBER FUNCTION : print_data
// DO : print data we need.
//==========================================================
	void print_data
	(
		const string &fn
	)const;
	
//==========================================================
// MEMEBER FUNCTION : fir_wf 
// DO : output wave functions.
//==========================================================
	void fir_wf
	(
		ComplexMatrix *wf,// address of wave functions. 
		const int npsi, // number of wave funcs to output.
		const string &fn // file address.
	);

//==========================================================
// MEMBER FUNCTION : sec_wf (2version)
// DO : readin wave functions
//==========================================================
	bool sec_wf
	(
		std::complex<double> ***phi,// address of wave functions. 
		const int npsi, // how many wave functions read in.
		const string &fn // file address.
	)const;

	bool sec_wf
	(
		ComplexMatrix *psi,// address of wave functions.
		const int npsi,// how many wave functions read in.
		const string &fn// file address.
	)const;

//==========================================================
// MEMBER FUNCTION : out_gspcae_wan
// DO : output gspace wave functions.
//==========================================================
	void out_gspace_wan
	(
		const ComplexMatrix *psi,// address of wave funcs.
		const int iw, // index of wave funs to export.
		const string &file_name// file address.
	);

//==========================================================
// MEMBER FUNCTION : 
// DO :
//==========================================================
	void nscf_chgfile
	(
		const string &chr_file// file address.
	);

//==========================================================
// MEMBER FUNCTION : nscf_band
// DO : output bands in nscf case.
//==========================================================
	static void nscf_band
	(
		const string &out_band_dir,// file address.
		const int &nband
	);

//==========================================================
// MEMBER FUNCTION : calculate_dos
// DO : output DOS in nscf case.
// RETURN : = TRUE (seccess)
//  	    = FALSE (failure)
//==========================================================
	bool calculate_dos
	(
 		const string &fa, // file address.
		const double de_ev, // delta energy in ev
		const double emax_ev, // maximal energy in ev.
		const double emin_ev, // minimal energy in ev.
		const int nks, // number of k points included.
		const int nbands, // number of nbands included.
		const matrix &et// store energy for each k point
		                    // and each band
	);


#ifdef __MPI
//==========================================================
// MEMBER FUNCTION : out_charge_mpi
// DO : output charge density in parallel case.
//
// MEMBER FUNCTION : in_charge_mpi
// DO : input charge density in parallel case.
//==========================================================
	void out_charge_mpi
	(
		const string &dir,
		double* rho_in
	);

	void in_charge_mpi
	(
		const string &dir
	);
#endif

//private:

	string name;

//==========================================================
// MEMBER FUNCTION  : out_wannier
// DO :  output wannier parameters.
//
// MEMBER FUNCTION  : in_wannier
// DO :  read in corresponding wannier parameters.
//
// MEMBER VARIABLES : 
// NAME : basis( which type of basis we use )
//==========================================================
	void out_wannier(ofstream &out);
	void in_wannier(ifstream &in);

	string basis;

//==========================================================
// MEMBER FUNCTION : out_input
// DO : out put information about the control data in this 
//      run.
//
// MEMBER FUNCTION : in_input
// DO : read in the corresponding control data.
//
// MEMBER VARIABLES : 
// NAME : latname(lattice name).
// NAME : calculation("scf" or "nscf" or "relax")
// NAME : ecutwfc(energy cutoff for wave functions).
// NAME : nband(number of bands).
// NAME : tr2(convergence of electronic relaxation).
// NAME : nx(points in x axis of wave function FFT grid).
// NAME : ny(points in y axis of wave function FFT grid).
// NAME : nz(points in z axis of wave function FFT grid).
// NAME : nxyz(total points of wave functions FFT grid).
// NAME : startingpot(starting potential used).
// NAME : Mixing_beta(parameter for charge mixing);
//==========================================================
	void out_input(ofstream &out);
	void in_input(ifstream &in);
	
	string latname;
	string calculation;
	double ecutwfc;
	int nband;
	double tr2;
	int nx;
	int ny;
	int nz;
	int nxyz;
	string startingpot;
	double Mixing_beta;


//==========================================================
// MEMBER FUNCTION : out_kpoints 
// DO : out put data about k points.
//
// MEMBER FUNCTION : in_kpoints
// DO : read in data about k points.
//
// MEMBER VARIABLES :
// NAME : nks(number of k points)
// NAME : ngk(number of G points for each k)
// NAME : kvector(cartesian coordinate of each k point)
// NAME : qtot(number of total q points, q = K+G)
//==========================================================
	void out_kpoints(ofstream &out);
	void out_planewave(ofstream &out);
	void out_igk(ofstream &out);
	void in_kpoints(ifstream &in);

	int nks;
	int *ngk;
	double **kvector;
	int qtot;

//==========================================================
// MEMBER FUNCTION : out_unitcell 
// DO : output information about unitcell.
//
// MEMBER FUNCTION : in_unitcell
// DO : readin corresponding information about unticell.
//
// MEMBER VARIABLES :
// NAME : lat0(scaling factor for dimension 1)
// NAME : latvec(3 real space lattice vectors)
// NAME : ntype(number of atom species)
// NAME : na(number of atoms for each type)
//==========================================================
	void out_unitcell(ofstream &out);
	void in_unitcell(ifstream &in);

	double lat0;	
	double **latvec;
	int ntype;
	int *na;

//==========================================================
// MEMBER FUNCTION : out_charge
// DO : output information about charge density.
//
// MEMBER FUNCTION : in_charge
// DO : input corresponding charge density.
//==========================================================	
	void out_charge(ofstream &out);
	void in_charge(ifstream &in);

//==========================================================
// MEMBER FUNCTION : 
// NAME : out_band(output information about bands calculation)
// NAME : in_band(input corresponding bands calculation 
// information)
//
// MEMBER VARIABLES :
// NAME : band_energy()
// NAME : omega(volumn of cell)
// NAME : ncx(points in x axis of charge FFT grid)
// NAME : ncy(points in y axis of charge FFT grid)
// NAME : ncz(points in z axis of charge FFT grid)
// NAME : rho_nr(number of spin,first dimension of charge)
// NAME : rho_nc(total points of charge FFT grid,rho_nc = ncxyz)
// NAME : rho(charge density)
//==========================================================
	void out_band(ofstream &out);
	void in_band(ifstream &in);

	double **band_energy;
	double omega;
	int ncx;
	int ncy;
	int ncz;
	int rho_nr;
	int rho_nc;
	double *rho;
	
//==========================================================
// MEMBER FUNCTION : out_energy  
// DO : output information about energy,like total energy,
//      band energy , Hartree energy , exchange-correlation
//      energy, ewald energy...
//
// MEMBER FUNCTION : in_energy
// DO : readin corresponding energy information as mentioned
//      in out_energy function.
//
// MEMBER VARIABLES :
// NAME : iter(number of iteration)
// NAME : etot(total energy in Ry)
// NAME : one_electron( one electron energy in Ry)
// NAME : hartree( hartree energy in Ry)
// NAME : xc( exchange-correlation energy in Ry)
// NAME : ewald( ewald energy in Ry)
//==========================================================
	void out_energy(ofstream &out);
	void in_energy(ifstream &in);
	
	int iter;
	double etot;
	double eband;
	double one_electron;
	double hartree;
	double xc;
	double ewald;

//==========================================================
// MEMBER FUNCTION : out_evc
// DO : output information about wave functions.
//
// MEMBER FUNCTION : in_evc
// DO : input data mentioned in out_evc.
//
// MEMBER VARIABLES :
// NAME : natomwfc(number fo states calcualted)
// NAME : evc(wave functions)
//==========================================================
	void out_evc(ofstream &out);
	void in_evc(ifstream &in);

	int natomwfc;
	std::complex <double> ***evc;
};

#endif
