//==========================================================
// AUTHOR : mohan
// DATE : 2011-06-12
//==========================================================
#ifndef DC_DRIV_H
#define DC_DRIV_H


// DC stands for 'Divide and Conqure' and related methods.
// A DC class should contain the following important informations:
//
// (1) A 'reading' routine which read all the needed information,
// especially the dividsion method

// (2) A 'divide_frag' routine which is in charge of divided pieces,
// setup the information for each pieces, and let the information
// can be easily accessed.
 
// (3) A 'setup_frag' routine to setup the information for each
// fragment or others.

// (4) A 'solve_each' routine which can be used to solve the
// Kohn-Sham equation in each 'piece', to get either charge density,
// potential, or even wave functions as new basis, or part of
// the Hamiltonian or Overlap matrix.

// (4) A 'solve_interaction' to solve the interaction betweeen 
// 'pieces', the results may turn out to be part of the Hamiltonian
// matrix.

// (5) A 'solve_whole' routine which can be used solve the terms
// which must be calculated in the whole systen, such as the 
// Hartree tern, the global plane wave can be genreated here.

// (6) A 'conqure' routine to deal with the final information,
// either making the final H and S matrix, or patching each
// part of charge density together to calculate the error,
// Global Kohn-Sham equation may be solved here or not,
// depends which method we are using.

class DC_Driv
{
	public:
	
	DC_Driv();
	~DC_Driv();

	void init();

	private:

	void reading();
	void divide_frag();
	void setup_frag();
	void solve_eachf();
	void solve_interaction();
	void solve_whole();
	void conqure();


};

#endif
