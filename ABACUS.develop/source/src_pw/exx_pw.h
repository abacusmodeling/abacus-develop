#ifndef EXX_PW_H
#define EXX_PW_H

#include "tools.h"
#include "exx_div.h"

class Exx_pw
{
	public:

	Exx_pw();
	~Exx_pw();

	void init(bool start_in);
	void get_exx(void);
	void get_exx2(void);
	void vxx_psi(complex<double>* psi_in, complex<double> *hpsi);
	int ik_now;

	private:

	Exx_Divergence exxdiv;
	double alpha; // 0.25 for PBE0.
	bool start;

};

// If we ignore the divergence term, the error intruduced in total energy
// is proportional to 1/(NQS*Omega)

// other tricks:
// dexx = fock1 - 0.5 * ( fock + fock2 ) must be converged.
// fock0 = sum <phi|Vx(phi)|phi>
// fock1 = sum <psi|Vx(phi)|psi>
// fock2 = sum <psi|Vx(psi)|psi>

// -------------- Level 1 -----------------------
// Target : Implement the H|psi> including exact
// exchange Vx | psi >

// * need to overview of CG method, see how to do
// h|psi>

// * need to calculate Vx, which is non-local.
// for eack k point.

// * calculate the pair density for a k point
// n_{nk,mq}(r) = phi_{mq}(r) * chi_{nk}(r)

// * solve the singularity problem.

// -------------- Level 2 -----------------------
// Target : Can calculate the total energy
// include Ex.

// -------------- Level 3 -----------------------
// Target : Implement 0.25 Ex, and PBE0.
// in potential, alphs Ex + (1-alpha) Ex(PBE)
// Ec. And the potential.

// -------------- Level 4 -----------------------
// Target : Parallel PBE0.

// -------------- Level 5 -----------------------
// Target : Force including PBE0.

// -------------- Level 6 -----------------------
// Target : consider HSE03 or B3LYP.
#endif
