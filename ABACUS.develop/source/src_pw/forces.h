#ifndef FORCES_H
#define FORCES_H

#include "tools.h"

class Forces
{
public:
	friend class Force_LCAO;
    /* This routine is a driver routine which compute the forces
     * acting on the atoms, the complete forces in plane waves
     * is computed from 4 main parts
     * (1) cal_force_loc: contribution due to local potential.
     * (2) cal_foce_ew: contribution due to ewald potential.
     * (3) cal_force_cc: contributino due to NLCC.
     * (4) cal_nl: contribution due to the non-local pseudopotential.
     * (4) cal_scc: contributino due to incomplete SCF calculation.
     */
    Forces(const int natom);
    ~Forces();

    void init();

    matrix force;
private:

    int nat;
	static double output_acc;

    matrix forcelc;
    matrix forceion;
    matrix forcecc;
    matrix forcenl;
    matrix forcescc;

    void cal_force_loc();
    void cal_force_ew();
    void cal_force_cc();
    void cal_force_nl();
    void cal_force_scc();

    static void print( const string &name, const matrix &f, bool rv=true );
    static void print_to_files( ofstream &ofs, const string &name, const matrix &f );
};

#endif
