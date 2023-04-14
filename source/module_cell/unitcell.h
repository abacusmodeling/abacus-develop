#ifndef UNITCELL_H
#define UNITCELL_H

#include "module_base/global_function.h"
#include "module_base/global_variable.h"
#include "module_base/matrix3.h"
#include "module_base/intarray.h"
#include "module_io/output.h"
#include "module_elecstate/magnetism.h"
#include "atom_spec.h"

#ifdef __LCAO
#include "module_basis/module_ao/ORB_read.h"
#include "setup_nonlocal.h"
#endif

// provide the basic information about unitcell.
class UnitCell
{
public:
    Atom *atoms;

    bool set_atom_flag;//added on 2009-3-8 by mohan
    Magnetism magnet;  // magnetism Yu Liu 2021-07-03
    void cal_ux();
    bool judge_parallel(double a[3],ModuleBase::Vector3<double> b);
	double *atom_mag;
	int n_mag_at;

    int ntype;// number of atom species in UnitCell
    int nat; // total number of atoms of all species in unitcell
    std::string Coordinate; // "Direct" or "Cartesian" or "Cartesian_angstrom"
    std::string latName; // Lattice name
    double lat0; // Lattice constant(bohr)(a.u.)
    double lat0_angstrom;// Lattice constant(angstrom)
    double tpiba;// 2*pi / lat0;
    double tpiba2; // tpiba ^ 2
    double omega;// the volume of the unit cell

    ModuleBase::Matrix3 latvec; // Unitcell lattice vectors
	int *lc;  // Change the lattice vectors or not
	ModuleBase::Vector3<double> a1,a2,a3; // Same as latvec, just at another form.
	ModuleBase::Vector3<double> latcenter; // (a1+a2+a3)/2 the center of vector
    ModuleBase::Matrix3 latvec_supercell; // Supercell lattice vectors
    ModuleBase::Matrix3 G; // reciprocal lattice vector (2pi*inv(R) )
    ModuleBase::Matrix3 GT; // traspose of G
    ModuleBase::Matrix3 GGT; // GGT = G*GT
    ModuleBase::Matrix3 invGGT; // inverse G

    //========================================================
    // relationship between:
    // ntype, it
    // nat, iat
    // atoms[it].na, ia,
    // atoms[it].nw, iw
    //
    // if know it ==> atoms[it].na; atoms[it].nw
    // if know iat ==> it; ia;
    // if know ia, mush have known it ==> iat
    // if know iwt, must have known it, ia ==> iwt
    //========================================================
    int namax;// the max na among all atom species
    int nwmax;// the max nw among all atom species

    int *iat2it; //iat==>it, distinguish a atom belong to which type
    int *iat2ia; //iat==>ia
	int *iwt2iat; // iwt ==> iat.
	int *iwt2iw; // iwt ==> iw, Peize Lin add 2018-07-02
    ModuleBase::IntArray itia2iat;//(it, ia)==>iat, the index in nat, add 2009-3-2 by mohan
    //atom index iat to the first global index for orbital of this atom
    std::vector<int> iat2iwt;
    // indexing tool for find orbital global index from it,ia,iw
    template<typename Tiait>
    inline Tiait itiaiw2iwt(const Tiait &it, const Tiait &ia, const Tiait &iw) const
    {
        return Tiait(this->iat2iwt[this->itia2iat(it, ia)] + iw);
    }

    //========================================================
    // indexing tools for ia and it
    // return true if the last out is reset
    //========================================================
    template<typename Tiat, typename Tiait>
    inline bool iat2iait(const Tiat iat, Tiait *ia, Tiait *it) const
    {
        if (iat >= nat)
        {
            *ia = 0;
            *it = ntype;
            return false;
        }
        *ia = (Tiait)iat2ia[iat];
        *it = (Tiait)iat2it[iat];
        return true;
    }

    template<typename Tiat, typename Tiait>
    inline bool ijat2iaitjajt(const Tiat ijat, Tiait *ia, Tiait *it, Tiait *ja, Tiait *jt) const
    {
        Tiat iat = ijat / nat;
        Tiat jat = ijat % nat;
        iat2iait(iat, ia, it);
        iat2iait(jat, ja, jt);
        return true;
    }

    template<typename Tiait>
    inline bool step_it(Tiait *it) const
    {
        if (++(*it) >= ntype) {
            *it = 0;
            return true;
        }
        return false;
    }

    template<typename Tiait>
    inline bool step_ia(const Tiait it, Tiait *ia) const
    {
        if (++(*ia) >= atoms[it].na) {
            *ia = 0;
            return true;
        }
        return false;
    }

    template<typename Tiait>
    inline bool step_iait(Tiait *ia, Tiait *it) const
    {
        if (step_ia(*it, ia)) {
            return step_it(it);
        }
        return false;
    }

    template<typename Tiait>
    inline bool step_jajtiait(Tiait *ja, Tiait *jt, Tiait *ia, Tiait *it) const
    {
        if (step_iait(ja, jt)) {
            return step_iait(ia, it);
        }
        return false;
    }

    //LiuXh add 20180515
    ModuleBase::Matrix3 G0;
    ModuleBase::Matrix3 GT0;
    ModuleBase::Matrix3 GGT0;
    ModuleBase::Matrix3 invGGT0;
	
    //I'm doing a bad thing here! Will change later
    bool ionic_position_updated = false; //whether the ionic position has been updated
    bool cell_parameter_updated = false; //whether the cell parameters are updated

    //============================================================
	// meshx : max number of mesh point in pseudopotential file
	// natomwfc : number of starting wavefunctions
	// lmax  : Max L used for localized orbital.
	// nmax  : Max N used for localized orbital.
	// lmax_ppwf : Max L of pseudo wave functinos
	// nelec : total number of electrons
	// lmaxmax : revert from INPUT
	//============================================================
	int meshx;
	int natomwfc;
	int lmax;
	int nmax;
	int nmax_total;//mohan add 2009-09-10
	int lmax_ppwf;
	int lmaxmax; // liuyu 2021-07-04
	bool init_vel; // liuyu 2021-07-15
	// double nelec;

private:
    ModuleBase::Matrix3 stress; //calculate stress on the cell

public:
    UnitCell();
    ~UnitCell();
    void print_cell(std::ofstream &ofs)const;
    void print_cell_xyz(const std::string &fn)const;
    void print_cell_cif(const std::string &fn)const;

    void update_pos_tau(const double* pos);
    void update_pos_taud(const ModuleBase::Vector3<double>* posd_in);
    void update_pos_taud(double* posd_in);
    void update_vel(const ModuleBase::Vector3<double>* vel_in);
    void periodic_boundary_adjustment();
    void bcast_atoms_tau();
    bool judge_big_cell(void)const;

    void update_stress(ModuleBase::matrix &scs); //updates stress
    void update_force(ModuleBase::matrix &fcs); //updates force in Atom

    double *atom_mass;
    std::string *atom_label;
    std::string *pseudo_fn;
    std::string *pseudo_type; // pseudopotential types for each elements, sunliang added 2022-09-15. 
    std::string *orbital_fn;   // filenames of orbitals, liuyu add 2022-10-19
    std::string descriptor_file;  // filenames of descriptor_file, liuyu add 2023-04-06

#ifdef __MPI
    void bcast_unitcell(void);
    void bcast_unitcell2(void);
#endif

	void set_iat2itia(void);

    void setup_cell(const std::string &fn, std::ofstream &log);

#ifdef __LCAO
	InfoNonlocal infoNL;//store nonlocal information of lcao, added by zhengdy 2021-09-07
    void read_orb_file(int it, std::string &orb_file, std::ofstream &ofs_running, Atom *atom);
#endif
	int read_atom_species(std::ifstream &ifa, std::ofstream &ofs_running); // read in the atom information for each type of atom
	bool read_atom_positions(std::ifstream &ifpos, std::ofstream &ofs_running, std::ofstream &ofs_warning); // read in atomic positions

    void read_pseudo(ofstream &ofs);
	int find_type(const std::string &label);
	void print_tau(void)const;
	void print_stru_file(const std::string &fn, const int &type=1, const int &level=0)const; // mohan add 2011-03-22
	void check_dtau(void);
    void setup_cell_after_vc(std::ofstream &log); //LiuXh add 20180515
	
	//for constrained vc-relaxation where type of lattice
	//is fixed, adjust the lattice vectors
	void remake_cell();

	// read in pseudopotential from files for each type of atom
	void read_cell_pseudopots(const std::string &fn, std::ofstream &log);

	//================================================================
	// cal_natomwfc : calculate total number of atomic wavefunctions
	// cal_nwfc     : calculate total number of local basis and lmax
	// cal_meshx	: calculate max number of mesh points in pp file
	//================================================================
	void cal_nwfc(std::ofstream &log);
	void cal_meshx();
	void cal_natomwfc(std::ofstream &log); 
	void print_unitcell_pseudo(const std::string &fn);
	bool check_tau(void)const; //mohan add 2011-03-03
	bool if_atoms_can_move()const;
	bool if_cell_can_change()const;
	void setup(const std::string &latname_in,
            const int &ntype_in,
			const int &lmaxmax_in,
			const bool &init_vel_in,
			const std::string &fixed_axes_in);

	void check_structure(double factor);
};

#endif //unitcell class
