#ifndef VL_IN_PW_H 
#define VL_IN_PW_H 

#include "../module_base/global_function.h"
#include "../module_base/global_variable.h"
#include "../module_base/matrix.h"

class pseudopot_cell_vl
{
public:

	pseudopot_cell_vl();
	~pseudopot_cell_vl();

	void init_vloc(const int &nggm, ModuleBase::matrix &vloc_in);

	ModuleBase::matrix vloc;   //(ntype,ngl),the local potential for each atom type(ntype,ngl)
	bool *numeric; //[ntype], =true

private:

	double *zp;   // (npsx),the charge of the pseudopotential

	void allocate(void);

	// generate vloc for a particular atom type.
	void vloc_of_g( const int &msh, const double *rab, const double *r, const double *vloc_at,
	               const double &zp, double *vloc) const;

	void print_vloc(void) const;

};

#endif // VL_IN_PW 
