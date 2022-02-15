#ifndef LOOP_CELL_H
#define LOOP_CELL_H
#include<vector>
#include "module_base/matrix.h"
#include "module_base/complexmatrix.h"

class LOOP_cell
{
	public:

	LOOP_cell();
	~LOOP_cell();

	void opt_cell(std::vector<ModuleBase::matrix> &wfc_gamma,
    std::vector<ModuleBase::ComplexMatrix> &wfc_k);

};

#endif
