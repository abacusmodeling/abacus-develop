#ifndef LOOP_CELL_H
#define LOOP_CELL_H
#include<vector>
#include "module_base/matrix.h"
#include "module_base/complexmatrix.h"
#include "module_orbital/ORB_control.h"

class LOOP_cell
{
	public:

	LOOP_cell();
	~LOOP_cell();

	void opt_cell(ORB_control &orb_con);

};

#endif
