#include "ions_move_methods.h"
#include "ions_move_basic.h"
#include "../src_pw/global.h"

Ions_Move_Methods::Ions_Move_Methods(){}
Ions_Move_Methods::~Ions_Move_Methods(){}

void Ions_Move_Methods::allocate()
{
	Ions_Move_Basic::dim = ucell.nat * 3;

	if(MOVE_IONS=="bfgs")
	{
		this->bfgs.allocate();
	}
	else if(MOVE_IONS=="sd")
	{
		this->sd.allocate();
	}
        else if(MOVE_IONS=="cg")
        {
                this->cg.allocate();
        }
        else if(MOVE_IONS=="cg_bfgs")
	{
		this->cg.allocate();
                this->bfgs.allocate();           // added by pengfei  13-8-8
	}
	else
	{
		WARNING("Ions_Move_Methods::init","the parameter MOVE_IONS is not correct.");
	}
	return;
}

//void Ions_Move_Methods::cal_movement(const int &istep, const matrix &f, const double &etot)
void Ions_Move_Methods::cal_movement(const int &istep, const int &force_step, const matrix &f, const double &etot)
{
	TITLE("Ions_Move_Methods","init");

	//Ions_Move_Basic::istep = istep;
    Ions_Move_Basic::istep = force_step;

	if(MOVE_IONS=="bfgs")
	{
		// move_ions
		// output tau
		// check all symmery
		bfgs.start(f, etot);
	}
	else if(MOVE_IONS=="sd")
	{
		sd.start(f, etot);	
	}
	else if(MOVE_IONS=="cg")
	{
		cg.start(f, etot );
	}
        else if(MOVE_IONS=="cg_bfgs")
        {
                cg.start(f, etot );                                                // added by pengfei 13-8-10
        }
	else
	{
		WARNING("Ions_Move_Methods::init","the parameter MOVE_IONS is not correct.");
	}

	// print the atom positions for convinience.
	stringstream ss;

	ss << global_out_dir << "STRU_ION";
	
	if(Ions_Move_Basic::out_stru)
	{	
		ss << istep; 
	}
	
	//xiaohui modify 2015-03-15, cancel outputfile "STRU_ION"
	//ucell.print_stru_file(ss.str(),1);
	ss << "_D";
	ucell.print_stru_file(ss.str(),2);

	return;
}


