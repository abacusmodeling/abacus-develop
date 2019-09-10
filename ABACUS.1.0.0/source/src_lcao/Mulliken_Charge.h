#ifndef Mulliken_Charge_H
#define Mulliken_Charge_H

//#include "../src_pw/tools.h"
#include "../src_pw/tools.h"
#include "use_overlap_table.h"
#include "sltk_grid_driver.h"
#include "lcao_matrix.h"
#include "lcao_matrix.h"
#include "../src_lcao/global_fp.h"
#include "../src_lcao/wfc_dm_2d.h"
#include "../src_global/lapack_connector.h"
#include "../src_global/scalapack_connector.h"
#include "../src_global/matrix.h"
#include "../src_global/complexmatrix.h"
#include <vector>
#include "../src_external/src_pdiag/pdiag_double.h"
#include "../src_external/src_pdiag/GenELPA.h"



class Mulliken_Charge
{
public:

	 Mulliken_Charge();
	~Mulliken_Charge();
              
            double**  DecMulP ;
	          double***  ADecMulP ;
		  std::vector<matrix> Mullik;
                         double**  MecMulP ;
                     std::vector<ComplexMatrix> Mulliken;
                      //PDOS_k.resize(kv.nks);
                         
                     
                                 
	 void  cal_mulliken(void);

	void   stdout_mulliken(void);
	

private:

};
#endif