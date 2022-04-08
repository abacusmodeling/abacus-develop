#include "iterdiagcon.h"

namespace ModuleHSolver
{

double IterDiagControl::avg_iter = 0.0;
int IterDiagControl::PW_DIAG_NMAX = 30;
double IterDiagControl::PW_DIAG_THR = 1.0e-2;
int IterDiagControl::ntry = 0;
int IterDiagControl::notconv = 0;

}