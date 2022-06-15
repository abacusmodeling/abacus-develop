#ifndef STO_FUNC_H
#define STO_FUNC_H

template <typename REAL>
class Sto_Func
{
public:
    Sto_Func();
    ~Sto_Func(){};
    REAL tem; //temperature
    REAL mu; //chemical potential
    REAL Emin, Emax;

public:
    REAL root_fd(REAL e);
    REAL fd(REAL e);
    REAL nroot_fd(REAL e);
    REAL nfd(REAL e);
    REAL nxfd(REAL e);
    REAL fdlnfd(REAL e);
    REAL nfdlnfd(REAL e);
};

#endif