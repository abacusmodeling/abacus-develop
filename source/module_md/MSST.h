#ifndef MSST_H
#define MSST_H

#include "verlet.h"

class MSST : public Verlet
{
public:
    MSST(MD_parameters& MD_para_in, UnitCell_pseudo &unit_in);
    ~MSST();

    void setup();
    void first_half();
    void second_half();
    void outputMD();
};

#endif