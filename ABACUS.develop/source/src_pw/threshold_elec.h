#ifndef THRESHOLD_ELEC_H
#define THRESHOLD_ELEC_H

#include "tools.h"

class Threshold_Elec
{

public:

    // constructor and deconstructor
    Threshold_Elec();
    ~Threshold_Elec() {};

protected:

    double dr2;
    bool conv_elec;
    void set_ethr() const;
    void update_ethr(const int &iter);
    void iter_end(ofstream &ofs);
    void print_eigenvalue(ofstream &ofs);
};

#endif
