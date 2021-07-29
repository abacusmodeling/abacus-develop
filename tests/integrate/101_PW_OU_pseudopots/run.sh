#!/bin/bash

output_file="OUT.autotest/running_scf.log"

# 1. Test purposes

# 2. Reach for global variables, such as the path of ABACUS,
# "threshold" for comparison, the 'number of processors" for tests

# 3. Run ABACUS, can run multiple tests

# 4. Get key information: energy, occupations, etc.
ewald=`grep E_Ewald $output_file | awk '{print $3}'`
echo "ewald $ewald"

force1=`grep -A 5 "TOTAL-FORCE (eV/Angstrom)" $output_file | tail -n 1 | awk '{print $2}'`
force2=`grep -A 5 "TOTAL-FORCE (eV/Angstrom)" $output_file | tail -n 1 | awk '{print $3}'`
force3=`grep -A 5 "TOTAL-FORCE (eV/Angstrom)" $output_file | tail -n 1 | awk '{print $4}'`
echo "atom1_force_x $force1"
echo "atom1_force_y $force2"
echo "atom1_force_z $force3"

# 5. Compare to the benchmark results

# 6. Give conclusions
