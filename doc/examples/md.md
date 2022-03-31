# Molecular dynamics

[back to main page](../../README.md)

A typical INPUT file for MD simulation is given
below, which contains additional keywords that need to added.
```
INPUT_PARAMETERS
#Parameters (1.General)
pseudo_dir          ./
ntype               1
nbands              600
pseudo_type         upf
gamma_only          1

#Parameters (2.Methods)
calculation         md
symmetry            0

out_level           m
relax_method         cg

smearing_method            gaussian
smearing_sigma               0.02
#Parameters (3.PW)
ecutwfc             30
scf_thr_rho                 1e-5
scf_nmax               100

#Parameters (5.LCAO)
basis_type          lcao
mixing_beta         0.4
chg_extrap       second-order

md_nstep          10   // md steps
md_type           1    //choose ensemble
md_dt               1    //time step
md_tfirst           700  //the first target temperature
md_restart            0    //whether restart md
md_dumpfreq       10   //The period to dump MD information
```

These MD parameters means that ABACUS will use NVT ensemble with Nosé-hoover themostat; the time step is 1fs, and target temperature is 700K; start renew without restart file, set the mass of themostat as 1g/mol, and calculate the MSD and diffusion coefficent from first step.

Note: *Please turn off symmetry when do MD simulation.*

- md_type : -1, FIRE; 0, NVE; 1, NHC; 2, LGV; 3, ADS; 4, MSST
- md_dt : time step in md simulation (fs)
- md_tfirst : target temperature in md simulation(K), you should set parameter md_tlast when you want to change temperature during md simulation.
- md_restart : 0, no need of restart ; 1, restart with restart file, you must repalce STRU file with STRU_MD before you run the restart task.
- md_dumpfreq : frequency for output consequence of md simulation

The STRU file is:
```
ATOMIC_SPECIES
Sn 118.71  ./Sn.pz-bhs.UPF 

LATTICE_CONSTANT
23.45378

NUMERICAL_ORBITAL
./Sn_pz-bhs_8.0au_16Ry_2s2p1d

LATTICE_VECTORS
1.00000       0.0000000      0.000000
0.000000      1.0000000      0.00000
0.00000       0.00000000     1.00000

ATOMIC_POSITIONS
Direct 

Sn
0.0
64
0	0	0	1	1	1
0	0.25	0.25	1	1	1
0.25	0	0.25	1	1	1
0.25	0.25	0	1	1	1
0.375	0.125	0.375	1	1	1
0.125	0.125	0.125	1	1	1
0.125	0.375	0.375	1	1	1
0.375	0.375	0.125	1	1	1
0.5	0	0	1	1	1
0.5	0.25	0.25	1	1	1
0.75	0	0.25	1	1	1
0.75	0.25	0	1	1	1
0.875	0.125	0.375	1	1	1
0.625	0.125	0.125	1	1	1
0.625	0.375	0.375	1	1	1
0.875	0.375	0.125	1	1	1
0	0.5	0	1	1	1
0	0.75	0.25	1	1	1
0.25	0.5	0.25	1	1	1
0.25	0.75	0	1	1	1
0.375	0.625	0.375	1	1	1
0.125	0.625	0.125	1	1	1
0.125	0.875	0.375	1	1	1
0.375	0.875	0.125	1	1	1
0.5	0.5	0	1	1	1
0.5	0.75	0.25	1	1	1
0.75	0.5	0.25	1	1	1
0.75	0.75	0	1	1	1
0.875	0.625	0.375	1	1	1
0.625	0.625	0.125	1	1	1
0.625	0.875	0.375	1	1	1
0.875	0.875	0.125	1	1	1
0	0	0.5	1	1	1
0	0.25	0.75	1	1	1
0.25	0	0.75	1	1	1
0.25	0.25	0.5	1	1	1
0.375	0.125	0.875	1	1	1
0.125	0.125	0.625	1	1	1
0.125	0.375	0.875	1	1	1
0.375	0.375	0.625	1	1	1
0.5	0	0.5	1	1	1
0.5	0.25	0.75	1	1	1
0.75	0	0.75	1	1	1
0.75	0.25	0.5	1	1	1
0.875	0.125	0.875	1	1	1
0.625	0.125	0.625	1	1	1
0.625	0.375	0.875	1	1	1
0.875	0.375	0.625	1	1	1
0	0.5	0.5	1	1	1
0	0.75	0.75	1	1	1
0.25	0.5	0.75	1	1	1
0.25	0.75	0.5	1	1	1
0.375	0.625	0.875	1	1	1
0.125	0.625	0.625	1	1	1
0.125	0.875	0.875	1	1	1
0.375	0.875	0.625	1	1	1
0.5	0.5	0.5	1	1	1
0.5	0.75	0.75	1	1	1
0.75	0.5	0.75	1	1	1
0.75	0.75	0.5	1	1	1
0.875	0.625	0.875	1	1	1
0.625	0.625	0.625	1	1	1
0.625	0.875	0.875	1	1	1
0.875	0.875	0.625	1	1	1
```

The KPT file is:
```
K_POINTS
0
Gamma
1 1 1 0 0 0
```

Run the program, and see results in the output directory. The following files are about MD:
- STRU_MD_$num: optimized structures in direct coordinate
- Restart_md.dat: output the information of md for restart
- If you want to restart md, you must replace the STRU with STRU_MD_$num.

MD information can be found in file running_md.log or in file MD_OUT

```
--------------------------------------------------
Molecular Dynamics (NVT) STEP +10
--------------------------------------------------
--------------------------------------------------
SUMMARY OF NVT CALCULATION
--------------------------------------------------
NVT Conservation     : -450.943151 (Rydberg) //total energy of system
NVT Temperature      : +689.931183 (K) //temperature at this second
NVT Kinetic energy   :  +0.419483 (Rydberg) //kinetic energy of system
NVT Potential energy : -452.902847 (Rydberg) //potential energy of system
```

Check these information to confirm whether temperature and conservation is steady.

[back to top](#molecular-dynamics)