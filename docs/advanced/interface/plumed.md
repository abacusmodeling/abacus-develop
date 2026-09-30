# PLUMED

## Introduction

[PLUMED](https://www.plumed.org) is an open-source library for molecular
dynamics that provides collective variables (CVs), biasing potentials
(metadynamics, umbrella sampling, walls, ...) and free-energy analysis tools.
It is driven through a small C API, so the same `plumed.dat` input file can be
reused across different MD engines.

ABACUS supports PLUMED during *ab initio* molecular dynamics (`calculation md`):
on every MD step the atomic positions, cell, masses, potential energy and (when
available) the virial are handed over to PLUMED, and the bias forces returned
by PLUMED are added to the nuclear forces. The interface is driven from the
shared MD loop, so all the algorithms available through `md_type` (`nve`,
`nvt`, ...) are supported.

> **Note:** the current implementation requires a **single MPI rank** (OpenMP
> parallelization within that rank is fully supported). Running PLUMED together
> with multi-rank MPI is not supported yet.

## Installation

Install PLUMED 2.x following the instructions on https://www.plumed.org/doc
(`./configure --prefix=<plumed-install> && make install`), then configure
ABACUS with the prefix of that installation:

```bash
cmake -B build -DPLUMED_ROOT=<plumed-install>
cmake --build build -j
```

When `PLUMED_ROOT` is not given, the feature is disabled and nothing else
changes, so PLUMED is a completely optional dependency.

## How to use

Enable PLUMED in the `INPUT` file and point to the PLUMED input file:

```
plumed                 1
plumed_file            plumed.dat
```

A minimal `plumed.dat` for well-tempered metadynamics along a distance:

```
UNITS LENGTH=A TIME=fs ENERGY=eV
d: DISTANCE ATOMS=1,2
metad: METAD ARG=d PACE=100 HEIGHT=0.05 SIGMA=0.1 BIASFACTOR=10 TEMP=300 FILE=HILLS
PRINT ARG=d,metad.bias FILE=COLVAR STRIDE=10
```

The output files requested by PLUMED (`COLVAR`, `HILLS`, ...) are written into
the working directory of the run and are flushed at every MD step, so they can
be followed while the simulation is running. The accumulated bias can be
converted into a free-energy surface with the PLUMED command line tool, e.g.

```bash
plumed sum_hills --hills HILLS --mintozero --outfile fes.dat
```

### Units

The interface declares the ABACUS atomic units (Bohr, Hartree and the atomic
unit of time) to PLUMED, and all the arrays are exchanged in these units;
PLUMED converts them internally, so no manual conversion is needed. Use the
`UNITS` directive in `plumed.dat` to select the units of the values written in
the file (for example `UNITS LENGTH=A TIME=fs ENERGY=eV`).
