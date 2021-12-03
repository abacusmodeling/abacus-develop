# STRU file
- [Examples](#examples)
    - [no latname](#no-latname)
    - [latname fcc](#latname-fcc)
- [Structure of the file](#structure-of-the-file)

    [back to main page](../README.md)
## Examples

The `STRU` file contains the information about the lattice geometry, the name(s) and/or location(s) of the pseudopotential and numerical orbital files, as well as the structural information about the system.
We supply two ways of specifying the lattice geometry. Below are two examples of the `STRU` file for the same system:

### No latname
For this example, no need to supply any input to the variable `latname` in the INPUT file. (See [input parameters](input-main.md#latname).)
```
ATOMIC_SPECIES
Si 28.00 Si_ONCV_PBE-1.0.upf // label; mass; pseudo_file

NUMERICAL_ORBITAL
Si_gga_8au_60Ry_2s2p1d.orb //numerical_orbital_file

LATTICE_CONSTANT
10.2 // lattice scaling factor (Bohr)

LATTICE_VECTORS
0.5 0.5 0.0 // latvec1
0.5 0.0 0.5 // latvec2
0.0 0.5 0.5 // latvec3

ATOMIC_POSITIONS
Direct //Cartesian or Direct coordinate.
Si // Element type
0.0 // magnetism(Start magnetization could also be defined in INPUT more conveniently.If so,there should not be any input of magnetism here)
2 // number of atoms
0.00 0.00 0.00 0 0 0
0.25 0.25 0.25 1 1 1
```
[back to top](#stru-file)

### latname fcc
We see that this example is a silicon fcc lattice. Apart from setting the lattice vectors manually, we also provide another solution where only the Bravais lattice type is required, and the lattice vectors will be generated automatically. For this example, we need to set `latname="fcc"` in the INPUT file. (See [input parameters](input-main.md#latname).)
And the `STRU` file becomes:
```
ATOMIC_SPECIES
Si 28.00 Si_ONCV_PBE-1.0.upf // label; mass; pseudo_file

NUMERICAL_ORBITAL
Si_gga_8au_60Ry_2s2p1d.orb //numerical_orbital_file

LATTICE_CONSTANT
10.2 // lattice scaling factor (Bohr)

ATOMIC_POSITIONS
Direct //Cartesian or Direct coordinate.
Si // Element type
0.0 // magnetism
2 // number of atoms
0.00 0.00 0.00 0 0 0//the position of atoms and other parameter specify by key word
0.25 0.25 0.25 1 1 1
```
The LATTICE_VECTORS section is removed.



[back to top](#stru-file)


## Structure of the file
The `STRU` file contains several sections, and each section must start with a keyword like
`ATOMIC_SPECIES`, `NUMERICAL_ORBITAL`, or `LATTICE_CONSTANT`, etc. to signify what type of
information that comes below.

- ATOMIC_SPECIES

    This section provides information about the type of chemical elements contained the unit cell. Each line defines one type of element. The user should specify the name, the mass, and the pseudopotential file used for each element. The mass of the elment is only used in molecular dynamics simulations. For electronic-structure calculations, the actual mass value isn’t important.
    In the above example, we see information is provided for the element `Si`:
    ```
    Si 28.00 Si_ONCV_PBE-1.0.upf // label; mass; pseudo_file
    ```
    Here `Si_ONCV_PBE-1.0.upf` is the pseudopotential file. When the path is not specified, the file is assumed to be located in work directory. Otherwise, please explicitly specify the location of the pseudopotential files.
-   NUMERICAL_ORBITAL

    Numerical atomic orbitals are only needed for `LCAO` calculations. Thus this section will be neglected in calcultions with plane wave basis. In the above example, numerical atomic orbitals is specified for the element `Si`:
    ```
    Si_gga_8au_60Ry_2s2p1d.orb //numerical_orbital_file
    ```
    ‘Si_gga_8au_60Ry_2s2p1d.orb’ is name of the numerical orbital file. Again here the path is not specified, which means that this file is located in the work directory.
- LATTICE_CONSTANT

    The lattice constant of the system in unit of Bohr.
- LATTICE_VECTORS

    The lattice vectors of the unit cell. It is a 3by3 matrix written in 3 lines. Please note that *the lattice vectors given here are scaled by the lattice constant*. This section must be removed if the type Bravais lattice is specified using the input parameter `latname`. (See [input parameters](input-main.md#latname).)
- LATTICE_PARAMETERS

    This section is only relevant when `latname` (see [input parameters](input-main.md#latname)) is used to specify the Bravais lattice type. The example above is a fcc lattice, where no additional information except the lattice constant is required to determine the geometry of the lattice.
    
    However, for other types of Bravais lattice, other parameters might be necessary. In that case, the section `LATTICE_PARAMETERS` must be present. It contains **one single line** with some parameters (separated by blank space if multiple parameters are needed), where the number of parameters required depends on specific type of lattice.

    The three lattice vectors v1, v2, v3 (in units of lattice constant) are generated in the following way:

    - latname = "sc": the LATTICE_PARAMETERS section is not required:
    ```
        v1 = (1, 0, 0)
        v2 = (0, 1, 0)
        v3 = (0, 0, 1)
    ```
    - latname = "fcc": the LATTICE_PARAMETERS section is not required:
    ```
        v1 = (-0.5, 0, 0.5)
        v2 = (0, 0.5, 0.5)
        v3 = (-0.5, 0.5, 0)
    ```
    - latname = "bcc" : the LATTICE_PARAMETERS section is not required:
    ```
        v1 = (0.5, 0.5, 0.5)
        v2 = (-0.5, 0.5, 0.5)
        v3 = (-0.5, -0.5, 0.5)
    ```
    - latname = "hexagonal" : One single parameter is required in the LATTICE_PARAMETERS section, being the ratio between axis length c/a. Denote it by x then:
    ```
        v1 = (1.0, 0, 0)
        v2 = (-0.5, sqrt(3)/2, 0)
        v3 = (0, 0, x)
    ```
    
    - latname = "trigonal" : One single parameter is required in the LATTICE_PARAMETERS section, which specifies cos&gamma; with &gamma; being the angle between any pair of crystallographic vectors. Denote it by x then:
    ```
        v1 = (tx, -ty, tz)
        v2 = (0, 2ty, tz)
        v3 = (-tx, -ty, tz)
    ```
    where tx=sqrt((1-x)/2), ty=sqrt((1-x)/6), and tz=sqrt((1+2x)/3).

    - latname = "st" (simple tetragonal) : One single parameter is required in the LATTICE_PARAMETERS section, which gives ratio between axis length c/a. Denote it by x then:
    ```
        v1 = (1, 0, 0)
        v2 = (0, 1, 0)
        v3 = (0, 0, x)
    ```

    - latname = "bct" (body-centered tetragonal) : One single parameter is required in the LATTICE_PARAMETERS section, which gives ratio between axis length c/a. Denote it by x then:
    ```
        v1 = (0.5, -0.5, x)
        v2 = (0.5, 0.5, x)
        v3 = (-0.5, -0.5, x)
    ```
    - latname = "so" (simple orthorhombic) : Two parameters are required in the LATTICE_PARAMETERS section, which gives ratios between axis length b/a and c/a. Denote them by x, y then:
    ```
        v1 = (1, 0, 0)
        v2 = (0, x, 0)
        v3 = (0, 0, y)
    ```
    - latname = "baco" (base-centered orthorhombic) : Two parameters are required in the LATTICE_PARAMETERS section, which gives ratios between axis length b/a and c/a. Denote them by x, y then:
    ```
        v1 = (0.5, x/2, 0)
        v2 = (-0.5, x/2, 0)
        v3 = (0, 0, y)
    ```

    - latname = "fco" (face-centered orthorhombic) : Two parameters are required in the LATTICE_PARAMETERS section, which gives ratios between axis length b/a and c/a. Denote them by x, y then:
    ```
        v1 = (0.5, 0, y/2)
        v2 = (0.5, x/2, 0)
        v3 = (0, x/2, y/2)
    ```

    - latname = "bco" (body-centered orthorhombic) : Two parameters are required in the LATTICE_PARAMETERS section, which gives ratios between lattice vector length b/a and c/a. Denote them by x, y then:
    ```
        v1 = (0.5, x/2, y/2)
        v2 = (-0.5, x/2, y/2)
        v3 = (-0.5, -x/2, y/2)
    ```

    - latname = "sm" (simple monoclinic) : Three parameters are required in the LATTICE_PARAMETERS section, which are the ratios of lattice vector length b/a, c/a as well as the cosine of angle between axis a and b. Denote them by x, y, z then:
    ```
        v1 = (1, 0, 0)
        v2 = (x*z, x*sqrt(1-z^2, 0)
        v3 = (0, 0, y)
    ```
    - latname = "bacm" (base-centered monoclinic) : Three parameters are required in the LATTICE_PARAMETERS section, which are the ratios of lattice vector length b/a, c/a as well as the cosine of angle between axis a and b. Denote them by x, y, z then:
    ```
        v1 = (0.5, 0, -y/2)
        v2 = (x*z, x*sqrt(1-z^2), 0)
        v3 = (0.5, 0, y/2)
    ```

    - latname = "triclinic" : Five parameters are required in the LATTICE_PARAMETERS section, namely the ratios b/a, c/a; the cosines of angle ab, ac, bc. Denote them by x,y,m,n,l, then:
    ```
        v1 = (1, 0, 0)
        v2 = (x*m, x*sqrt(1-m^2), 0)
        v3 = (y*n, y*(l-n*m/sqrt(1-m^2)), y*fac)
    ```
    where fac=sqrt(1+2*m\*n\*l-m<sup>2</sup>-n<sup>2</sup>-l<sup>2</sup>)/sqrt(1-m<sup>2</sup>)
- ATOMIC_POSITIONS
    
    This section specifies the positions and other information of individual atoms. The first line signifies whether atom positions are given in `Cartesian` or `Direct` coordinates.
    
    The following three lines tells the elemental type (`Si`), the initial magnetic moment (`0.0`), and the number of atoms for this particular element (`2`) repsectively. Notice this  magnetic moment will be a default value for every atom of this type but will be overrided if one define it for each atom by keyword(see below).
    
    The last two lines in this example are the coordinates of atomic positions. There are three numbers in each line, which specifies the atomic positions, following by other parameters marked by keywords.

    Several other parameters could be defined after the atom position using key word :
    - `m` or NO key word: three numbers, which take value in 0 or 1, control how the atom move in geometry relaxation calculations. The numbers `0 0 0` following the coordinates of the first atom means this atom are *not allowed* to move in all three directions, and the numbers `1 1 1` following the coordinates of the second atom means this atom *can* move in all three directions.  
    - `v` or `vel` or `velocity`: set the three components of initial velocity of atoms in geometry relaxation calculations.
    - `mag` or `magmom` : set the start magnetization for each atom. In colinear case only one number should be given. In non-colinear case one have two choice:either set one number for the norm of magnetization here and specify two polar angle later, or set three number for the xyz commponent of magnetization here.
    - `angle1`: in non-colinear case, specify the angle between c-axis and real spin, in angle measure instead of radian measure
    - `angle2`: in non-colinear case, specify angle between a-axis and real spin in projection in ab-plane , in angle measure instead of radian measure

    e.g.:
    ```
    Fe
    1
    2
    0.0 0.0 0.0 0 0 0 angle1 90 angle2 0
    0.5 0.5 0.5 0 0 0 angle1 90 angle2 180
    ```

[back to top](#stru-file)