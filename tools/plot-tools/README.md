<!--
 * @Date: 2021-08-21 21:58:06
 * @LastEditors: jiyuyang
 * @LastEditTime: 2022-01-03 17:21:08
 * @Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
-->

# Plotting tool for ABACUS

## Requirements
- [NumPy](https://numpy.org/)
- [Matplotlib](https://matplotlib.org/)
- [lxml](https://lxml.de/)
- [setuptools](https://setuptools.pypa.io/en/latest/index.html)

## Installation
Use the following command to install this tool:
```shell
python setup.py install
```

## Usage
There are two ways to use this tool:
1. Specify parameters in `band.py` or `dos.py` directly, and then `python band.py` or `python dos.py`. And you can also import module in your own script e.g. `from abacus_plot.band import Band`. (Recommend)
2. Command-line tools are also supported in this tool. In this way, you need prepare an input file and execute some commands (see below). You can use `abacus-plot -h` to check command-line information

### Band Structure
First, prepare a file named 'config.json' in json format:
```json
{
	"bandfile" : "BANDS_1.dat",
	"efermi" : 6.585653952007503,
	"energy_range"	: [-1.5, 6],
	"kptfile" : "KPT"
}
```
|    Property    |                                    Type                                     |                                                                                             Note                                                                                             |
| :------------: | :-------------------------------------------------------------------------: | :------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------: |
|   *bandfile*   |                            `str` or `List[str]`                             |                                       Bands data file output from ABACUS. 'BANDS_*.dat' for band structure and 'PBANDS_*' for projected band structure                                       |
|    *efermi*    |                          `float` or `List[float]`                           |                                                                                      Fermi level in eV                                                                                       |
| *energy_range* |                                   `list`                                    |                                                                                    Range of energy in eV                                                                                     |
|    *shift*     |                                   `bool`                                    |                                           If set `'true'`, it will evaluate band gap (only for semiconductors and insulators). Default: `'false'`                                            |
|    *index*     | `List[int]` or `Dict[str, List[int]]` or `Dict[str, Dict[int, List[int]]]]` |                     Extract PBANDS of each atom from 'PBANDS_*' file e.g. [index_0, index_1...] or {index_0:[l_0, l_1, ...], ...}, [index_0:{l_0:[m_0, m_1, ...], ...}, ...]                      |
|  *atom_index*  | `List[int]` or `Dict[str, List[int]]` or `Dict[str, Dict[int, List[int]]]]` | Extract PBANDS of each atom with same atom_index from 'PBANDS_*' file e.g. [atom_index_0, atom_index_1...] or {atom_index_0:[l_0, l_1, ...], ...}, [atom_index_0:{l_0:[m_0, m_1, ...], ...}, ...] |
|   *species*    | `List[str]` or `Dict[str, List[int]]` or `Dict[str, Dict[int, List[int]]]]` |              Extract PBANDS of each atom with same species from 'PBANDS_*' file e.g. [elem_0, elem_1...] or {elem_0:[l_0, l_1, ...], ...}, [elem_0:{l_0:[m_0, m_1, ...], ...}, ...]               |
|   *kptfile*    |                                    `str`                                    |                                                                       K-points file with `Line` mode in ABACUS format                                                                        |
|    *label*     |                            `str` or `List[str]`                             |                                                                                   Label of band structure                                                                                    |
|    *color*     |                            `str` or `List[str]`                             |                                                                                   Color of band structure                                                                                    |

If you want to plot band structure of `nspin=2` or compare two band structure on same k-path, you can set *bandfile*, *efermi*, *label*, *color* in `list` type. 

The *kptfile* should be as follows, and notes after `#` will be set as k-points label automatically.
```shell
K_POINTS
7
Line
0.5     0.5     0.5     20      #L
0.0     0.0     0.0     20      #G
0.5     0.0     0.5     10      #X
0.5     0.25    0.75    15      #W
0.5     0.5     0.5     10      #L
0.375   0.375   0.75    20      #K
0.0     0.0     0.0     1       #G
```
Then, use the following command to plot band structure:
```shell
abacus-plot -b 
```

Then, the following command will plot projected band structure and figures output to directory `PBAND*_FIG` ($*=1$ for $nspin=1$, $*=2$ for $nspin=2$):
```shell
abacus-plot -d -p 
```

Then, the following command will output parsed partial DOS to directory `PBAND*_FILE` ($*=1$ for $nspin=1$, $*=2$ for $nspin=2$):
```shell
abacus-plot -d -o
```

### DOS
First, prepare a file named 'config.json' in json format:
```json
{
	"pdosfile": "PDOS",
	"efermi": 6.585653952007503,
	"energy_range": [
		-5,
		7
	],
	"dos_range": [
		0,
		5
	],
	"species": {
			"C": {
				"0": [
					0
				],
				"1": [
					0,
					1
				]
			},
			"Si": [
				0,
				1
			]
		},
	"pdosfig": "pdos.png"
}
```
If you only want to plot total DOS, you can modify `pdosfile` to `tdosfile` and do not set `species` and `pdosfig`.
|    Property    |                                    Type                                     |                                                                                 Note                                                                                  |
| :------------: | :-------------------------------------------------------------------------: | :-------------------------------------------------------------------------------------------------------------------------------------------------------------------: |
|   *tdosfile*   |                                    `str`                                    |                                                                Total DOS data file output from ABACUS                                                                 |
|   *pdosfile*   |                                    `str`                                    |                                                        Partial DOS data file output from ABACUS in xml format                                                         |
|    *efermi*    |                                   `float`                                   |                                                                           Fermi level in eV                                                                           |
| *energy_range* |                                   `list`                                    |                                                                         Range of energy in eV                                                                         |
|    *shift*     |                                   `bool`                                    |                                                    If set `'true'`, it will evaluate band gap. Default: `'false'`                                                     |
|    *index*     | `List[int]` or `Dict[str, List[int]]` or `Dict[str, Dict[str, List[int]]]]` |                     Extract PDOS of each atom e.g. [index_0, index_1...] or {index_0:[l_0, l_1, ...], ...}, [index_0:{l_0:[m_0, m_1, ...], ...}, ...]                      |
|  *atom_index*  | `List[int]` or `Dict[str, List[int]]` or `Dict[str, Dict[str, List[int]]]]` | Extract PDOS of each atom with same atom_index e.g. [atom_index_0, atom_index_1...] or {atom_index_0:[l_0, l_1, ...], ...}, [atom_index_0:{l_0:[m_0, m_1, ...], ...}, ...] |
|   *species*    | `List[str]` or `Dict[str, List[int]]` or `Dict[str, Dict[str, List[int]]]]` |              Extract PDOS of each atom with same species e.g. [elem_0, elem_1...] or {elem_0:[l_0, l_1, ...], ...}, [elem_0:{l_0:[m_0, m_1, ...], ...}, ...]               |
|   *tdosfig*    |                                    `str`                                    |                                                                      Output picture of total DOS                                                                      |
|   *pdosfig*    |                                    `str`                                    |                                                                     Output picture of partial DOS                                                                     |

Then, the following command will plot total DOS:
```shell
abacus-plot -d 
```

Then, the following command will plot partial DOS:
```shell
abacus-plot -d -p 
```

Then, the following command will output parsed partial DOS to directory `PDOS_FILE`:
```shell
abacus-plot -d -o
```