<!--
 * @Date: 2021-08-21 21:58:06
 * @LastEditors: jiyuyang
 * @LastEditTime: 2021-08-21 23:14:43
 * @Mail: jiyuyang@mail.ustc.edu.cn, 1041176461@qq.com
-->

# Plotting tool for ABACUS

## Requirements
- [NumPy](https://numpy.org/)
- [Matplotlib](https://matplotlib.org/)
- [lxml](https://lxml.de/)

## Installation
Use the following command to install this tool:
```shell
python setup.py install
```

## Usage
### Band Structure
First, prepare a json file e.g. band-input.json:
```json
{
	"filename" : "BANDS_1.dat",
	"efermi" : 10.39537843955395,
	"energy_range"	: [-1.5, 6],
	"kptfile" : "KPT"
}
```
|Property|Type|Note|
|:--:|:--:|:--:|
|*filename*|`str | List[str]`|Bands data file output from ABACUS|
|*efermi*|`float | List[float]`|Fermi level in eV|
|*energy_range*|`list`| Range of energy in eV|
|*kptfile*|`str`|K-points file in ABACUS format|
|*label*|`str | List[str]`|Label of band structure|
|*color*|`str | List[str]`|Color of band structure|
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
abacus-plot -b band-input.json
```

### DOS
First, prepare a json file e.g. dos-input.json:
```json
{
	"pdosfile": "PDOS",
	"efermi": 10.39537843955395,
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
	"pdosfig": "pdos.png",
	"tdosfig": "tdos.png"
}
```
If you only want to plot total DOS, you can modify `pdosfile` to `tdosfile` and do not set `species` and `pdosfig`.
|Property|Type|Note|
|:--:|:--:|:--:|
|*tdosfile*|`str`|Total DOS data file output from ABACUS|
|*pdosfile*|`str`|Partial DOS data file output from ABACUS in xml format|
|*efermi*|`float`|Fermi level in eV|
|*energy_range*|`list`| Range of energy in eV|
|*species*|`List[str] | Dict[str, List[int]] | Dict[str, Dict[str, List[int]]]]`| Three ways to plot partial DOS e.g. ["C", "Si"], {"C":[0, 1], "Si":[0]}, ["C":{"0":[0]}, "Si":{"1":[0, 1]}]|
|*tdosfig*|`str`|Output picture of total DOS|
|*pdosfig*|`str`|Output picture of partial DOS|
Then, the following command will plot both total DOS and partial DOS:
```shell
abacus-plot -b dos-input.json
```