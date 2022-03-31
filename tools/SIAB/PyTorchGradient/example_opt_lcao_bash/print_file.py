import json
import utils
import textwrap

def print_file_pw(info,dis):

	with open("INPUT","w") as file:
		file.write(textwrap.dedent(f"""\
			INPUT_PARAMETERS
			pseudo_dir          {info["input"]["pseudo_dir"]}
			wannier_card        INPUTw
			calculation         scf
			ntype               1
			nspin               1
			lmaxmax             {len(info["orbital"])-1}

			symmetry            0
			nbands             	{info["input"]["nbands"]} 

			ecutwfc             {info["input"]["ecut"]}
			scf_thr_rho                 2.0e-8  // about iteration
			scf_nmax               1000

			smearing_method            gauss
			smearing_sigma               {info["input"]["smearing_sigma"]}

			mixing_type         pulay       // about charge mixing
			mixing_beta         0.4
			mixing_ndim         8
			printe				1
			"""))



	with open("INPUTw","w") as file:
		file.write(textwrap.dedent(f"""\
			WANNIER_PARAMETERS
			rcut 10
			out_spillage 2
			"""))



	with open("INPUTs","w") as file:
		file.write(textwrap.dedent(f"""\
			INPUT_ORBITAL_INFORMATION
			<SPHERICAL_BESSEL>
			1           // smooth or not
			0.1         // smearing_sigma
			{info["input"]["ecut"]}       // energy cutoff for spherical bessel functions(Ry)
			{info["input"]["rcut"]}       // cutoff of wavefunctions(a.u.)
			1.0e-12     // tolerence
			</SPHERICAL_BESSEL>
			"""))



	with open("STRU","w") as file:
		file.write("ATOMIC_SPECIES\n")
		file.write(f'{info["input"]["element"]} 1 {info["input"]["pseudo"]}\n')
		file.write(textwrap.dedent(f"""\
			LATTICE_CONSTANT
			{utils.lat0}\n
			LATTICE_VECTORS
			1 0 0
			0 1 0
			0 0 1\n
			ATOMIC_POSITIONS
			Cartesian_angstrom
			"""))
		if info["input"]["element"] in ["Na","Li","K","Ca"]:
			file.write(textwrap.dedent(f"""\
				{info["input"]["element"]}
				0.0
				3
				0 0		0		0 0 0
				0 0		{dis}		0 0 0			
				0 {dis*0.86603}		{dis*0.5}		0 0 0			
				"""))
		else:
			file.write(textwrap.dedent(f"""\
				{info["input"]["element"]}
				0.0
				2
				0 0 0		0 0 0
				0 0 {dis}		0 0 0
				"""))



	with open("KPT","w") as file:
		file.write(textwrap.dedent(f"""\
			K_POINTS
			0
			Gamma
			1 1 1 0 0 0
			"""))



	if utils.sub=="qsub":
		with open("sub.sh","w") as file:
			core = info["exe"]["qsub"][0]*info["exe"]["qsub"][1]
			file.write(textwrap.dedent(f"""\
				#!/bin/bash
				#PBS -q gold5120
				#PBS -l nodes={info["exe"]["qsub"][0]}:ppn={info["exe"]["qsub"][1]}
				#PBS -l walltime=1:00:00
				#PBS -o job.log
				#PBS -e job.err
				ulimit -s unlimited
				cd $PBS_O_WORKDIR
				EXEC={info["exe"]["exe_pw"]}
				mpirun -n {core} -env OMP_NUM_THREADS=1 $EXEC
				"""))
	elif utils.sub=="bsub":
		with open("sub.sh","w") as file:
			core = info["exe"]["qsub"][0]*info["exe"]["qsub"][1]
			file.write(textwrap.dedent(f"""\
				#!/bin/sh
				#BSUB -q renxg
				#BSUB -o job.log -e job.err
				#BSUB -n {core}
				EXEC={info["exe"]["exe_pw"]}
				mpirun -n {core} -env OMP_NUM_THREADS=1 $EXEC
				"""))
	else:
		raise KeyError("utils.sub = ",utils.sub)





def print_file_opt(info,dis):

	with open("input.json","w") as file:
		input = {
			"file_list": [ f'../{info["input"]["element"]}-{info["input"]["rcut"]}-{distance}/test.{utils.lat0}.dat' for distance in dis[info["input"]["element"]] ],
			"info": {
				"Nt_all":	[info["input"]["element"]],
				"Nu":		{info["input"]["element"] : info["orbital"]},
				"Nb_true":	info["input"]["ref_bands"] if isinstance(info["input"]["ref_bands"],list) else [info["input"]["ref_bands"]] * len(dis[info["input"]["element"]]),
				"weight":	[1] * len(dis[info["input"]["element"]]),
				"Rcut":		{info["input"]["element"] : info["input"]["rcut"]},
				"dr":		{info["input"]["element"] : utils.dr},
				"Ecut":		{info["input"]["element"] : info["input"]["ecut"]},
				"lr":		utils.lr
			},
			"C_init_info": {
				"init_from_file":	False
			},
			"V_info": {
				"same_band":		True,
				"init_from_file":	False
			}
		}
		file.write(json.dumps(input,indent=4))
		
		
	if utils.sub=="qsub":
		with open("sub.sh","w") as file:
			file.write(textwrap.dedent(f"""\
				#!/bin/bash
				#PBS -q gold5120
				#PBS -l nodes=1:ppn=1
				#PBS -l walltime=1:00:00
				#PBS -o job.log
				#PBS -e job.err
				ulimit -s unlimited
				cd $PBS_O_WORKDIR
				export OMP_NUM_THREADS=1
				EXEC={info["exe"]["exe_orbital"]}
				python3 $EXEC
				"""))
	elif utils.sub=="bsub":
		with open("sub.sh","w") as file:
			file.write(textwrap.dedent(f"""\
				#!/bin/sh
				#BSUB -q renxg
				#BSUB -o job.log -e job.err
				#BSUB -n 1
				export OMP_NUM_THREADS=1
				EXEC={info["exe"]["exe_orbital"]}
				python3 $EXEC
				"""))
	else:
		raise KeyError("utils.sub = ",utils.sub)