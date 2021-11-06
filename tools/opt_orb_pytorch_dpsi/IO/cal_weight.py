import IO.read_istate
import numpy as np
import re

def cal_weight(info_weight, stru_file_list=None):
	""" weight[ist][ib] """

	if "bands_file" in info_weight.keys():
		if "bands_range" in info_weight.keys():
			raise IOError('"bands_file" and "bands_range" only once')

		weight = []
		for weight_stru, file_name in zip(info_weight["stru"], info_weight["bands_file"]):
			occ = IO.read_istate.read_istate(file_name)
			weight += [occ_k * weight_stru for occ_k in occ]
		return weight

	elif "bands_range" in info_weight.keys():
		k_weight = read_k_weight(stru_file_list)

		st_weight = []
		for weight_stru, bands_range in zip(info_weight["stru"], info_weight["bands_range"]):
			st_weight.append( np.ones((bands_range,)) * weight_stru )

		weight = []
		for ist,_ in enumerate(k_weight):
			for ik,_ in enumerate(k_weight[ist]):
				weight.append(st_weight[ist] * k_weight[ist][ik])
		return weight

	else:
		raise IOError('"bands_file" and "bands_range" must once')


def read_k_weight(stru_file_list):
	""" weight[ist][ik] """
	weight = []
	for file_name in stru_file_list:
		weight_k = []
		with open(file_name,"r") as file:
			data = re.compile(r"<WEIGHT_OF_KPOINTS>(.+)</WEIGHT_OF_KPOINTS>", re.S).search(file.read()).group(1).split("\n")
			for line in data:
				line = line.strip()
				if line:
					weight_k.append(float(line.split()[-1]))
		weight.append(weight_k)
	return weight