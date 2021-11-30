from util import *
import torch
import itertools
import numpy as np
import re
import copy

def read_file(info,file_list,V_info):
	""" QI[ist][it][il][ib*ia*im,ie]	<\psi|jY> """
	""" SI[ist][it1][it2][il1][il2][ie1,ia1,im1,ia2,im2,ie2]	<jY|jY> """
	""" VI[ist][ib]		<\psi|\psi> """
	info_true = copy.deepcopy(info)
	info_true.Nst = len(file_list)
	info_true.Nt = ND_list(info_true.Nst,element="list()")
	info_true.Na = ND_list(info_true.Nst,element="dict()")
	info_true.Nb = ND_list(info_true.Nst)
	info_true.Nk = ND_list(info_true.Nst)
	info_true.Ne = dict()
	QI=[];	SI=[];	VI=[]

	for ist_true,file_name in enumerate(file_list):
		print(file_name)
		with open(file_name,"r") as file:

			ignore_line(file,4)
			Nt_tmp = int(file.readline().split()[0])
			for it in range(Nt_tmp):
				t_tmp = file.readline().split()[0]
				assert t_tmp in info.Nt_all
				info_true.Nt[ist_true].append( t_tmp )
				info_true.Na[ist_true][t_tmp] = int(file.readline().split()[0])
				ignore_line( file, info_true.Na[ist_true][t_tmp] )
			ignore_line(file,6)
			Nl_ist = int(file.readline().split()[0])+1
			for it,Nl_C in info.Nl.items():
				print(it,Nl_ist,Nl_C)
				assert Nl_ist>=Nl_C
				info_true.Nl[it] = Nl_ist
			info_true.Nk[ist_true] = int(file.readline().split()[0])
			info_true.Nb[ist_true] = int(file.readline().split()[0])
			ignore_line(file,1)
#			Ne_tmp = list(map(int,file.readline().split()[:Nt_tmp]))
#			for it,Ne in zip(info_true.Nt[ist_true],Ne_tmp):
#				assert info_true.Ne.setdefault(it,Ne)==Ne
			Ne_tmp = int(file.readline().split()[0])
			for it in info_true.Nt[ist_true]:
				info_true.Ne[it] = Ne_tmp

	for ist_true,file_name in enumerate(file_list):
		with open(file_name,"r") as file:
			data = re.compile(r"<OVERLAP_Q>(.+)</OVERLAP_Q>", re.S).search(file.read())
			data = map(float,data.group(1).split())
			for ik in range(info_true.Nk[ist_true]):
				print("read QI:",ist_true,ik)
				qi = read_QI(info_true,ist_true,data)
				QI.append( qi )
		with open(file_name,"r") as file:
			data = re.compile(r"<OVERLAP_Sq>(.+)</OVERLAP_Sq>", re.S).search(file.read())
			data = map(float,data.group(1).split())
			for ik in range(info_true.Nk[ist_true]):
				print("read SI:",ist_true,ik)
				si = read_SI(info_true,ist_true,data)
				SI.append( si )
		if V_info["init_from_file"]:
			with open(file_name,"r") as file:
				data = re.compile(r"<OVERLAP_V>(.+)</OVERLAP_V>", re.S).search(file.read())
				data = map(float,data.group(1).split())
		else:
			data = ()
		for ik in range(info_true.Nk[ist_true]):
			print("read VI:",ist_true,ik)
			vi = read_VI(info_true,V_info,ist_true,data)
			VI.append( vi )
	print()

	info_all = copy.deepcopy(info)
	info_all.Nst = sum(info_true.Nk,0)
	repeat_Nk = lambda x: list( itertools.chain.from_iterable( map( lambda x:itertools.repeat(*x), zip(x,info_true.Nk) ) ) )
	info_all.Nt = repeat_Nk(info_true.Nt)
	info_all.Na = repeat_Nk(info_true.Na)
	info_all.Nb = repeat_Nk(info_true.Nb)
	info_all.Ne = info_true.Ne

	return QI,SI,VI,info_all




def read_QI(info,ist,data):
	""" QI[it][il][ib*ia*im,ie]	<\psi|jY> """
	QI = dict()
	for it in info.Nt[ist]:
		QI[it] = ND_list(info.Nl[it])
		for il in range(info.Nl[it]):
			QI[it][il] = torch.zeros((info.Nb[ist],info.Na[ist][it],info.Nm(il),info.Ne[it]), dtype=torch.complex128)
	for ib in range(info.Nb[ist]):
		for it in info.Nt[ist]:
			for ia in range(info.Na[ist][it]):
				for il in range(info.Nl[it]):
					for im in range(info.Nm(il)):
						for ie in range(info.Ne[it]):
							QI[it][il][ib,ia,im,ie] = complex(next(data), next(data))
	for it in info.Nt[ist]:
		for il in range(info.Nl[it]):
			QI[it][il] = QI[it][il].view(-1,info.Ne[it]).conj()
	return QI



def read_SI(info,ist,data):
	""" SI[it1,it2][il1][il2][ie1,ia1,im1,ia2,im2,ie2]	<jY|jY> """
	SI = dict()
	for it1,it2 in itertools.product( info.Nt[ist], info.Nt[ist] ):
		SI[it1,it2] = ND_list(info.Nl[it1],info.Nl[it2])
		for il1,il2 in itertools.product( range(info.Nl[it1]), range(info.Nl[it2]) ):
			SI[it1,it2][il1][il2] = torch.zeros((info.Na[ist][it1],info.Nm(il1),info.Ne[it1],info.Na[ist][it2],info.Nm(il2),info.Ne[it2]), dtype=torch.complex128)
	for it1 in info.Nt[ist]:
		for ia1 in range(info.Na[ist][it1]):
			for il1 in range(info.Nl[it1]):
				for im1 in range(info.Nm(il1)):
					for it2 in info.Nt[ist]:
						for ia2 in range(info.Na[ist][it2]):
							for il2 in range(info.Nl[it2]):
								for im2 in range(info.Nm(il2)):
									for ie1 in range(info.Ne[it1]):
										for ie2 in range(info.Ne[it2]):
											SI[it1,it2][il1][il2][ia1,im1,ie1,ia2,im2,ie2] = complex(next(data), next(data))
#	for it1,it2 in itertools.product( info.Nt[ist], info.Nt[ist] ):
#		for il1,il2 in itertools.product( range(info.Nl[it1]), range(info.Nl[it2]) ):	
#			SI[it1,it2][il1][il2] = torch_complex.ComplexTensor(
#				torch.from_numpy(SI[it1,it2][il1][il2].real),
#				torch.from_numpy(SI[it1,it2][il1][il2].imag))
	return SI



def read_VI(info,V_info,ist,data):
	if V_info["same_band"]:
		""" VI[ib]	<psi|psi> """
		if V_info["init_from_file"]:
			VI = np.empty(info.Nb[ist],dtype=np.float64)
			for ib in range(info.Nb[ist]):
				VI.data[ib] = next(data)
		else:
			VI = np.ones(info.Nb[ist],dtype=np.float64)
	else:
		""" VI[ib1,ib2]	<psi|psi> """
		if V_info["init_from_file"]:
			VI = np.empty((info.Nb[ist],info.Nb[ist]),dtype=np.float64)
			for ib1,ib2 in itertools.product( range(info.Nb[ist]), range(info.Nb[ist]) ):
				VI[ib1,ib2] = next(data)
		else:
			VI = np.eye(info.Nb[ist],info.Nb[ist],dtype=np.float64)
	return torch.from_numpy(VI)
