import IO.read_QSV
import IO.print_QSV
import IO.func_C
import IO.read_json
import IO.print_orbital
import opt_orbital
import orbital
import torch
import numpy as np
import time
import torch_optimizer 
import IO.cal_weight
import util

def main():
	seed = int(1000*time.time())%(2**32)
	np.random.seed(seed)
	print("seed:",seed)
	time_start = time.time()

	file_list, info_true, weight_info, C_init_info, V_info = IO.read_json.read_json("input.json")

	weight = IO.cal_weight.cal_weight(weight_info, V_info["same_band"], file_list["origin"])

	QI,SI,VI,info = IO.read_QSV.read_file(info_true,file_list["origin"],V_info)
	print(info, flush=True)
	if "linear" in file_list.keys():
		QI_linear, SI_linear, VI_linear, info_linear = list(zip(*( IO.read_QSV.read_file(info_true,file,V_info) for file in file_list["linear"] )))

	if C_init_info["init_from_file"]:
		C, C_read_index = IO.func_C.read_C_init( C_init_info["C_init_file"], info )
	else:
		C = IO.func_C.random_C_init(info)
	E = orbital.set_E(info,info.Rcut)
	orbital.normalize( orbital.generate_orbital(info,C,E,info.Rcut,info.dr), info.dr,C,flag_norm_C=True)

	opt_orb = opt_orbital.Opt_Orbital()

	#opt = torch.optim.Adam(sum( ([c.real,c.imag] for c in sum(C,[])), []), lr=info.lr, eps=1e-8)
	#opt = torch.optim.Adam( sum(C.values(),[]), lr=info.lr, eps=1e-20, weight_decay=info.weight_decay)
	#opt = radam.RAdam( sum(C.values(),[]), lr=info.lr, eps=1e-20 )
	opt = torch_optimizer.SWATS( sum(C.values(),[]), lr=info.lr, eps=1e-20 )


	with open("Spillage.dat","w") as S_file:

		loss_old = np.inf
		for istep in range(10000):

			Q = opt_orb.change_index_Q(opt_orb.cal_Q(QI,C,info),info)
			S = opt_orb.change_index_S(opt_orb.cal_S(SI,C,info),info)
			V = opt_orb.cal_V(Q,S,info,V_info)

			if "linear" in file_list.keys():
				V_linear = [None] * len(file_list["linear"])
				for i in range(len(file_list["linear"])):
					Q_linear = opt_orb.change_index_Q(opt_orb.cal_Q(QI_linear[i],C,info),info)
					S_linear = opt_orb.change_index_S(opt_orb.cal_S(SI_linear[i],C,info),info)
					V_linear[i] = opt_orb.cal_V_linear(Q,S,Q_linear,S_linear,V,info,V_info)

			def cal_Spillage(V_delta):
				Spillage = torch.Tensor([0])
				for ist, Vi_delta in enumerate(V_delta):
					Spillage += (Vi_delta * weight[ist]).sum()
				return Spillage

			def cal_delta(VI, V):
				return ( ((VIi-Vi)/util.update0(VIi)).abs() for VIi,Vi in zip(VI,V) )		# abs or **2?
			
			Spillage = 2*cal_Spillage(cal_delta(VI,V))
			if "linear" in file_list.keys():
				for i in range(len(file_list["linear"])):
					Spillage += cal_Spillage(cal_delta(VI_linear[i],V_linear[i]))

			if info.cal_T:
				T = opt_orb.cal_T(C,E)
				if not "TSrate" in vars():	TSrate = torch.abs(0.002*Spillage/T).data[0]
				Loss = Spillage + TSrate*T
			else:
				Loss = Spillage

			if info.cal_T:
				print_content = [istep, Spillage.item(), T.item(), Loss.item()]
			else:
				print_content = [istep, Spillage.item()]
			print(*print_content, sep="\t", file=S_file, flush=True)
			if not istep%100:
				print(*print_content, sep="\t", flush=True)

			if Loss.item() < loss_old:
				loss_old = Loss.item()
				C_old = IO.func_C.copy_C(C,info)
				flag_finish = 0
			else:
				flag_finish += 1
				if flag_finish > 50:
					break

			opt.zero_grad()
			Loss.backward()		
			if C_init_info["init_from_file"] and not C_init_info["opt_C_read"]:
				for it,il,iu in C_read_index:
					C[it][il].grad[:,iu] = 0
			opt.step()
	#		orbital.normalize( orbital.generate_orbital(info,C,E,info.Rcut,info.dr), info.dr,C,flag_norm_C=True)

		orb = orbital.generate_orbital(info,C_old,E,info.Rcut,info.dr)
		if info.cal_smooth:
			orbital.smooth_orbital(orb,info.Rcut,info.dr,0.1)
		orbital.orth(orb,info.dr)
		IO.print_orbital.print_orbital(orb,info)
		IO.print_orbital.plot_orbital(orb,info.Rcut,info.dr)

		IO.func_C.write_C("C.dat",info,C_old)

		print("time:\t",time.time()-time_start)


if __name__=="__main__":
	import sys
	np.set_printoptions(threshold=sys.maxsize, linewidth=10000)
	
	main()