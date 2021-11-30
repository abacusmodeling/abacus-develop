import addict

def change_info(info_old, weight_old):
	info_stru = [None] * info_old.Nst
	for ist in range(len(info_stru)):
		info_stru[ist] = addict.Dict()
	for ist,Na in enumerate(info_old.Na):
		info_stru[ist].Na = Na
	for ist,weight in enumerate(weight_old):
		info_stru[ist].weight = weight
		info_stru[ist].Nb = weight.shape[0]

	info_element = addict.Dict()
	for it,Nu in info_old.Nu.items():
		info_element[it].Nu = Nu
		info_element[it].Nl = len(Nu)
	for it,Rcut in info_old.Rcut.items():
		info_element[it].Rcut = Rcut
	for it,dr in info_old.dr.items():
		info_element[it].dr = dr
	for it,Ecut in info_old.Ecut.items():
		info_element[it].Ecut = Ecut
	for it,Ne in info_old.Ne.items():
		info_element[it].Ne = Ne

	info_opt = addict.Dict()
	info_opt.lr = info_old.lr
	info_opt.cal_T = info_old.cal_T
	info_opt.cal_smooth = info_old.cal_smooth

	return info_stru, info_element, info_opt