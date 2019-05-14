import os
import sys
from read_info import read_info
from print_file import print_file_pw, print_file_opt
from dis import dis

def cal_pw():
	info = read_info("info.json")
	for distance in dis[info["input"]["element"]]:
		folder = f'{info["input"]["element"]}-{info["input"]["rcut"]}-{distance}'
		os.mkdir(folder)
		os.chdir(folder)
		print_file_pw(info,distance)
		if utils.sub=="qsub":
			os.system("qsub sub.sh")
		elif utils.sub=="bsub":
			os.system("bsub < sub.sh")
		else:
			raise KeyError("utils.sub = ",utils.sub)		
		os.chdir("../")
		
def cal_opt():
	info = read_info("info.json")
	os.mkdir("opt_orb")
	os.chdir("opt_orb")
	print_file_opt(info,dis)
	os.system("qsub sub.sh")
	if utils.sub=="qsub":
		os.system("qsub sub.sh")
	elif utils.sub=="bsub":
		os.system("bsub < sub.sh")
	else:
		raise KeyError("utils.sub = ",utils.sub)
	os.chdir("../")
		
if __name__=="__main__":
	if sys.argv[1]=='1':
		cal_pw()
	elif sys.argv[1]=='2':
		cal_opt()