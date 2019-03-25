#ifndef TOWannier90_H
#define TOWannier90_H

#include <iostream>
using namespace std;
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include "src_pw/tools.h"
#include "src_global/lapack_connector.h"
#include "src_pw/global.h"
#include "src_pw/wavefunc_in_pw.h"
#include "src_lcao/local_orbital_wfc.h"


class toWannier90
{
public:
	//const int k_supercell = 5;                                                              // default the k-space supercell
	//const int k_cells = (2 * k_supercell + 1)*(2 * k_supercell + 1)*(2 * k_supercell + 1);  // the primitive cell number in k-space supercell
	//const int k_shells = 12;                                                                // default the shell numbers
	//const double large_number = 99999999.0;
	//const double small_number = 0.000001;
	//vector<Vector3<double>> lmn;                                                            //每个k点原胞序号
	//vector<double> dist_shell;                                                              //每一层shell的近邻k点距离
	//vector<int> multi;                                                                      //每一层shell的近邻k点数目
	//int num_shell_real;                                                                     //真正的满足B1条件的shell数目，最终结果。(注从1开始计数)
	//int *shell_list_real;                                                                   //1到12层shell中不平行不等价的shell标签，长度为num_shell_real
	//double *bweight;                                                                        //每个shell的bweight，长度为num_shell_real
	
	int num_kpts;                                                                           // k点的数目
	int cal_num_kpts;                                                                       // 需要计算的k的数目，对于nspin=2时有用处
	Matrix3 recip_lattice;
	vector<vector<int>> nnlist;                                                             //每个k点的近邻k点序号
	vector<vector<Vector3<double>>> nncell;                                                 //每个k点的近邻k点所在的原胞序号
	int nntot = 0;                                                                          //每个k点的近邻k点数目   
	int num_wannier;																		//想要计算wannier函数的个数
	int *L;																					//试探轨道的角量子数指定,长度为num_wannier
	int *m;																					//试探轨道的磁量子数指定,长度为num_wannier
	int *rvalue;																			//试探轨道的径向部分函数形式,只有三种形式,长度为num_wannier
	double *alfa;																			//试探轨道的径向部分函数中的调节参数,长度为num_wannier
	Vector3<double> *R_centre;																//试探轨道函数中心,长度为num_wannier,cartesian坐标
	string wannier_file_name = "seedname";                                                  // .mmn,.amn文件名
	int num_exclude_bands = 0;																// 排除计算的能带数目，-1表示没有需要排除的能带
	int *exclude_bands;                                                                     // 排除能带的index
	bool *tag_cal_band;																		// 判断NBANDS能带那一条需要计算
	int num_bands;																		   	// wannier90 中的num_bands
	bool gamma_only_wannier = false;														// 只用gamma点生成wannier函数
	string wannier_spin = "up";                                                             // spin参数，有up,down两个参数
	int start_k_index = 0;                                                                  // 用于for循环寻找k的指标，spin=2时开始的index是不一样的

	
	// 以下是lcao基组下的wannier90所需参数
	realArray table_local;
	ComplexMatrix *unk_inLcao;                                                             // lcao基组下波函数的周期部分unk



	toWannier90(int num_kpts,Matrix3 recip_lattice);
	~toWannier90();

	//void kmesh_supercell_sort(); //按照与原点的距离从小到大排序lmn
	//void get_nnkpt_first();      //计算获得12层shell的近邻k点的距离和个数
	//void kmesh_get_bvectors(int multi, int reference_kpt, double dist_shell, vector<Vector3<double>>& bvector);  //获取指定shell层，指定参考k点的近邻k点的bvector
	//void get_nnkpt_last(); //获取最终的shell数目和bweight
    //void get_nnlistAndnncell();

	void init_wannier();
	void read_nnkp();
	void outEIG();
	void cal_Amn(const ComplexMatrix *wfc_pw);
	void cal_Mmn(const ComplexMatrix *wfc_pw);
	void produce_trial_in_pw(const int &ik, ComplexMatrix &trial_orbitals_k);
	void get_trial_orbitals_lm_k(const int wannier_index, const int orbital_L, const int orbital_m, matrix &ylm, 
										matrix &dr, matrix &r, matrix &psir, const int mesh_r, 
										Vector3<double> *gk, const int npw, ComplexMatrix &trial_orbitals_k);
	void integral(const int meshr, const double *psir, const double *r, const double *rab, const int &l, double* table);
	void writeUNK(const ComplexMatrix *wfc_pw);
	void ToRealSpace(const int &ik, const int &ib, const ComplexMatrix *evc, complex<double> *psir, const Vector3<double> G);
	complex<double> unkdotb(const complex<double> *psir, const int ikb, const int bandindex, const ComplexMatrix *wfc_pw);
	complex<double> gamma_only_cal(const int &ib_L, const int &ib_R, const ComplexMatrix *wfc_pw, const Vector3<double> G);
	
	// lcao部分
	void lcao2pw_basis(const int ik, ComplexMatrix &orbital_in_G);
	void getUnkFromLcao();
	void get_lcao_wfc_global_ik(complex<double> **ctot, complex<double> **cc);

};

#endif