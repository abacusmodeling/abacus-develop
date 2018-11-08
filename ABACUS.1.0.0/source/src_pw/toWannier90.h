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


class toWannier90
{
public:
	const int k_supercell = 5;                                                              // default the k-space supercell
	const int k_cells = (2 * k_supercell + 1)*(2 * k_supercell + 1)*(2 * k_supercell + 1);  // the primitive cell number in k-space supercell
	const int k_shells = 12;                                                                // default the shell numbers
	const double large_number = 99999999.0;
	const double small_number = 0.000001;
	int num_kpts;                                                                           // k点的数目
	Matrix3 recip_lattice;
	vector<Vector3<double>> lmn;                                                            //每个k点原胞序号
	vector<double> dist_shell;                                                              //每一层shell的近邻k点距离
	vector<int> multi;                                                                      //每一层shell的近邻k点数目
	int num_shell_real;                                                                     //真正的满足B1条件的shell数目，最终结果。(注从1开始计数)
	int *shell_list_real;                                                                   //1到12层shell中不平行不等价的shell标签，长度为num_shell_real
	double *bweight;                                                                        //每个shell的bweight，长度为num_shell_real
	vector<vector<int>> nnlist;                                                             //每个k点的近邻k点序号
	vector<vector<Vector3<double>>> nncell;                                                 //每个k点的近邻k点所在的原胞序号
	int nntot = 0;                                                                          //每个k点的近邻k点数目   
	int num_wannier;																		//想要计算wannier函数的个数
	int *L;																					//试探轨道的角量子数指定,长度为num_wannier
	int *m;																					//试探轨道的磁量子数指定,长度为num_wannier
	int *rvalue;																			//试探轨道的径向部分函数形式,只有三种形式,长度为num_wannier
	double *alfa;																			//试探轨道的径向部分函数中的调节参数,长度为num_wannier
	Vector3<double> *R_centre;																//试探轨道函数中心,长度为num_wannier,cartesian坐标






	toWannier90(int num_kpts,Matrix3 recip_lattice);
	~toWannier90();

	void kmesh_supercell_sort(); //按照与原点的距离从小到大排序lmn
	void get_nnkpt_first();      //计算获得12层shell的近邻k点的距离和个数
	void kmesh_get_bvectors(int multi, int reference_kpt, double dist_shell, vector<Vector3<double>>& bvector);  //获取指定shell层，指定参考k点的近邻k点的bvector
	void get_nnkpt_last(); //获取最终的shell数目和bweight
    void get_nnlistAndnncell();

	void init_wannier();
	void cal_Amn();
	void cal_Mmn();
	void produce_trial_in_pw(const int &ik, ComplexMatrix &trial_orbitals_k);
	void integral(const int meshr, const double *psir, const double *r, const double *rab, const int &l, double* table);
	void writeUNK();
	void ToRealSpace(const int &ik, const int &ib, const ComplexMatrix *evc, complex<double> *psir, const Vector3<double> G);
	//void ToReciSpace(const complex<double> *psir, complex<double> *psik, const int ib);
	complex<double> unkdotb(const complex<double> *psir, const int ikb, const int bandindex);

};

#endif