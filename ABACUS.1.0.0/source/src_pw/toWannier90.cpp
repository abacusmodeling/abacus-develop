#include "toWannier90.h"



toWannier90::toWannier90(int num_kpts, Matrix3 recip_lattice)
{
	this->num_kpts = num_kpts;
	this->recip_lattice = recip_lattice;

	//k空间supercell中每个小原胞的序号指定
	lmn.resize(k_cells);
	int counters = 0;
	for (int l = -k_supercell; l <= k_supercell; l++)
	{
		for (int m = -k_supercell; m <= k_supercell; m++)
		{
			for (int n = -k_supercell; n <= k_supercell; n++)
			{
				lmn[counters].set(l, m, n);
				counters++;
			}
		}
	}
	
	//初始化vector容器
	dist_shell.resize(k_shells);
	multi.resize(k_shells);


}


toWannier90::~toWannier90()
{
	delete[] shell_list_real;
	delete[] bweight;
}

void toWannier90::kmesh_supercell_sort()
{
	int size = lmn.size();
	vector<Vector3<double>> lmn_cpy;
	vector<pair<double, int>> dist;
	lmn_cpy.resize(size);
	dist.resize(size);
	for (int i = 0; i < size; i++)
	{
		dist[i].first = (lmn[i] * recip_lattice).norm2();
		dist[i].second = i;
	}
	sort(dist.begin(), dist.end(), [](pair<double, int> a, pair<double, int> b) -> bool {return a.first < b.first; });

	for (int i = 0; i < size; i++)
	{
		lmn_cpy[i] = lmn[dist[i].second];
	}

	for (int i = 0; i < size; i++)
	{
		lmn[i] = lmn_cpy[i];
	}
}

void toWannier90::get_nnkpt_first()
{
	int size = lmn.size();
	double dist = 0;
	double dist_0 = 0;
	double dist_1 = large_number;
	int counter = 0;

	

	// 计算获得12层shell的近邻k点的距离和个数，(以第一个k点为参考点)
	for (int shell = 0; shell < k_shells; shell++)
	{

		for (int nk = 0; nk < num_kpts; nk++)
		{
			for (int loop = 0; loop < size; loop++)
			{
				dist = (kv.kvec_c[nk] + lmn[loop] * recip_lattice - kv.kvec_c[0]).norm();

				if ((dist > small_number) && (dist > dist_0 + small_number))
				{
					if (dist < (dist_1 - small_number))
					{
						dist_1 = dist;
						counter = 0;
					}

					if ((dist >(dist_1 - small_number)) && (dist < (dist_1 + small_number))) counter++;

				}
			}
		}

		dist_shell[shell] = dist_1;
		multi[shell] = counter;
		dist_0 = dist_1;
		dist_1 = large_number;


	}

}

void toWannier90::kmesh_get_bvectors(int multi, int reference_kpt, double distshell, vector<Vector3<double>>& bvector)
{
	bvector.resize(multi);
	double dist = 0;
	int counter = -1;
	const double pi2 = 2 * 3.141592654;
	bool breakloop = true;
	for (int loop = 0; loop < k_cells && breakloop; loop++)
	{
		for (int nk = 0; nk < num_kpts && breakloop; nk++)
		{
			dist = (kv.kvec_c[nk] + lmn[loop] * recip_lattice - kv.kvec_c[reference_kpt]).norm();
			if ((dist > distshell * (1 - small_number)) && (dist < distshell * (1 + small_number)))
			{
				counter++;
				bvector[counter] = ( kv.kvec_c[nk] + lmn[loop] * recip_lattice - kv.kvec_c[reference_kpt] ) * pi2;
			}

			if ((counter + 1) == multi) breakloop = false;

		}



	}

}

void toWannier90::get_nnkpt_last()
{

	
	int num_shell = -1;
	vector<int> shell_list(k_shells);  // 记录每一层shell是否有效，无效值为0，否则为1
	double delta = 0.0;
	vector<vector<Vector3<double>>> bvector;
	bvector.resize(k_shells);
	matrix A(6, 1);  // 需要SVD分解的矩阵A
	matrix U(6, 6);	 // SVD分解中的U矩阵
	matrix VT(1, 1); // SVD分解中的V矩阵的转置
	matrix result(6, 1);
	result(0, 0) = 1.0; result(1, 0) = 1.0; result(2, 0) = 1.0; result(3, 0) = 0.0; result(4, 0) = 0.0; result(5, 0) = 0.0;
	bool B1_correct = false;
	

	
	for (int shell = 0; shell < k_shells; shell++)
	{
		kmesh_get_bvectors(multi[shell], 0, dist_shell[shell], bvector[shell]);
		bool parallel = false;
		if (shell > 0)
		{
			// 判断当前shell的bvector是否与前面所有shell平行

			for (int loop_shell = 0; loop_shell < shell; loop_shell++)
			{
				for (int loop_currentShell = 0; loop_currentShell < multi[shell]; loop_currentShell++)
				{
					for (int loop_beforeShell = 0; loop_beforeShell < multi[loop_shell]; loop_beforeShell++)
					{
						delta = (bvector[shell][loop_currentShell] * bvector[loop_shell][loop_beforeShell]) / (bvector[shell][loop_currentShell].norm()*bvector[loop_shell][loop_beforeShell].norm());
						if (abs(abs(delta) - 1.0) < 1.0e-6) parallel = true;
					}
				}
			}

		}
		
		

		if (!parallel)
		{
			num_shell++;
			shell_list[num_shell] = shell;
		}
		else
		{
			cout << " the shell " << shell << " is parallel !!! " << endl; 
			continue;
		}

		A.create(6, num_shell + 1);
		U.create(6, 6);
		VT.create(num_shell + 1, num_shell + 1);

		for (int loop_shell = 0; loop_shell <= num_shell; loop_shell++)
		{
			for (int loop = 0; loop < multi[shell_list[loop_shell]]; loop++)
			{
				A(0, loop_shell) += bvector[shell_list[loop_shell]][loop].x * bvector[shell_list[loop_shell]][loop].x;
				A(1, loop_shell) += bvector[shell_list[loop_shell]][loop].y * bvector[shell_list[loop_shell]][loop].y;
				A(2, loop_shell) += bvector[shell_list[loop_shell]][loop].z * bvector[shell_list[loop_shell]][loop].z;
				A(3, loop_shell) += bvector[shell_list[loop_shell]][loop].x * bvector[shell_list[loop_shell]][loop].y;
				A(4, loop_shell) += bvector[shell_list[loop_shell]][loop].y * bvector[shell_list[loop_shell]][loop].z;
				A(5, loop_shell) += bvector[shell_list[loop_shell]][loop].z * bvector[shell_list[loop_shell]][loop].x;

			}
			

		}

		


		double *S = new double[min(6, num_shell + 1)];  // SVD分解中的S数组，为一系列非负本征值
		int lwork = 6 * 10;
		double *work = new double[lwork];
		int info = 0;
		LapackConnector::dgesvd('A', 'A', 6, num_shell+1, A, 6, S, U, 6, VT, num_shell+1, work, lwork, info);

		if (info < 0) cout << "dgesvd is incorrect" << endl;
		else if (info > 0) cout << "dgesvd did not converge" << endl;

	
		
		bool effect = false;
		for (int i = 0; i < min(6, num_shell + 1); i++)
		{
			if (abs(S[i]) < 1.0e-5) effect = true;
		}

		if (effect)
		{
			if (num_shell == 0)
			{
				cout << "SVD has found a very small singular value" << endl;
			}
			else
			{
				cout << "SVD found small singular value, Rejecting this shell and trying the next" << endl;
				num_shell--;
				delete[] S;
				delete[] work;
				continue;
			}
		}
		
		matrix S_mat(num_shell+1, 6);
		for (int i = 0; i < min(6, num_shell + 1); i++)
		{
			S_mat(i, i) = 1.0 / S[i];
		}
		matrix bweight_mat( num_shell + 1, 1);


		//计算bweight矩阵
		bweight_mat = transpose(VT)*(S_mat * (transpose(U) * result));
		

		//检查 B1 条件是否成立，如果成立说明之前的SVD分解是正确的
		B1_correct = true;
		// lapack运算后A矩阵被重置
		A.create(6, num_shell + 1);		
		for (int loop_shell = 0; loop_shell <= num_shell; loop_shell++)
		{
			for (int loop = 0; loop < multi[shell_list[loop_shell]]; loop++)
			{
				A(0, loop_shell) += bvector[shell_list[loop_shell]][loop].x * bvector[shell_list[loop_shell]][loop].x;
				A(1, loop_shell) += bvector[shell_list[loop_shell]][loop].y * bvector[shell_list[loop_shell]][loop].y;
				A(2, loop_shell) += bvector[shell_list[loop_shell]][loop].z * bvector[shell_list[loop_shell]][loop].z;
				A(3, loop_shell) += bvector[shell_list[loop_shell]][loop].x * bvector[shell_list[loop_shell]][loop].y;
				A(4, loop_shell) += bvector[shell_list[loop_shell]][loop].y * bvector[shell_list[loop_shell]][loop].z;
				A(5, loop_shell) += bvector[shell_list[loop_shell]][loop].z * bvector[shell_list[loop_shell]][loop].x;

			}
		}
		matrix check_b1(6, 1);
		check_b1 = A*bweight_mat;
		for (int i = 0; i < 6; i++)
		{
			
			if ( abs(check_b1(i, 0) - result(i, 0)) > small_number )  
			{
				B1_correct = false;
				break;
			}
		}

		if (!B1_correct)
		{
			cout << " the B1 condition is not satisfy and check your KPT or the SVD is wrong" << endl;
			delete[] S;
			delete[] work;
			continue;
		}
		else
		{
			num_shell_real = num_shell+1;
			shell_list_real = new int[num_shell_real];
			bweight = new double[num_shell_real];
			for (int i = 0; i < num_shell_real; i++)
			{
				shell_list_real[i] = shell_list[i];
				bweight[i] = bweight_mat(i, 0);

			}

			delete[] S;
			delete[] work;
			break;

		}


	}




}


void toWannier90::get_nnlistAndnncell()
{
	for(int loop_s = 0; loop_s < num_shell_real; loop_s++)
	{
		nntot = nntot + multi[shell_list_real[loop_s]];
	}
	nnlist.resize(num_kpts);
	nncell.resize(num_kpts);
	for (int ik = 0; ik < num_kpts; ik++)
	{
		nnlist[ik].resize(nntot);
		nncell[ik].resize(nntot);
		int nn = -1;
		for (int loop_s = 0; loop_s < num_shell_real; loop_s++)
		{
			double dist = 0;
			int counter = -1;
			bool breakloop = true;
			for (int loop = 0; loop < k_cells && breakloop; loop++)
			{
				for (int jk = 0; jk < num_kpts && breakloop; jk++)
				{
					dist = (kv.kvec_c[jk] + lmn[loop] * recip_lattice - kv.kvec_c[ik]).norm();
					if ((dist > dist_shell[shell_list_real[loop_s]] * (1 - small_number)) && (dist < dist_shell[shell_list_real[loop_s]] * (1 + small_number)))
					{
						counter++;
						nn++;
						nnlist[ik][nn] = jk;
						nncell[ik][nn] = lmn[loop];
				
					}

					if ((counter + 1) == multi[shell_list_real[loop_s]]) breakloop = false;

				}


			}
			
		}
		
	}
	
	
}

void toWannier90::init_wannier()
{
	get_nnkpt_first();
	get_nnkpt_last();
	get_nnlistAndnncell();

	
	// read *.nnkp file
	ifstream nnkp_read(INPUT.NNKP.c_str());
	
	if( SCAN_BEGIN(nnkp_read,"projections") )
	{
		READ_VALUE(nnkp_read, num_wannier);
		// test
		cout << "num_wannier = " << num_wannier << endl;
		// test
		if(num_wannier < 0)
		{
			WARNING_QUIT("init_wannier","wannier number is lower than 0");
		}
		
		R_centre = new Vector3<double>[num_wannier];
		L = new int[num_wannier];
		m = new int[num_wannier];
		rvalue = new int[num_wannier];
		Vector3<double>* z_axis = new Vector3<double>[num_wannier];
		Vector3<double>* x_axis = new Vector3<double>[num_wannier];
		alfa = new double[num_wannier];
		
		
		for(int count = 0; count < num_wannier; count++)
		{
			nnkp_read >> R_centre[count].x >> R_centre[count].y >> R_centre[count].z;
			nnkp_read >> L[count] >> m[count];
			READ_VALUE(nnkp_read,rvalue[count]);
			nnkp_read >> z_axis[count].x >> z_axis[count].y >> z_axis[count].z;
			nnkp_read >> x_axis[count].x >> x_axis[count].y >> x_axis[count].z;
			READ_VALUE(nnkp_read,alfa[count]);			
		}
		
	}	
	
	nnkp_read.close();
	
	for(int i = 0; i < num_wannier; i++)
	{
		R_centre[i] = R_centre[i] * ucell.latvec;
		m[i] = m[i] - 1; // ABACUS and wannier90 对磁角动量m的定义不一样，ABACUS是从0开始的，wannier90是从1开始的
	}
	
	

	cout << "num_wannier is " << num_wannier << endl;
	for(int i = 0; i < num_wannier; i++)
	{
		cout << "num_wannier" << endl;
		cout << L[i] << " " << m[i] << " " << rvalue[i] << " " << alfa[i] << endl;
	}

	
	if(MY_RANK==0)
	{
		cal_Amn();
		cal_Mmn();
		writeUNK();
	}
}


void toWannier90::writeUNK()
{

	
	complex<double> *porter = new complex<double>[pw.nrxx];
	
	for(int ik = 0; ik < num_kpts; ik++)
	{
		stringstream name;
		name << "UNK" << setw(5) << setfill('0') << ik+1 << ".1" ;
		ofstream unkfile(name.str());
		
		unkfile << setw(12) << pw.nx << setw(12) << pw.ny << setw(12) << pw.nz << setw(12) << ik+1 << setw(12) << NBANDS << endl;
		
		for(int ib = 0; ib < NBANDS; ib++)
		{
			//complex<double> *porter = UFFT.porter;
			//  u_k in real space
			ZEROS(porter, pw.nrxx);
			for (int ig = 0; ig < kv.ngk[ik]; ig++)
			{
				porter[pw.ig2fftw[wf.igk(ik, ig)]] = wf.evc[ik](ib, ig);
			}
			pw.FFT_wfc.FFT3D(porter, 1);
			
			for(int k=0; k<pw.nz; k++)
			{
				for(int j=0; j<pw.ny; j++)
				{
					for(int i=0; i<pw.nx; i++)
					{
						unkfile << "       " << setw(12) << setprecision(9) << setiosflags(ios::scientific) << porter[i*pw.ncy*pw.ncz + j*pw.ncz + k].real()
								<< "       " << setw(12) << setprecision(9) << porter[i*pw.ncy*pw.ncz + j*pw.ncz + k].imag() 
								//jingan test
								//<< "       " << setw(12) << setprecision(9) << setiosflags(ios::scientific) << abs(porter[i*pw.ncy*pw.ncz + j*pw.ncz + k])
						        << endl;
					}
				}
			}
			
			
		}

		
		unkfile.close();
		
	}
	
	delete[] porter;
}







void toWannier90::cal_Amn()
{
	// 第一步：建立实球谐函数lm在某个k点下的平面波基组下的表格（矩阵）	
	// 第二步：将试探轨道的径向部分向某个k点下平面波投影
	// 第三步：获取试探轨道在某个k点下平面波基组下的投影
	string fileaddress = global_out_dir + "seedname.amn";
	ofstream Amn_file( fileaddress.c_str() );
	const int pwNumberMax = wf.npwx;
	ComplexMatrix *trial_orbitals = new ComplexMatrix[num_kpts];
	for(int ik = 0; ik < num_kpts; ik++)
	{
		trial_orbitals[ik].create(num_wannier,pwNumberMax);
		produce_trial_in_pw(ik,trial_orbitals[ik]);
	}
	
	
	time_t  time_now = time(NULL);
	Amn_file << " Created on " << ctime(&time_now);
	Amn_file << setw(12) << NBANDS << setw(12) << num_kpts << setw(12) << num_wannier << endl;
	
	for(int ik = 0; ik < num_kpts; ik++)
	{
		for(int iw = 0; iw < num_wannier; iw++)
		{
			for(int ib = 0; ib < NBANDS; ib++)
			{
				complex<double> amn(0.0,0.0);
				for(int ig = 0; ig < pwNumberMax; ig++)
				{
					amn = amn + conj( wf.evc[ik](ib,ig) ) * trial_orbitals[ik](iw,ig);
				}
				
				Amn_file << "   " << setw(5) << setiosflags(ios::left) << ib+1 << setw(5) << iw+1 << setw(5) << ik+1 
						 << "   " << setw(18) << setprecision(13) << amn.real() 
						 << "   " << setw(18) << setprecision(13) << amn.imag()
						 //jingan test
						 //<< "   " << setw(18) << setprecision(13) << abs(amn)
						 << endl;
			}
		}
	}
	
	Amn_file.close();
	
}



void toWannier90::cal_Mmn()
{
	
	string fileaddress = global_out_dir + "seedname.mmn";
	ofstream mmn_file( fileaddress.c_str() );	
	
	time_t  time_now = time(NULL);
	mmn_file << " Created on " << ctime(&time_now);
	mmn_file << setw(12) << NBANDS << setw(12) << num_kpts << setw(12) << nntot << endl;
	
	for(int ik = 0; ik < num_kpts; ik++)
	{
		for(int ib = 0; ib < nntot; ib++)
		{
			int ikb = nnlist[ik][ib];             // ik+b : ik的近邻k点		
			Vector3<double> phase_G = nncell[ik][ib];
			int G_number = kv.ngk[ib];
			mmn_file << "   " << setw(5) << ik+1 << setw(5) << ikb+1 << setw(5) 
					 << int(phase_G.x) << setw(5) << int(phase_G.y) << setw(5) << int(phase_G.z) 
					 << endl;
			
			for(int n = 0; n < NBANDS; n++)
			{
				for(int m = 0; m < NBANDS; m++)
				{
					complex<double> mmn(0.0,0.0);
					
					complex<double> *unk_L_r = new complex<double>[pw.nrxx];
					//complex<double> *unk_L_k;
					ToRealSpace(ik,m,wf.evc,unk_L_r,phase_G);
					//ToReciSpace(unk_L_r,unk_L_k,ikb);					
					
					mmn = unkdotb(unk_L_r,ikb,n);
					
					mmn_file << "    " << setw(12) << setprecision(9) << setiosflags(ios::scientific) << mmn.real() 
							 << "    " << setw(12) << setprecision(9) << mmn.imag()
							 // jingan test
							 //<< "    " << setw(12) << setprecision(9) << abs(mmn)
							 << endl;
					
					//delete[] unk_L_k;
					delete[] unk_L_r;
					
				}
			}
		}
		
	}
	
	mmn_file.close();
	
}


void toWannier90::produce_trial_in_pw(const int &ik, ComplexMatrix &trial_orbitals_k)
{
	// 检查参数是否正确
	for(int i =0; i < num_wannier; i++)
	{
		if(L[i] < -5 || L[i] > 3) cout << "toWannier90::produce_trial_in_pw() your L angular momentum is wrong , please check !!! " << endl;
	
		if(L[i] >= 0) 
		{
			if(m[i] < 0 || m[i] > 2*L[i]) cout << "toWannier90::produce_trial_in_pw() your m momentum is wrong , please check !!! " << endl;
		}
		else
		{
			if(m[i] < 0 || m[i] > -L[i]) cout << "toWannier90::produce_trial_in_pw() your m momentum is wrong , please check !!! " << endl;
		
		}
	}
	
	const int npw = kv.ngk[ik];
	const int total_lm = 16;
	matrix ylm(total_lm,npw);               //所有类型的球谐函数
	matrix wannier_ylm(num_wannier,npw);    //要试探轨道的使用的球谐函数
	double bs2, bs3, bs6, bs12;
	bs2 = 1.0/sqrt(2.0);
	bs3 = 1.0/sqrt(3.0);
	bs6 = 1.0/sqrt(6.0);
	bs12 = 1.0/sqrt(12.0);
	
	Vector3<double> *gk = new Vector3<double>[npw];
	for(int ig = 0; ig < npw; ig++)
	{
		gk[ig] = wf.get_1qvec_cartesian(ik, ig);  // k+G矢量
	}
	
	Mathzone::Ylm_Real(total_lm, npw, gk, ylm);
	
	// 1.生成 wannier_ylm
	for(int i = 0; i < num_wannier; i++)
	{
		if(L[i] >= 0)
		{
			for(int ig = 0; ig < npw; ig++)
			{
				int index = L[i] * L[i] + m[i];
				if(index == 2 || index == 3 || index == 5 || index == 6 || index == 14 || index == 15)
				{
					wannier_ylm(i,ig) = -1 * ylm(index,ig);
				}
				else
				{
					wannier_ylm(i,ig) = ylm(index,ig);
				}
			}
		}
		else
		{
			if(L[i] == -1 && m[i] == 0)
			{
				for(int ig = 0; ig < npw; ig++)
				{
					wannier_ylm(i,ig) = bs2 * ylm(0,ig) + bs2 * -ylm(2,ig);
				}
			}
			else if(L[i] == -1 && m[i] == 1)
			{
				for(int ig = 0; ig < npw; ig++)
				{
					wannier_ylm(i,ig) = bs2 * ylm(0,ig) - bs2 * -ylm(2,ig);
				}				
			}
			else if(L[i] == -2 && m[i] == 0)
			{
				for(int ig = 0; ig < npw; ig++)
				{
					wannier_ylm(i,ig) = bs3 * ylm(0,ig) - bs6 * -ylm(2,ig) + bs2 * -ylm(3,ig);
				}				
			}
			else if(L[i] == -2 && m[i] == 1)
			{
				for(int ig = 0; ig < npw; ig++)
				{
					wannier_ylm(i,ig) = bs3 * ylm(0,ig) - bs6 * -ylm(2,ig) - bs2 * -ylm(3,ig);
				}				
			}			
			else if(L[i] == -2 && m[i] == 2)
			{
				for(int ig = 0; ig < npw; ig++)
				{
					wannier_ylm(i,ig) = bs3 * ylm(0,ig) + 2 * bs6 * -ylm(2,ig);
				}				
			}			
			else if(L[i] == -3 && m[i] == 0)
			{
				for(int ig = 0; ig < npw; ig++)
				{
					wannier_ylm(i,ig) = 0.5*(ylm(0,ig) + -ylm(2,ig) + -ylm(3,ig) + ylm(1,ig));
				}				
			}			
			else if(L[i] == -3 && m[i] == 1)
			{
				for(int ig = 0; ig < npw; ig++)
				{
					wannier_ylm(i,ig) = 0.5*(ylm(0,ig) + -ylm(2,ig) - -ylm(3,ig) - ylm(1,ig));
				}				
			}			
			else if(L[i] == -3 && m[i] == 2)
			{
				for(int ig = 0; ig < npw; ig++)
				{
					wannier_ylm(i,ig) = 0.5*(ylm(0,ig) - -ylm(2,ig) + -ylm(3,ig) - ylm(1,ig));
				}				
			}			
			else if(L[i] == -3 && m[i] == 3)
			{
				for(int ig = 0; ig < npw; ig++)
				{
					wannier_ylm(i,ig) = 0.5*(ylm(0,ig) - -ylm(2,ig) - -ylm(3,ig) + ylm(1,ig));
				}				
			}			
			else if(L[i] == -4 && m[i] == 0)
			{
				for(int ig = 0; ig < npw; ig++)
				{
					wannier_ylm(i,ig) = bs3 * ylm(0,ig) - bs6 * -ylm(2,ig) + bs2 * -ylm(3,ig);
				}				
			}			
			else if(L[i] == -4 && m[i] == 1)
			{
				for(int ig = 0; ig < npw; ig++)
				{
					wannier_ylm(i,ig) = bs3 * ylm(0,ig) - bs6 * -ylm(2,ig) - bs2 * -ylm(3,ig);
				}				
			}			
			else if(L[i] == -4 && m[i] == 2)
			{
				for(int ig = 0; ig < npw; ig++)
				{
					wannier_ylm(i,ig) = bs3 * ylm(0,ig) - 2 * bs6 * -ylm(2,ig);
				}				
			}			
			else if(L[i] == -4 && m[i] == 3)
			{
				for(int ig = 0; ig < npw; ig++)
				{
					wannier_ylm(i,ig) = bs2 * ylm(1,ig) + bs2 * ylm(4,ig);
				}				
			}			
			else if(L[i] == -4 && m[i] == 4)
			{
				for(int ig = 0; ig < npw; ig++)
				{
					wannier_ylm(i,ig) = -1.0 * bs2 * ylm(1,ig) + bs2 * ylm(4,ig);
				}				
			}			
			else if(L[i] == -5 && m[i] == 0)
			{
				for(int ig = 0; ig < npw; ig++)
				{
					wannier_ylm(i,ig) = bs6 * ylm(0,ig) - bs2 * -ylm(2,ig) - bs12 * ylm(4,ig) + 0.5 * ylm(7,ig);
				}				
			}			
			else if(L[i] == -5 && m[i] == 1)
			{
				for(int ig = 0; ig < npw; ig++)
				{
					wannier_ylm(i,ig) = bs6 * ylm(0,ig) + bs2 * -ylm(2,ig) - bs12 * ylm(4,ig) + 0.5 * ylm(7,ig);
				}				
			}			
			else if(L[i] == -5 && m[i] == 2)
			{
				for(int ig = 0; ig < npw; ig++)
				{
					wannier_ylm(i,ig) = bs6 * ylm(0,ig) - bs2 * -ylm(3,ig) - bs12 * ylm(4,ig) - 0.5 * ylm(7,ig);
				}				
			}			
			else if(L[i] == -5 && m[i] == 3)
			{
				for(int ig = 0; ig < npw; ig++)
				{
					wannier_ylm(i,ig) = bs6 * ylm(0,ig) + bs2 * -ylm(3,ig) - bs12 * ylm(4,ig) - 0.5 * ylm(7,ig);
				}				
			}			
			else if(L[i] == -5 && m[i] == 4)
			{
				for(int ig = 0; ig < npw; ig++)
				{
					wannier_ylm(i,ig) = bs6 * ylm(0,ig) - bs2 * ylm(1,ig) + bs3 * ylm(4,ig);
				}				
			}			
			else if(L[i] == -5 && m[i] == 5)
			{
				for(int ig = 0; ig < npw; ig++)
				{
					wannier_ylm(i,ig) = bs6 * ylm(0,ig) + bs2 * ylm(1,ig) + bs3 * ylm(4,ig);
				}				
			}			
			
			
		}	
		

	}
	

	// 2.生成径向轨道在某个k点平面波基组的投影
	const int mesh_r = 333; 		//描述径向函数所需要的格点数
	const double dx = 0.025; 		//固定间隔，用于生成非固定间隔的dr来提高精度,这个值很巧妙
	const double x_min = -6.0;  	// 用于生成dr和r的起始点
	matrix r(num_wannier,mesh_r);   //不同alfa的径向函数的r
	matrix dr(num_wannier,mesh_r);  //不同alfa的径向函数的每个r点的间隔
	matrix psi(num_wannier,mesh_r); //径向函数psi in 实空间
	matrix psir(num_wannier,mesh_r);// psi * r in 实空间
	matrix psik(num_wannier,npw);   //径向函数在某个k点下倒空间的投影
	
	// 生成r,dr
	for(int i = 0; i < num_wannier; i++)
	{
		double x = 0;
		for(int ir = 0; ir < mesh_r; ir++)
		{
			x = x_min + ir * dx;
			r(i,ir) = exp(x) / alfa[i];
			dr(i,ir) = dx * r(i,ir);
		}
		
	}
	
	// 生成psi
	for(int i = 0; i < num_wannier; i++)
	{
		double alfa32 = pow(alfa[i],3.0/2.0);
		double alfa_new = alfa[i];
		int wannier_index = i;
		
		if(rvalue[i] == 1)
		{
			for(int ir = 0; ir < mesh_r; ir++)
			{
				psi(wannier_index,ir) = 2.0 * alfa32 * exp( -alfa_new * r(wannier_index,ir) );
			}
		}
	
		if(rvalue[i] == 2)
		{
			for(int ir = 0; ir < mesh_r; ir++)
			{
				psi(wannier_index,ir) = 1.0/sqrt(8.0) * alfa32
										* (2.0 - alfa_new * r(wannier_index,ir))
										* exp( -alfa_new * r(wannier_index,ir) * 0.5 );
			}
		}

		if(rvalue[i] == 3)
		{
			for(int ir = 0; ir < mesh_r; ir++)
			{
				psi(wannier_index,ir) = sqrt(4.0/27.0) * alfa32
										* ( 1.0 - 2.0/3.0 * alfa_new * r(wannier_index,ir) + 2.0/27.0 * pow(alfa_new,2.0) * r(wannier_index,ir) * r(wannier_index,ir) )
										* exp( -alfa_new * r(wannier_index,ir) * 1.0/3.0 );
			}
		}
		
	}

	// 生成psir
	for(int i = 0; i < num_wannier; i++)
	{
		for(int ir = 0; ir < mesh_r; ir++)
		{
			psir(i,ir) = psi(i,ir) * r(i,ir);
		}
	}
	
	//计算径向函数在某个k点下倒空间的投影
	double *psir_tem = new double[mesh_r];
	double *r_tem = new double[mesh_r];
	double *dr_tem = new double[mesh_r];
	double *psik_tem = new double[NQX];    //径向函数在固定k空间的投影（临时使用的数组）
	ZEROS(psir_tem,mesh_r);
	ZEROS(r_tem,mesh_r);
	ZEROS(dr_tem,mesh_r);
	
	for(int wannier_index = 0; wannier_index < num_wannier; wannier_index++)
	{

		for(int ir = 0; ir < mesh_r; ir++)
		{
			psir_tem[ir] = psir(wannier_index,ir);
			r_tem[ir] = r(wannier_index,ir);
			dr_tem[ir] = dr(wannier_index,ir);
		}
		
		toWannier90::integral(mesh_r,psir_tem,r_tem,dr_tem,L[wannier_index],psik_tem);
		
		for(int ig = 0; ig < npw; ig++)
		{
			psik(wannier_index,ig) = Mathzone::Polynomial_Interpolation(psik_tem, NQX, DQ, gk[ig].norm() * ucell.tpiba);
		}
	
		ZEROS(psir_tem,mesh_r);
		ZEROS(r_tem,mesh_r);
		ZEROS(dr_tem,mesh_r);
		ZEROS(psik_tem,NQX);

	}
	
	
	// 3.计算与原点选择（即轨道中心）而产生的相位在平面波基组下

	ComplexMatrix sk(num_wannier,npw);
	for(int wannier_index = 0; wannier_index < num_wannier; wannier_index++)
	{ 
		for(int ig = 0; ig < npw; ig++)
		{
			const double arg = ( gk[ig] * R_centre[wannier_index] ) * TWO_PI;
			sk(wannier_index,ig) = complex <double> ( cos(arg),  -sin(arg) );
		}
	}
	
	
	// 4.计算最终试探轨道在某个k点下平面波基组的投影
	
	//ComplexMatrix trial_orbitals_k(num_wannier,npw);
	for(int wannier_index = 0; wannier_index < num_wannier; wannier_index++)
	{
		complex<double> lphase = pow(NEG_IMAG_UNIT, L[wannier_index]);
		//test by jingan
		//cout << "the lphase in wannier_index " << wannier_index << " value is " << lphase << endl;
		for(int ig = 0; ig < wf.npwx; ig++)
		{
			if(ig < npw)
			{
				trial_orbitals_k(wannier_index,ig) = lphase * sk(wannier_index,ig) * wannier_ylm(wannier_index,ig) * psik(wannier_index,ig);
			}
			else trial_orbitals_k(wannier_index,ig) = complex<double>(0.0,0.0);
		}
	}
	
	// 5.归一化
	for(int wannier_index = 0; wannier_index < num_wannier; wannier_index++)
	{
		complex<double> anorm(0.0,0.0);
		for(int ig = 0; ig < wf.npw; ig++)
		{
			anorm = anorm + conj(trial_orbitals_k(wannier_index,ig)) * trial_orbitals_k(wannier_index,ig);
		}
		
		for(int ig = 0; ig < wf.npw; ig++)
		{
			trial_orbitals_k(wannier_index,ig) = trial_orbitals_k(wannier_index,ig) / sqrt(anorm);
		}
	}
	
}


void toWannier90::integral(const int meshr, const double *psir, const double *r, const double *rab, const int &l, double* table)
{
	const double pref = FOUR_PI / sqrt(ucell.omega);
	
	double *inner_part = new double[meshr];
	for(int ir=0; ir<meshr; ir++)
	{
		inner_part[ir] = psir[ir] * psir[ir];
	}
	
	double unit = 0.0;
	Mathzone::Simpson_Integral(meshr, inner_part, rab, unit);
	delete[] inner_part;
	OUT(ofs_running,"normalize unit",unit);

	double *aux = new double[meshr];
	double *vchi = new double[meshr];
	for (int iq=0; iq<NQX; iq++)
	{
		const double q = DQ * iq;
		Mathzone::Spherical_Bessel(meshr, r, q, l, aux);
		for (int ir = 0;ir < meshr;ir++)
		{
			vchi[ir] = psir[ir] * aux[ir] * r[ir];
		}
		
		double vqint = 0.0;
		Mathzone::Simpson_Integral(meshr, vchi, rab, vqint);

		table[iq] =  vqint * pref;
	}
	delete[] aux;
	delete[] vchi;
	return;
}


void toWannier90::ToRealSpace(const int &ik, const int &ib, const ComplexMatrix *evc, complex<double> *psir, const Vector3<double> G)
{
	// (1) set value
	complex<double> *phase = UFFT.porter;
    ZEROS( psir, pw.nrxx );
	ZEROS( phase, pw.nrxx);

    for (int ig = 0; ig < kv.ngk[ik]; ig++)
    {
        psir[ pw.ig2fftw[ wf.igk(ik,ig) ] ] = evc[ik](ib, ig);
    }
	
	// get the phase value in realspace
	for (int ig = 0; ig < pw.ngmw; ig++)
	{
		if (pw.gdirect[ig] == G)
		{
			phase[ pw.ig2fftw[ig] ] = complex<double>(1.0,0.0);
			break;
		}
	}
	// (2) fft and get value
    pw.FFT_wfc.FFT3D(psir, 1);
	pw.FFT_wfc.FFT3D(phase, 1);
	

	
	for (int ir = 0; ir < pw.nrxx; ir++)
	{
		psir[ir] = psir[ir] * phase[ir];
	}
    return;
}

complex<double> toWannier90::unkdotb(const complex<double> *psir, const int ikb, const int bandindex)
{
	complex<double> result(0.0,0.0);
	int knumber = kv.ngk[ikb];
	complex<double> *porter = UFFT.porter;
	ZEROS( porter, pw.nrxx);
	for (int ir = 0; ir < pw.nrxx; ir++)
	{
		porter[ir] = psir[ir];
	}
	pw.FFT_wfc.FFT3D( porter, -1);
	
	
	for (int ig = 0; ig < knumber; ig++)
	{
		result = result + conj( porter[ pw.ig2fftw[wf.igk(ikb, ig)] ] ) * wf.evc[ikb](bandindex,ig);	
		
	}
	return result;
}


