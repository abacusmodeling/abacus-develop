#include "Interpolation_3D.h"
#include "../src_pw/tools.h"
#include "random_def.h"
IPL_3D :: IPL_3D(): psi(NULL), result(NULL), pos(NULL)
{
	cout << " 3D_Interpolation " << endl;

	this->test = 3;
	this->flag = 0;

	if (flag)
	{
		if (test > 0)
		{
			cout << "Prepare Interpolation Infomation and write them to the hard disk!" << endl;
		}

		write();
	}
	else
	{
		if (test > 0)
		{
			cout << "Read Interpolation Information from hard disk and ready for Interpolate!" << endl;
		}

		read();
	}
}

IPL_3D :: ~IPL_3D()
{
//	delete[] psi;
//	delete[] pos;
};

void IPL_3D :: write()
{
	this -> init();

	this -> Cal_Coord();

	this -> QSort_Indata(psi, pos, 0, psi_num - 1);

	if (test > 1)
	{
		cout << "Correctness of Sort Function" << endl;

		for (int i = 0; i < psi_num; i++)
		{
			Vector3<int> temp = Get_Index(pos[i], v_min, cell_length);
			cout << "pos =  " << pos[i].x << setw(5) << pos[i].y << setw(5) << pos[i].z << setw(5);
			cout << "index =   (" << temp.x << ", " << temp.y << ",  " << temp.z << ")" << endl;
			int key = temp.x * length.y * length.z + temp.y * length.z + temp.z;
			cout << " key = " << key << "   i =  " << i << endl;
		}
	}

	ofstream outf("IPL_3D.dat");

	if (outf)
	{
		outf << IPL_nbnd << endl;
		outf << psi_num << endl;
		outf << old_lat0 << endl;
		outf << setw(12) << v_min.x
		<< setw(12) << v_min.y
		<< setw(12) << v_min.z << endl;
		outf << setw(12) << v_max.x
		<< setw(12) << v_max.y
		<< setw(12) << v_max.z << endl;
		outf << setw(12) << length.x
		<< setw(12) << length.y
		<< setw(12) << length.z << endl;
		outf << setw(12) << cell_length.x
		<< setw(12) << cell_length.y
		<< setw(12) << cell_length.z << endl;

		for (int i = 0; i < IPL_nbnd; i++)
		{
			for (int j = 0; j < psi_num;j++)
			{
				outf << setw(12) << psi[i][j].real() << setw(12) << psi[i][j].imag() << endl;
			}
		}

		for (int i = 0; i < psi_num; i++)
		{
			outf << setw(12) << pos[i].x
			<< setw(12) << pos[i].y
			<< setw(12) << pos[i].z << endl;
		}

		outf << setw(12) << old_latvec.e11 << setw(12) << old_latvec.e12 << setw(12) << old_latvec.e13 << endl

		<< setw(12) << old_latvec.e21 << setw(12) << old_latvec.e22 << setw(12) << old_latvec.e23 << endl
		<< setw(12) << old_latvec.e31 << setw(12) << old_latvec.e32 << setw(12) << old_latvec.e33 << endl;
	}
	else
	{
		cout << "Cannot open file IPL_3D!!" << endl;
		exit(0);
	}

	outf.close();

	if (test > 0)
		cout << "Finish the write process!" << endl;
}

void IPL_3D :: read()
{
	ifstream input("IPL_3D.dat");

	if (!input)
	{
		cout << "IPL_3D has not been saved yet!!" << endl;
		exit(0);
	}

	input >> IPL_nbnd;

	input >> this->psi_num;
	input >> this->old_lat0;
	input >> this->v_min.x >> this->v_min.y >> this->v_min.z;
	input >> this->v_max.x >> this->v_max.y >> this->v_max.z;
	input >> this->length.x >> this->length.y >> this->length.z;
	input >> this->cell_length.x >> this->cell_length.y >> this->cell_length.z;

	this->psi = new complex<double>* [IPL_nbnd]();

	for (int i = 0; i < IPL_nbnd; i++)
		this->psi[i] = new complex<double> [psi_num];

	for (int i = 0; i < nbnd; i++)
	{
		for (int j = 0; j < psi_num; j++)
		{
			double real, img;
			input >> real >> img;
			psi[i][j] = complex<double> (real, img);
		}
	}


	this->pos = new Vector3<double> [psi_num];

	for (int i = 0; i < psi_num; i++)
	{
		input >> this->pos[i].x >> this->pos[i].y >> this->pos[i].z;
	}

	input >> this->old_latvec.e11 >> this->old_latvec.e12 >> this->old_latvec.e13

	>> this->old_latvec.e21 >> this->old_latvec.e22 >> this->old_latvec.e23
	>> this->old_latvec.e31 >> this->old_latvec.e32 >> this->old_latvec.e33;
	input.close();

	if (test > 0)
	{
		cout << "Finish Reading!" << endl;
	}

//	this->Cal_Result();
}

void IPL_3D::init()
{
	if (test > 0)
	{
		cout << "\n init " << endl;

		cout << "\n kv.nks = " << kv.nks << endl;
		cout << "\n st.kgxyz = " << st.kgxyz << endl;
		cout << "\n st.qtot = " << st.qtot << endl;
	}

	this->length.x = st.kg[0];

	this->length.y = st.kg[1];
	this->length.z = st.kg[2];
	this->psi_num = st.kgxyz;
	this->old_latvec = pw.latvec;
	this->old_lat0 = pw.lat0;

	if (nbnd > pw.natomwfc)
	{
		IPL_nbnd = nbnd;
	}
	else IPL_nbnd = pw.natomwfc;

	this->psi = new complex<double>* [IPL_nbnd];

	for (int i = 0; i < IPL_nbnd; i++)
		this->psi[i] = new complex<double> [psi_num];

	for (int iw = 0;iw < IPL_nbnd; iw++)
	{
		int iq = 0;

		for (int ik = 0;ik < kv.nks; ik++)
		{
			for (int ig = 0;ig < kv.ngk[ik];ig++)
			{
				const int num = st.kg2fftw[iq];
				psi[iw][num] = wf.evc[ik](iw, ig);
				/*		if(test > 1)
						{
							cout<<"psi["<<num<<"] =   "<<psi[i][num].real()<<",  "<<psi[num][i].imag()<<endl;
						}
				*/		iq ++;
			}
		}
	}

	for (int i = 0; i < IPL_nbnd; i++)
	{
		//	fftwan.FFT3D(psi[i],1);
	}

	return;
}

void IPL_3D::Cal_Coord()
{
	if (test > 0)
		cout << "\n Cal_Coord : " << endl;

	this->pos = new Vector3<double> [psi_num];

	this->v_min = Vector3<double>(0, 0, 0);

	this->v_max = Vector3<double>(0, 0, 0);

	int x = st.kg[0];

	int y = st.kg[1];

	int z = st.kg[2];

	this->length = Vector3<int>(x, y, z);

	int xy = x * y;

	int xz = x * z;

	int yz = y * z;

	int x2 = x / 2;

	int y2 = y / 2;

	int z2 = z / 2;

	for (int i = 0;i < st.kgxyz; i++)
	{
		Vector3<double> temp;
		temp.x = i / xy;
		temp.y = i / x - temp.x * z;
		temp.z = i - temp.y * y - temp.x * yz;

		if (test > 1)
		{
			cout << "\ntemp  = ( " << temp.x
			     << ", " << temp.y << ", " << temp.z
			     << ")" << setw(5);
		}

		const int num = static_cast<int>(temp.x * yz + temp .y * z + temp.z);

		if (temp.x > x2) temp.x -= x;

		if (temp.y > y2) temp.y -= y;

		if (temp.z > z2) temp.z -= z;

		pos[num] = temp;

		if (test > 1)
		{
			cout << "pos[ " << num << "] = ( " << pos[num].x
			     << ", " << pos[num].y << ", " << pos[num].z
			     << ")" << endl;
		}

		if (pos[num].x < v_min.x) v_min.x = pos[num].x;

		if (pos[num].y < v_min.y) v_min.y = pos[num].y;

		if (pos[num].z < v_min.z) v_min.z = pos[num].z;

		if (pos[num].x > v_max.x) v_max.x = pos[num].x;

		if (pos[num].y > v_max.y) v_max.y = pos[num].y;

		if (pos[num].x > v_max.z) v_max.z = pos[num].z;

	}

	cell_length.x = (v_max.x - v_min.x) / (length.x - 1);

	cell_length.y = (v_max.y - v_min.y) / (length.y - 1);
	cell_length.z = (v_max.z - v_min.z) / (length.z - 1);

	if (test > 0)
	{
		cout << "\nv_min = ( " << v_min.x << " , " << v_min.y << " , " << v_min.z << " )\n"
		     << "v_max = ( " << v_max.x << " , " << v_max.y << " , " << v_max.z << " )\n"
		     << "\nlength.x =  " << length.x
		     << "\nlength.y =  " << length.y
		     << "\nlength.z =  " << length.z
		     << "\ncell_length.x =     " << cell_length.x
		     << "\ncell_length.y =     " << cell_length.y
		     << "\ncell_length.z =     " << cell_length.z << endl;
	}
}

Vector3<int> IPL_3D :: Get_Index(Vector3<double> pos, Vector3<double> Origin, Vector3<double> length)
{
	Vector3<int> Index;
	Index.x = static_cast<int>(std::floor((pos.x - Origin.x) / length.x));
	Index.y = static_cast<int>(std::floor((pos.y - Origin.y) / length.y));
	Index.z = static_cast<int>(std::floor((pos.z - Origin.z) / length.z));
	return Index;
}

int IPL_3D :: Sort_Partition(complex <double> **val, Vector3<double> *pos, int low, int high)
{
	if (test > 1)
		cout << "Begin Sort_Partition." << endl;

	Vector3<double> temp_pos = pos[low];

	complex <double>* temp_psi ;

	temp_psi = new complex<double> [IPL_nbnd];

	for (int i = 0;i < IPL_nbnd; i++)
	{
		temp_psi[i] = val[i][low];
	}

	Vector3<int> low_index = Get_Index(pos[low], v_min, cell_length);

	Vector3<int> high_index = Get_Index(pos[high], v_min, cell_length);

	int low_key = low_index.x * length.y * length.z + low_index.y * length.z + low_index.z;
	int pivotkey = low_key;
	int high_key = high_index.x * length.y * length.z + high_index.y * length.z + high_index.z;
	int times = 0;

	while (low < high)
	{
		if (test > 2)
		{
			cout << "Loop times:  " << times << endl;
			cout << "Low = " << low << endl;
			cout << "high = " << high << endl;
			cout << "low_key =  " << low_key << endl;
			cout << "high_key =  " << high_key << endl;
			cout << "pivotkey = " << pivotkey << endl;
			cout << "pos[low] = (" << pos[low].x
			     << " ," << pos[low].y
			     << " ," << pos[low].z << ")" << endl;
			cout << "pos[high] = (" << pos[high].x
			     << " ," << pos[high].y
			     << " ," << pos[high].z << ")" << endl;
		}

		times++;

		while (low < high && high_key >= pivotkey)
		{
			--high;
			high_index = Get_Index(pos[high], v_min, cell_length);
			high_key = high_index.x * length.y * length.z + high_index.y * length.z + high_index.z;
		}

		pos[low] = pos[high];

		for (int i = 0; i < IPL_nbnd; i++)
		{
			val[i][low] = val[i][high];
		}

		if (test > 2)
		{
			cout << "pos[low] = (" << pos[low].x
			     << " ," << pos[low].y
			     << " ," << pos[low].z << ")" << endl;
			cout << "low_key =  " << low_key << endl;
		}

		low_index = Get_Index(pos[low], v_min, cell_length);

		low_key = low_index.x * length.y * length.z + low_index.y * length.z + low_index.z;

		while (low < high && low_key <= pivotkey)
		{
			++low;
			low_index = Get_Index(pos[low], v_min, cell_length);
			low_key = low_index.x * length.y * length.z + low_index.y * length.z + low_index.z;
		}

		pos[high] = pos[low];

		for (int i = 0; i < IPL_nbnd; i++)
		{
			val[i][high] = val[i][low];
		}

		high_index = Get_Index(pos[high], v_min, cell_length);

		high_key = high_index.x * length.y * length.z + high_index.y * length.z + high_index.z;
	}

	if (test > 1)
		cout << "Loop times:  " << times << endl;

	pos[low] = temp_pos;

	for (int i = 0; i < IPL_nbnd; i++)
	{
		val[i][low] = temp_psi[i];
	}

	return low;
}

int IPL_3D::count = 0;
void IPL_3D::QSort_Indata(complex <double> **val, Vector3<double> *pos, int low, int high)
{
	if (test > 1)
	{
		cout << "\n<==========================Ready for the " << count
		     << " times Qsort=============================" << endl;

		if (test > 2)
		{
			for (int i = 0; i < psi_num; i++)
			{
				Vector3<int> temp = Get_Index(pos[i], v_min, cell_length);
				cout << "pos =  (" << pos[i].x << setw(5) << pos[i].y << setw(5) << pos[i].z << ")  " << setw(5);
				cout << "index =   (" << temp.x << ", " << temp.y << ",  " << temp.z << ")" << endl;
				int key = temp.x * length.y * length.z + temp.y * length.z + temp.z;
				cout << " key = " << key << "   i =  " << i << endl;
			}
		}

		cout << "\nlength of the Sqlist" << high - low + 1 << endl;
	}

	if (low < high)
	{
		if (test > 1)
			cout << "Iterate Time =              " << count << endl;

		count++;

		int pivotloc = this->Sort_Partition(val, pos, low, high);

		if (test > 1)
			cout << "pivotloc = " << pivotloc << endl;

//		if(count == 1)exit(0);
		QSort_Indata(val, pos, low, pivotloc - 1);

		QSort_Indata(val, pos, pivotloc + 1, high);
	}
}

complex<double> IPL_3D::Interpolation_3D(complex<double> **psi, const double &x, const double &y, const double &z, const int seq)
{
//	transform the coordinate
	Vector3<double> cor = Vector3<double>(x, y, z);
	cor = cor * this->old_latvec;
	cor = (old_lat0 / pw.lat0) * cor;
	cor = cor * pw.GT;

	if (test > 1)
	{
		cout << "band number =  " << seq << endl;
		cout << "Direct Coordinate in old Coordinates:  ("
		     << cor.x << ", " << cor.y << ", " << cor.z << ")" << endl;
	}

	Vector3<int> index = Get_Index(cor, v_min, cell_length);

	if (test > 1)
	{
		cout << "index in the old system:  (" << index.x << ", " << index.y << ", " << index.z << ")" << endl;
		cout << length.x << setw(5) << length.y << setw(5) << length.z << endl;
	}

	if (index.x > length.x || index.y > length.y || index.z > length.z)
	{
		cout << "Input Coordinates of is beyond the bound of the box!!" << endl;
		exit(0);
	}

	int pos = index.x * length.y * length.z + index.y * length.z + index.z;

	complex <double> V_000 = psi[seq][pos];

	pos = index.x * length.y * length.z + index.y * length.z + index.z + 1;
	complex <double> V_001 = psi[seq][pos];

	pos = index.x * length.y * length.z + (index.y + 1) * length.z + index.z;
	complex <double> V_010 = psi[seq][pos];

	pos = index.x * length.y * length.z + (index.y + 1) * length.z + index.z + 1;
	complex <double> V_011 = psi[seq][pos];

	pos = (index.x + 1) * length.y * length.z + index.y * length.z + index.z;
	complex <double> V_100 = psi[seq][pos];

	pos = (index.x + 1) * length.y * length.z + index.y * length.z + index.z + 1;
	complex <double> V_101 = psi[seq][pos];

	pos = (index.x + 1) * length.y * length.z + (index.y + 1) * length.z + index.z;
	complex <double> V_110 = psi[seq][pos];

	pos = (index.x + 1) * length.y * length.z + (index.y + 1) * length.z + index.z + 1;
	complex <double> V_111 = psi[seq][pos];

	complex <double> V_IPL;
	V_IPL = (1 - cor.x) * (1 - cor.y) * (1 - cor.z) * V_000 +
	        cor.x * (1 - cor.y) * (1 - cor.z) * V_100 +
	        (1 - cor.x) * cor.y * (1 - cor.z) * V_010 +
	        (1 - cor.x) * (1 - cor.y) * cor.z * V_001 +
	        cor.x * (1 - cor.y) * cor.z * V_101 +
	        (1 - cor.x) * cor.y * cor.z * V_011 +
	        cor.x * cor.y * (1 - cor.z) * V_110 +
	        cor.x * cor.y * cor.z * V_111 ;

	if (test > 1)
	{
		cout << "Interpolation Result:   (" << V_IPL.real() << ", " <<
		     V_IPL.imag() << ")" << endl;
	}

	return V_IPL;
}

void IPL_3D :: Cal_Result()
{
	if (test > 0)
		cout << "Cai_Result" << endl;

	int x = st.kg[0];

	int y = st.kg[1];

	int z = st.kg[2];

	int xy = x * y;

	int xz = x * z;

	int yz = y * z;

	int x2 = x / 2;

	int y2 = y / 2;

	int z2 = z / 2;

	this->result = new complex<double>* [IPL_nbnd];

	for (int i = 0; i < IPL_nbnd; i++)
	{
		this->result[i] = new complex<double> [st.kgxyz]();
	}

	for (int i = 0; i < IPL_nbnd; i++)
	{
		for (int j = 0;j < st.kgxyz; j++)
		{
			Vector3<double> temp;
			temp.x = j / xy;
			temp.y = j / x - temp.x * z;
			temp.z = j - temp.y * y - temp.x * yz;
			const int num = static_cast<int>(temp.x * yz + temp .y * z + temp.z);

			if (temp.x > x2) temp.x -= x;

			if (temp.y > y2) temp.y -= y;

			if (temp.z > z2) temp.z -= z;


			if (test > 1)
			{
				cout << "Direct Coordinate in new Coordinates:  ("
				     << temp.x << ", " << temp.y << ", " << temp.z << ")" << endl;
			}

			result[i][num] = this->Interpolation_3D(psi, temp.x, temp.y, temp.z, i);

		}
	}

	return;
}
