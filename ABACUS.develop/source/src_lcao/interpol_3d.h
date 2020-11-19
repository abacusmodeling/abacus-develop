#ifndef INCLUDE_IPL_3D
#define INCLUDE_IPL_3D

#include "../src_pw/global.h"
#include "../src_pw/tools.h"

class IPL_3D
{
public:
	IPL_3D();
	~IPL_3D();

	complex<double> **psi;
	int psi_num;
	Vector3<double> *pos;
	complex<double> **result;
private:
	static int count ;//count iterater time in the Sort function
	Vector3<double> v_min;//origin of the box
	Vector3<double> v_max;//farthest points of the box
	Vector3<int> length;//number of cells on each dim
	Vector3<double> cell_length;// distance between two adjacent points in each dim
	Matrix3 old_latvec;
	double old_lat0;
	int IPL_nbnd;
	int test;//switch of output information
	bool flag;//"0" for read,"1"for write
public:
	void write();//Prepare information for Interpolation,write to the hard disk
	void read();//read information from file and conduct interpolation
	void init(void);//initial psi
	void Cal_Coord();//initial pos,and other parameters
	Vector3<int> Get_Index(Vector3<double> pos, Vector3<double> Origin, Vector3<double> length);
	int Sort_Partition(complex<double> **val, Vector3<double> *pos, int low, int high);
	void QSort_Indata(complex <double> **val, Vector3<double> *pos, int low, int high);
	complex<double> Interpolation_3D(complex<double> **psi, const double &x, const double &y, const double &z, const int seq);
	void Cal_Result();
};

#endif
