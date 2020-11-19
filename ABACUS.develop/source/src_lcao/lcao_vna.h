#ifndef LCAO_VNA_H
#define LCAO_VNA_H

class LCAO_Vna
{
	public:

	LCAO_Vna();
	~LCAO_Vna();

	static void smooth_vl1(void);

	void smooth_vl2(void);
	void dense_vna(const char &matrix_type);
	void two_center_vna(void);

	private:




};

#endif 
