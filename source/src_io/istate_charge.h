#ifndef ISTATE_CHARGE_H
#define ISTATE_CHARGE_H

class IState_Charge
{
	public:

	IState_Charge();
	~IState_Charge();

	void begin();

	private:

	int *bands_picked;

#ifdef __MPI
	void idmatrix(const int &ib);
#endif

};
#endif
