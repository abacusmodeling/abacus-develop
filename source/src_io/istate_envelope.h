#ifndef ISTATE_ENVELOPE_H
#define ISTATE_ENVELOPE_H
#include "src_lcao/local_orbital_wfc.h"
#include "src_lcao/gint_gamma.h"

class IState_Envelope
{
	public:
	IState_Envelope();
	~IState_Envelope();

	void begin(Local_Orbital_wfc &lowf, Gint_Gamma &gg);

	private:
	bool *bands_picked;


};
#endif
