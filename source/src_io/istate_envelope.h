#ifndef ISTATE_ENVELOPE_H
#define ISTATE_ENVELOPE_H

class IState_Envelope
{
	public:
	IState_Envelope();
	~IState_Envelope();

	void begin();

	private:
	bool *bands_picked;


};
#endif
