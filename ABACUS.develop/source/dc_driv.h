//==========================================================
// AUTHOR : mohan
// DATE : 2011-06-12
//==========================================================
#ifndef DC_DRIV_H
#define DC_DRIV_H

class DC_Driv
{
	public:
	
	DC_Driv();
	~DC_Driv();

	void init();

	private:

	// reading the parameters
	void reading();

    void prepare();

	// do stuff, have fun!
	void atomic_world();


};

#endif
