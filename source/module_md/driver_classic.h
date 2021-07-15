#ifndef DRIVER_CLASSIC_H
#define DRIVER_CLASSIC_H

class Driver_classic
{
	public:
	
	Driver_classic();
	~Driver_classic();

	private:

    void init();

	// reading the parameters
	void reading();

    // convert INPUT parameters for classic MD
    void convert();

	// classic MD
	void classic_world();


};

#endif