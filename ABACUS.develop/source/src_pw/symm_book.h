#ifndef SYMM_BOOK_H
#define SYMM_BOOK_H

#include <string>
using namespace std;

class Symm_Book
{
	public:
	Symm_Book();
	~Symm_Book();

	void Bravais();

	private:

	void write(const string &name);

};



#endif
