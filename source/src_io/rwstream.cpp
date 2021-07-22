#include <stdio.h>
#include "rwstream.h"

Rwstream::Rwstream(FILE* &ptr)
{
	fileptr=ptr;
}

Rwstream::Rwstream(const string filename,const char *op)
{
	fileptr=fopen(filename.c_str(),op);
}

void Rwstream:: setptr(FILE* &ptr)
{
	fileptr=ptr;
	return;
}

void Rwstream::open(const string filename,const char *op)
{
	fileptr=fopen(filename.c_str(),op);
}

void Rwstream:: close()
{
	fclose(fileptr);
	return;
}

bool Rwstream::operator!() const
{
	if (fileptr==NULL)
		return true;
	else
		return false;
}
Rwstream::operator bool() const
{
	if (fileptr==NULL)
		return false;
	else
		return true;
}
