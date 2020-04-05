#include "restart.h"
#include "src_global/global_function.h"
#include "src_pw/global.h"
#include "src_lcao/global_fp.h"
#include <fcntl.h> 
#include <unistd.h>
#include <fstream>
#include <stdexcept>

void Restart::write_file1(const std::string &file_name, const void*const ptr, const size_t size)
{
	ofstream ofs(file_name, ofstream::binary|ofstream::trunc);
	ofs.write(static_cast<const char*>(ptr),size);
}

void Restart::read_file1(const std::string &file_name, void*const ptr, const size_t size)
{
	ifstream ifs(file_name, ifstream::binary);
	ifs.read(static_cast<char*>(ptr),size);
}

void Restart::write_file2(const std::string &file_name, const void*const ptr, const size_t size)
{
	const int file = open(file_name.c_str(), O_WRONLY|O_CREAT|O_TRUNC, S_IRUSR|S_IWUSR);
	if(-1==file)	throw runtime_error("errno="+TO_STRING(errno)+".\n"+TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));
	write(file, ptr, size);
	close(file);
}

void Restart::read_file2(const std::string &file_name, void*const ptr, const size_t size)
{
	const int file = open(file_name.c_str(), O_RDONLY);
	if(-1==file)	throw runtime_error("errno="+TO_STRING(errno)+".\n"+TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));
	read(file, ptr, size);
	close(file);
}

void Restart::save_disk(const int i)
{
	if("charge"==info_save.mode)
		write_file2(info_save.mode+"_"+TO_STRING(MY_RANK)+"_"+TO_STRING(i), chr.rho[i], pw.nrxx*sizeof(double));
	else if("H"==info_save.mode)
	{
		if(GAMMA_ONLY_LOCAL)
			write_file2(info_save.mode+"_"+TO_STRING(MY_RANK)+"_"+TO_STRING(i), LM.Hloc, ParaO.nloc*sizeof(double));
		else
			write_file2(info_save.mode+"_"+TO_STRING(MY_RANK)+"_"+TO_STRING(i), LM.Hloc2, ParaO.nloc*sizeof(complex<double>));
	}
	else
		throw invalid_argument(TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));
}

void Restart::load_disk(const int i)
{
	if("charge"==info_load.mode)
		read_file2(info_load.mode+"_"+TO_STRING(MY_RANK)+"_"+TO_STRING(i), chr.rho[i], pw.nrxx*sizeof(double));
	else if("H"==info_load.mode)
	{
		if(GAMMA_ONLY_LOCAL)
			read_file2(info_load.mode+"_"+TO_STRING(MY_RANK)+"_"+TO_STRING(i), LM.Hloc, ParaO.nloc*sizeof(double));
		else
			read_file2(info_load.mode+"_"+TO_STRING(MY_RANK)+"_"+TO_STRING(i), LM.Hloc2, ParaO.nloc*sizeof(complex<double>));
	}
	else
		throw invalid_argument(TO_STRING(__FILE__)+" line "+TO_STRING(__LINE__));
}