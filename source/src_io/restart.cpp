#include "restart.h"
#include "../module_base/global_function.h"
#include "../src_pw/global.h"
#ifdef __LCAO
#include "../src_lcao/global_fp.h"
#endif
#include <fcntl.h> 
#include <unistd.h>
#include <fstream>
#include <stdexcept>

void Restart::write_file1(const std::string &file_name, const void*const ptr, const size_t size) const
{
	std::ofstream ofs(file_name, std::ofstream::binary|std::ofstream::trunc);
	ofs.write(static_cast<const char*>(ptr),size);
}

void Restart::read_file1(const std::string &file_name, void*const ptr, const size_t size) const
{
	std::ifstream ifs(file_name, std::ifstream::binary);
	ifs.read(static_cast<char*>(ptr),size);
}

void Restart::write_file2(const std::string &file_name, const void*const ptr, const size_t size) const
{
	const int file = open(file_name.c_str(), O_WRONLY|O_CREAT|O_TRUNC, S_IRUSR|S_IWUSR);
	if(-1==file)	throw std::runtime_error("can't open restart save file. \nerrno="+ModuleBase::GlobalFunc::TO_STRING(errno)+".\n"+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
	write(file, ptr, size);
	close(file);
}

void Restart::read_file2(const std::string &file_name, void*const ptr, const size_t size) const
{
	const int file = open(file_name.c_str(), O_RDONLY);
	if(-1==file)	throw std::runtime_error("can't open restart load file. \nerrno="+ModuleBase::GlobalFunc::TO_STRING(errno)+".\n"+ModuleBase::GlobalFunc::TO_STRING(__FILE__)+" line "+ModuleBase::GlobalFunc::TO_STRING(__LINE__));
	read(file, ptr, size);
	close(file);
}

void Restart::save_disk(const std::string mode, const int i) const
{
	if("charge"==mode)
		write_file2(folder+"charge_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK)+"_"+ModuleBase::GlobalFunc::TO_STRING(i), GlobalC::CHR.rho[i], GlobalC::rhopw->nrxx*sizeof(double));
}

void Restart::load_disk(const std::string mode, const int i) const
{
	if("charge"==mode)
		read_file2(folder+"charge_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK)+"_"+ModuleBase::GlobalFunc::TO_STRING(i), GlobalC::CHR.rho[i], GlobalC::rhopw->nrxx*sizeof(double));
}

#ifdef __LCAO
void Restart::save_disk(LCAO_Matrix &lm, const std::string mode, const int i) const
{
	if("charge"==mode)
		write_file2(folder+"charge_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK)+"_"+ModuleBase::GlobalFunc::TO_STRING(i), GlobalC::CHR.rho[i], GlobalC::rhopw->nrxx*sizeof(double));
	if("H"==mode)
	{
		if(GlobalV::GAMMA_ONLY_LOCAL)
			write_file2(folder+"Hgamma_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK)+"_"+ModuleBase::GlobalFunc::TO_STRING(i), lm.Hloc.data(), lm.ParaV->nloc*sizeof(double));
		else
			write_file2(folder+"Hk_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK)+"_"+ModuleBase::GlobalFunc::TO_STRING(i), lm.Hloc2.data(), lm.ParaV->nloc*sizeof(std::complex<double>));
	}
}

void Restart::load_disk(LCAO_Matrix &lm, const std::string mode, const int i) const
{
	if("charge"==mode)
		read_file2(folder+"charge_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK)+"_"+ModuleBase::GlobalFunc::TO_STRING(i), GlobalC::CHR.rho[i], GlobalC::rhopw->nrxx*sizeof(double));
	if("H"==mode)
	{
		if(GlobalV::GAMMA_ONLY_LOCAL)
			read_file2(folder+"Hgamma_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK)+"_"+ModuleBase::GlobalFunc::TO_STRING(i), lm.Hloc.data(), lm.ParaV->nloc*sizeof(double));
		else
			read_file2(folder+"Hk_"+ModuleBase::GlobalFunc::TO_STRING(GlobalV::MY_RANK)+"_"+ModuleBase::GlobalFunc::TO_STRING(i), lm.Hloc2.data(), lm.ParaV->nloc*sizeof(std::complex<double>));
	}
}
#endif