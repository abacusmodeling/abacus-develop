#ifndef RESTART_H
#define RESTART_H

#include <string>

class Restart
{
public:
	struct Info_Save
	{
		std::string mode = "no";
	};
	Info_Save info_save;
	
	struct Info_Load
	{
		std::string mode = "no";
		bool finish = false;
		bool exx = false;
	};
	Info_Load info_load;
	
	void save_disk(const int i);
	void load_disk(const int i);
	
private:
	void write_file1(const std::string &file_name, const void*const ptr, const size_t size);
	void read_file1(const std::string &file_name, void*const ptr, const size_t size);
	void write_file2(const std::string &file_name, const void*const ptr, const size_t size);
	void read_file2(const std::string &file_name, void*const ptr, const size_t size);
};

#endif