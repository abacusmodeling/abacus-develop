#ifndef READ_TXT_INIPUT_PROCESS_H
#define READ_TXT_INIPUT_PROCESS_H

#include "read_txt_input.h"
#include "read_txt_input_list.h"

namespace Read_Txt_Input
{
	class Input_Process
	{
	public:
		Input_Process(Input_List &input_in) :input(input_in) {}
		void read_and_convert(const std::string &file_name);

	private:
		Input_List &input;
		void read(const std::string &file_name);
		void check_transform();
		void default2();
		void out(const std::string &file_name) const;
		void bcast();
		void convert();		
	};
}

#endif