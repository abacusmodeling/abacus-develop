#ifndef READ_TXT_INIPUT_LIST_H
#define READ_TXT_INIPUT_LIST_H

#include "read_txt_input.h"

#include <map>
#include <string>

namespace Read_Txt_Input
{
	class Input_List
	{
	public:
		void read_and_convert(const std::string &file_name);

	private:	
		std::map<std::string, Input_Item> list;
		void add_item(const Input_Item &input_item);
		static void check_value_size(const Input_Item& input_item, const int size);
		static void check_value_size(const Input_Item& input_item, const int size_low, const int size_up);

		void set_items();
		void read(const std::string &file_name);
		void check_transform();
		void default2();
		void out(const std::string &file_name) const;
		void bcast();
		void convert();
		std::vector<std::string> add_order;

		void set_items_general();
		void set_items_pw();
	};
}

#endif