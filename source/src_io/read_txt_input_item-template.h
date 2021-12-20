//=======================
// AUTHOR : Peize Lin
// DATE :   2021-12-13
//=======================

#ifndef READ_TXT_INPUT_ITEM_TEMPLATE_H
#define READ_TXT_INPUT_ITEM_TEMPLATE_H

#include "read_txt_input_item.h"

namespace Read_Txt_Input
{
	template<typename T>
	void Input_Item::default_1(const T &value)
	{
		set_value(value);
		++this->values_size_lower_limit;
		++this->values_size_upper_limit;
	}

	template <typename T_head, typename... T_tail>
	void Input_Item::default_1( const T_head &value_head, const T_tail... values_tail)
	{
		set_value(value_head);
		default_1(values_tail...);
		++this->values_size_lower_limit;
		++this->values_size_upper_limit;
	}

	template<>
	void Input_Item::set_value(const bool &b)
	{
		this->values.push_back({});
		this->values.back().setb(b);
		this->values_type.push_back("b");
	}

	template<>
	void Input_Item::set_value(const int &i)
	{
		this->values.push_back({});
		this->values.back().seti(i);
		this->values_type.push_back("i");
	}

	template<>
	void Input_Item::set_value(const double &d)
	{
		this->values.push_back({});
		this->values.back().setd(d);
		this->values_type.push_back("d");
	}

	template<>
	void Input_Item::set_value(const std::string &s)
	{
		this->values.push_back({});
		this->values.back().sets(s);
		this->values_type.push_back("s");
	}

	template<>
	void Input_Item::set_value(const char*const &s)
	{
		this->values.push_back({});
		this->values.back().sets(std::string(s));
		this->values_type.push_back("s");
	}
}

#endif