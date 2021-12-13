#include "read_txt_input.h"

namespace Read_Txt_Input
{
	void Input_Item::default1(const std::vector<std::string> &value_in)
	{
		this->value.resize(value_in.size());
		for(size_t i=0; i<value_in.size(); ++i)
			this->value[i].s = value_in[i];
	}
}