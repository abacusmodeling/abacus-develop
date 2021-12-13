#include "read_txt_input_value.h"

#include <stdexcept>

namespace Read_Txt_Input
{
	void Input_Value::setb(const bool &b_in)
	{
		this->s = std::to_string(b_in);
		this->d = static_cast<double>(b_in);
		this->i = static_cast<int>(b_in);
		this->b = b_in;
	}
	
	void Input_Value::seti(const int &i_in)
	{
		this->s = std::to_string(i_in);
		this->d = static_cast<double>(i_in);
		this->i = i_in;
		this->b = static_cast<bool>(i_in);
	}
	
	void Input_Value::setd(const double &d_in)
	{
		this->s = std::to_string(d_in);
		this->d = d_in;
		this->i = static_cast<int>(d_in);
		this->b = static_cast<bool>(d_in);
	}
	
	void Input_Value::sets(const std::string &s_in)
	{
		try
		{
			this->s = s_in;
			this->d = std::stod(s_in);
			this->i = std::stoi(s_in);
			this->b = static_cast<bool>(std::stoi(s_in));
		}
		catch(const std::invalid_argument &e)
		{
			this->s = s_in;
		}
	}
}