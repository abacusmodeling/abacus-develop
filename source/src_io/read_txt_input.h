#ifndef READ_TXT_INPUT_H
#define READ_TXT_INPUT_H

#include <vector>
#include <map>
#include <string>
#include <functional>

#ifdef USE_CEREAL_SERIALIZATION
#include <cereal/access.hpp>
#endif

namespace Read_Txt_Input
{
	struct Input_Value
	{
		bool b;
		int i;
		double d;
		std::string s;
	};

	class Input_Item
	{
	public:
		// call these functions
		Input_Item(const std::string &label_in) :label(label_in) {}
		void default1(const std::vector<std::string> &value_in);

		// set these variables and functions
		std::string annotation;
		std::function<void(Input_Item&)> check_transform
			= [](Input_Item&self){};
		std::function<void(const std::map<std::string, Input_Item>&)> default2
			= [](const std::map<std::string, Input_Item>&list){};
		std::function<void(const Input_Item&)> convert
			= [](const Input_Item&item){};

		std::vector<Input_Value> value;

	private:
		std::string label;
		int value_read_size = -1;
		friend class Input_List;

#ifdef USE_CEREAL_SERIALIZATION
	public:
		Input_Item() = default;	
	private:
		friend class cereal::access;
		template <class Archive> void serialize( Archive & ar );
#endif
	};	
}

#endif