#include "read_txt_input_list.h"

#include "src_io/read_txt_tools.h"
#include "module_base/global_variable.h"

#include "src_lcao/serialization_cereal.h"
#include <cereal/archives/binary.hpp>
#include <mpi.h>
#include <sstream>

namespace Read_Txt_Input
{
	void Input_List::add_item(const Input_Item &input_item)
	{
		list.insert(make_pair(input_item.label, input_item));
		add_order.push_back(input_item.label);
	}

	void Input_List::check_value_size(const Input_Item& input_item, const int size)
	{
		if(input_item.value_read_size==-1)
			return;
		if(input_item.value_read_size!=size)
			throw std::out_of_range(std::to_string(input_item.value_read_size)+std::to_string(size));
	}	

	void Input_List::check_value_size(const Input_Item& input_item, const int size_low, const int size_up)
	{
		if(input_item.value_read_size==-1)
			return;
		if(input_item.value_read_size<size_low || input_item.value_read_size>size_up)
			throw std::out_of_range(std::to_string(input_item.value_read_size)+std::to_string(size_low)+std::to_string(size_up));
	}

	void Input_List::read_and_convert(const std::string &file_name)
	{
		this->set_items();
		if(GlobalV::MY_RANK==0)
		{
			this->read(file_name);
			this->check_transform();
			this->default2();
			this->out(GlobalV::global_out_dir + file_name);
		}
		this->bcast();
		this->convert();
	}

	void Input_List::set_items()
	{
		set_items_general();
		set_items_pw();
	}

	void Input_List::read(const std::string &file_name)
	{
		const std::map<std::string, std::vector<std::string>> inputs_read
			= Read_Txt_Tools::read_file_to_map(file_name, {"#","\\"});
		for(const auto & input_read : inputs_read)
		{
			const auto ptr = list.find(input_read.first);
			if(ptr==list.end())
				throw std::out_of_range("input_read.first");
			if(input_read.second.size() > ptr->second.value.size())
				throw std::out_of_range("size error");
			for(size_t i=0; i<input_read.second.size(); ++i)
				ptr->second.value[i].s = input_read.second[i];
			ptr->second.value_read_size = input_read.second.size();
		}
	}

	void Input_List::check_transform()
	{
		for(auto &item : this->list)
			item.second.check_transform(item.second);
	}

	void Input_List::default2()
	{
		for(auto &item : this->list)
			item.second.default2(this->list);
	}

	void Input_List::out(const std::string &file_name) const
	{
		std::ofstream ofs(file_name);
		for(const std::string &label : add_order)
		{
			ofs<<label<<"\t";
			for(const Input_Value &value : this->list.at(label).value)
				ofs<<value.s<<" ";
			ofs<<"\t# "<<this->list.at(label).annotation<<std::endl;
		}
	}

	void Input_List::bcast()
	{
#ifdef USE_CEREAL_SERIALIZATION
		if(GlobalV::MY_RANK==0)
		{
			std::stringstream ss;
			{
				cereal::BinaryOutputArchive ar(ss);
				ar(this->list);
			}
			const int size = ss.str().size();
			MPI_Bcast( const_cast<int*>(&size), 1, MPI_INT, 0, MPI_COMM_WORLD );
			MPI_Bcast( const_cast<char*>(ss.str().c_str()), size, MPI_CHAR, 0, MPI_COMM_WORLD ); 
		}
		else
		{
			int size;
			MPI_Bcast( &size, 1, MPI_INT, 0, MPI_COMM_WORLD );
			std::vector<char> c(size);
			MPI_Bcast( c.data(), size, MPI_CHAR, 0, MPI_COMM_WORLD );   
			std::stringstream ss;  
			ss.rdbuf()->pubsetbuf(c.data(),size);
			{
				cereal::BinaryInputArchive ar(ss);
				ar(this->list);
			}
		}
#else
#error Input_List::bcast() needs cereal
#endif
	}

	void Input_List::convert()
	{
		for(auto &item : this->list)
			item.second.convert(item.second);
	}
}