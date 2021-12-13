//=======================
// AUTHOR : Peize Lin
// DATE :   2021-12-13
//=======================

#include "read_txt_input_process.h"

#include "src_io/read_txt_tools.h"
#include "module_base/global_variable.h"

#include <mpi.h>
#include <sstream>

#ifdef USE_CEREAL_SERIALIZATION
#include "src_lcao/serialization_cereal.h"
#include <cereal/archives/binary.hpp>
#endif

namespace Read_Txt_Input
{
	void Input_Process::read_and_convert(const std::string &file_name)
	{
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

	void Input_Process::read(const std::string &file_name)
	{
		const std::map<std::string, std::vector<std::string>> inputs_read
			= Read_Txt_Tools::read_file_to_map(file_name, {"#","\\"});
		for(const auto & input_read : inputs_read)
		{
			const auto ptr = input.list.find(input_read.first);
			if(ptr==input.list.end())
				throw std::out_of_range("input_read.first");

			//if(input_read.second.size() > ptr->second.values.size())
			if(input_read.second.size() != ptr->second.values.size())
				throw std::out_of_range(
					"size of INPUT "+ptr->second.label+" needs "+std::to_string(ptr->second.values.size())
					+" but reads "+std::to_string(input_read.second.size()));

			for(size_t i=0; i<input_read.second.size(); ++i)
				ptr->second.values[i].sets(input_read.second[i]);
			//ptr->second.value_read_size = input_read.second.size();
		}
	}

	void Input_Process::check_transform()
	{
		for(auto &item : input.list)
		{
			item.second.check_transform(item.second);

			for(size_t i=0; i<item.second.values.size(); ++i)
			{
				if(item.second.values_type[i]=="b")
				{
					if(!Read_Txt_Tools::in_set( item.second.values[i].gets(), Read_Txt_Tools::Preset::Bool))
						throw std::invalid_argument("INPUT "+item.second.label+" must be bool");
				}
			}
		}
	}

	void Input_Process::default2()
	{
		for(auto &item : input.list)
			item.second.default2(input.list);
	}

	void Input_Process::out(const std::string &file_name) const
	{
		std::ofstream ofs(file_name);
		for(const std::string &label : input.add_order)
		{
			ofs<<label<<"\t";
			const Read_Txt_Input::Input_Item &item = input.list.at(label);
			for(size_t i=0; i<item.values.size(); ++i)
			{
				if(item.values_type[i]=="s")
					ofs<<item.values[i].gets()<<" ";
				else if(item.values_type[i]=="d")
					ofs<<std::to_string(item.values[i].getd())<<" ";
				else if(item.values_type[i]=="i")
					ofs<<std::to_string(item.values[i].geti())<<" ";
				else if(item.values_type[i]=="b")
				{
					if(item.values[i].getb())
						ofs<<"true"<<" ";
					else
						ofs<<"false"<<" ";
				}
				else
					throw invalid_argument("Input_Process::out() value_type["+std::to_string(i)+"]="+item.values_type[i]);
			}
			ofs<<"\t# "<<item.annotation<<std::endl;
		}
	}

	void Input_Process::bcast()
	{
#ifdef USE_CEREAL_SERIALIZATION
		if(GlobalV::MY_RANK==0)
		{
			std::stringstream ss;
			{
				cereal::BinaryOutputArchive ar(ss);
				ar(input.list);
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
				ar(input.list);
			}
		}
#else
#error Input_Process::bcast() needs cereal
#endif
	}

	void Input_Process::convert()
	{
		for(auto &item : input.list)
			item.second.convert(item.second);
	}
}