#include "read_txt_input_list.h"

#include "src_io/read_txt_tools.h"
#include <stdexcept>

// for convert
#include "module_base/global_variable.h"
#include "src_pw/global.h"

namespace Read_Txt_Input
{
	void Input_List::set_items_general()
	{
		{
			Input_Item item("nspin");
			item.default_1(1);
			item.annotation = "1: single spin; 2: up and down spin; 4: noncollinear spin";
			item.check_transform = [](Input_Item &self)
			{
				Input_List::check_value_size(self, 1);
				if(!Read_Txt_Tools::in_set(self.values[0].gets(), {"1","2","4"}))
					throw std::invalid_argument("nspin must be 1,2,4");
			};
			item.convert = [](const Input_Item &self)
			{
				GlobalV::NSPIN = self.values[0].geti();
			};
			this->add_item(item);
		}

		{
			Input_Item item("ntype");
			item.default_1(0);
			item.annotation = "atom species number";
			item.check_transform = [](Input_Item &self)
			{
				Input_List::check_value_size(self, 1);
				if(self.values[0].geti()<=0)
					throw std::invalid_argument("ntype must > 0");
			};
			item.convert = [](const Input_Item &self)
			{
				GlobalC::ucell.ntype = self.values[0].geti();
			};
			this->add_item(item);
		}
	}
}