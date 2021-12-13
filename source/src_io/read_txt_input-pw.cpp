#include "read_txt_input_list.h"

#include "module_base/constants.h"
#include <stdexcept>

// for convert
#include "module_base/global_variable.h"

namespace Read_Txt_Input
{
	void Input_List::set_items_pw()
	{
		{
			Input_Item item("ecutwfc");
			item.default1({"100","eV"});
			item.annotation = "energy cutoff for wave functions";
			item.check_transform = [](Input_Item &self)
			{
				Input_List::check_value_size(self, 1, 2);
				if(self.value[1].gets()=="eV"){}
				else if(self.value[1].gets()=="Ry")
					self.value[0].setd( self.value[0].getd() / ModuleBase::Ry_to_eV );
				else
					throw std::invalid_argument(self.value[1].gets());
				self.value[1].sets("eV");
				if(self.value[0].getd()<=0)
					throw std::invalid_argument("ecutwfc must > 0");
			};
			item.convert = [](const Input_Item &self)
			{
				// ?? GlobalC::pw.set()
			};
			this->add_item(item);
		}
	}
}