#include "xc_functional.h"
#include "../src_pw/global.h"
#include "../module_base/global_function.h"

XC_Functional::XC_Functional(){}

XC_Functional::~XC_Functional(){}

std::vector<int> XC_Functional::func_id(1);
int XC_Functional::func_type = 0;
bool XC_Functional::use_libxc = true;

int XC_Functional::get_func_type()
{
    return func_type;
}

// The setting values of functional id according to the index in LIBXC
// for detail, refer to https://www.tddft.org/programs/libxc/functionals/
void XC_Functional::set_xc_type(const std::string xc_func_in)
{
    func_id.clear();
    std::string xc_func = xc_func_in;
    std::transform(xc_func.begin(), xc_func.end(), xc_func.begin(), (::toupper));
	if( xc_func == "LDA" || xc_func == "PZ" || xc_func == "SLAPZNOGXNOGC") //SLA+PZ
	{
        // I should use XC_LDA_X and XC_LDA_C_PZ here,
        // but since we are not compiling with libxc as default,
        // I chose to set it manually instead.
        // Same for the rest.
        func_id.push_back(XC_LDA_X);
        func_id.push_back(XC_LDA_C_PZ);
        func_type = 1;
        use_libxc = false;
	}
    else if (xc_func == "PWLDA")
    {
        func_id.push_back(XC_LDA_X);
        func_id.push_back(XC_LDA_C_PW);
        func_type = 1;
        use_libxc = false;
    }
	else if ( xc_func == "PBE" || xc_func == "SLAPWPBXPBC") //PBX+PBC
	{
        func_id.push_back(XC_GGA_X_PBE);
        func_id.push_back(XC_GGA_C_PBE);
        func_type = 2;
        use_libxc = false;
	} 
	else if( xc_func == "revPBE" ) //rPBX+PBC
	{
		func_id.push_back(XC_GGA_X_RPBE);
        func_id.push_back(XC_GGA_C_PBE);
        func_type = 2;
        use_libxc = false;
	}
	else if ( xc_func == "WC") //WC+PBC
	{
        func_id.push_back(XC_GGA_X_WC);
        func_id.push_back(XC_GGA_C_PBE);
        func_type = 2;
        use_libxc = false;
	}	
	else if ( xc_func == "BLYP") //B88+LYP
	{
        func_id.push_back(XC_GGA_X_B88);
        func_id.push_back(XC_GGA_C_LYP);
        func_type = 2;
        use_libxc = false;
	}
	else if ( xc_func == "BP") //B88+P86
	{
        func_id.push_back(XC_GGA_X_B88);
        func_id.push_back(XC_GGA_C_P86);
        func_type = 2;
        use_libxc = false;
	} 
	else if ( xc_func == "PW91") //PW91_X+PW91_C
	{
        func_id.push_back(XC_GGA_X_PW91);
        func_id.push_back(XC_GGA_C_PW91);
        func_type = 2;
        use_libxc = false;
	} 
	else if ( xc_func == "HCTH") //HCTH_X+HCTH_C
	{
        func_id.push_back(XC_GGA_X_HCTH_A);
        func_id.push_back(XC_GGA_C_HCTH_A);
        func_type = 2;
        use_libxc = false;
	}
	else if ( xc_func == "OLYP") //OPTX+LYP
	{
        func_id.push_back(XC_GGA_X_OPTX);
        func_id.push_back(XC_GGA_C_LYP);
        func_type = 2;
        use_libxc = false;
	}
#ifdef USE_LIBXC
	else if ( xc_func == "SCAN")
	{
        func_id.push_back(XC_MGGA_X_SCAN);
        func_id.push_back(XC_MGGA_C_SCAN);
        func_type = 3;
	}
#endif
   	else if( xc_func == "PBE0")
	{
        func_id.push_back(XC_HYB_GGA_XC_PBEH);
        func_type = 4;
        use_libxc = false;
	}
    else if( xc_func == "HF" || xc_func == "OPT_ORB" ||  xc_func == "NONE" || xc_func == "NOX+NOC")
    {
        // not doing anything
        if(xc_func == "HF") func_type = 4;
    }
#ifdef USE_LIBXC
    else if( xc_func == "HSE")
    {
        func_id.push_back(XC_HYB_GGA_XC_HSE06);
        func_type = 4;
        use_libxc = true;
    }
#endif
    else
    {
#ifdef USE_LIBXC
        //see if it matches libxc functionals
        set_xc_type_libxc(xc_func);
#else
        std::cout << "functional name not recognized!" << std::endl;
#endif
    }

	if (func_id[0] == XC_GGA_X_OPTX)
	{
		std::cerr << "\n OPTX untested please test,";
	}

    if(func_type == 3 && GlobalV::BASIS_TYPE != "pw")
    {
        ModuleBase::WARNING_QUIT("set_xc_type","mGGA not realized for LCAO yet");
    }
    if(func_type == 4 && GlobalV::BASIS_TYPE == "pw")
    {
        ModuleBase::WARNING_QUIT("set_xc_type","hybrid functional not realized for planewave yet");
    }
    if(func_type == 3 && GlobalV::CAL_STRESS == 1 && GlobalV::NSPIN!=1)
    {
        ModuleBase::WARNING_QUIT("set_xc_type","mgga stress not implemented for polarized case yet");
    }

#ifndef USE_LIBXC
    if(xc_func == "SCAN" || xc_func == "HSE")
    {
        ModuleBase::WARNING_QUIT("set_xc_type","to use SCAN or HSE, LIBXC is required");
    }
    use_libxc = false;
#endif

}

#ifdef USE_LIBXC
void XC_Functional::set_xc_type_libxc(std::string xc_func_in)
{

    // determine the type (lda/gga/mgga)
    func_type = 1;
    if(xc_func_in.find("GGA") != std::string::npos) func_type = 2;
    if(xc_func_in.find("MGGA") != std::string::npos) func_type = 3;
    if(xc_func_in.find("HYB") != std::string::npos) func_type =4;

    // determine the id
    int pos = 0;
    std::string delimiter = "+";
    std::string token;
    while ((pos = xc_func_in.find(delimiter)) != std::string::npos)
    {
        token = xc_func_in.substr(0, pos);
        int id = xc_functional_get_number(token.c_str());
        std::cout << "func,id" << token << " " << id << std::endl;
        if (id == -1) ModuleBase::WARNING_QUIT("XC_Functional::set_xc_type_libxc","functional name not recognized!");
        func_id.push_back(id);
        xc_func_in.erase(0, pos + delimiter.length());
    }
    int id = xc_functional_get_number(xc_func_in.c_str());
    std::cout << "func,id" << xc_func_in << " " << id << std::endl;
    if (id == -1) ModuleBase::WARNING_QUIT("XC_Functional::set_xc_type_libxc","functional name not recognized!");
    func_id.push_back(id);

}
#endif

#ifdef USE_LIBXC
std::vector<xc_func_type> XC_Functional::init_func(const int xc_polarized)
{
	// 'funcs' is the return value
	std::vector<xc_func_type> funcs;

	//-------------------------------------------
	// define a function named 'add_func', which 
	// initialize a functional according to its ID
	//-------------------------------------------
	auto add_func = [&]( const int func_id )
	{
		funcs.push_back({});
		// 'xc_func_init' is defined in Libxc
		xc_func_init( &funcs.back(), func_id, xc_polarized );
	};

	for(int id : func_id)
	{
		if( id == 406 ) // PBE0
		{
			add_func( XC_HYB_GGA_XC_PBEH );		
			double parameter_hse[3] = { GlobalC::exx_global.info.hybrid_alpha, 
				GlobalC::exx_global.info.hse_omega, 
				GlobalC::exx_global.info.hse_omega };
			xc_func_set_ext_params(&funcs.back(), parameter_hse);	
		}
		else if( id == 428 ) // HSE06 hybrid functional
		{
			add_func( XC_HYB_GGA_XC_HSE06 );	
			double parameter_hse[3] = { GlobalC::exx_global.info.hybrid_alpha, 
				GlobalC::exx_global.info.hse_omega, 
				GlobalC::exx_global.info.hse_omega };
			xc_func_set_ext_params(&funcs.back(), parameter_hse);
		}
		else
		{
			add_func( id );
		}
	}
	return funcs;
}
#endif