void test_libxc_parameter()
{
	auto test_parameter_PBE0 = [&]()
	{
		TITLE("test_parameter_PBE0");
		xc_func_type func;
		xc_func_init( &func, XC_HYB_GGA_XC_PBEH, XC_POLARIZED );
		
		std::cout<<xc_hyb_exx_coef(&func)<<std::endl;						// 0.25
	};
	auto test_parameter_HSE = [&]()
	{
		TITLE("test_parameter_HSE");
		xc_func_type func;
		xc_func_init( &func, XC_HYB_GGA_XC_HSE06, XC_POLARIZED );
		
		double omega=-1,alpha=-1,beta=-1;
		std::cout<<xc_hyb_exx_coef(&func)<<std::endl;						// 0
		xc_hyb_cam_coef(&func, &omega, &alpha, &beta);
		std::cout<<omega<<"\t"<<alpha<<"\t"<<beta<<std::endl;					// 0.11, 0, 0.25
		
		double parameter_hse[3] = {0.123,0.456,0.789};				// beta, omega_HF, omega_PBE 
		xc_func_set_ext_params(&func, parameter_hse);
		
		omega=-2,alpha=-2,beta=-2;
		std::cout<<xc_hyb_exx_coef(&func)<<std::endl;						// 0
		xc_hyb_cam_coef(&func, &omega, &alpha, &beta);
		std::cout<<omega<<"\t"<<alpha<<"\t"<<beta<<std::endl;					// 0.456, 0, 0.123
	};	
}