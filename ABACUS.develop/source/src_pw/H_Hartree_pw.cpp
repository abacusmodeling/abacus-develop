//--------------------------------------------------------------------
// Transform charge density to hartree potential.
//--------------------------------------------------------------------
void potential::v_h(int NSPIN,double &ehart, matrix &v, double** rho)
{
    TITLE("potential","v_h");
    timer::tick("potential","v_hartree");

    complex<double> *Porter = UFFT.porter;

    //  Hartree potential VH(r) from n(r)
    ZEROS( Porter, pw.nrxx );
    int nspin0 = 1;
    if(NSPIN==2)nspin0 = NSPIN;
    for(int is=0; is<nspin0; is++)
    {
        for (int ir=0; ir<pw.nrxx; ir++) 
        {
            Porter[ir] += complex<double>( rho[is][ir], 0.0 );
        }
    }

    //=============================
    //  bring rho (aux) to G space
    //=============================
    pw.FFT_chg.FFT3D(Porter, -1);

    //double charge;
    //if (pw.gstart == 1)
    //{
    //    charge = ucell.omega * Porter[pw.ig2fftc[0]].real();
    //}
    //OUT(ofs_running, "v_h charge", charge);

    //=======================================================
    // calculate hartree potential in G-space (NB: V(G=0)=0 )
    //=======================================================
    ehart = 0.0;

    complex<double> *vh_g  = new complex<double>[pw.ngmc];
    ZEROS(vh_g, pw.ngmc);

    for (int ig = pw.gstart; ig<pw.ngmc; ig++)
    {
        const int j = pw.ig2fftc[ig];
        if(pw.gg[ig] >= 1.0e-12) //LiuXh 20180410
        {
            const double fac = e2 * FOUR_PI / (ucell.tpiba2 * pw.gg [ig]);

            ehart += ( conj( Porter[j] ) * Porter[j] ).real() * fac;
            vh_g[ig] = fac * Porter[j];
        }
    }

    Parallel_Reduce::reduce_double_pool( ehart );
    ehart *= 0.5 * ucell.omega;

    //cout << " ehart=" << ehart << endl;

    ZEROS(Porter, pw.nrxx);

    for (int ig = 0;ig < pw.ngmc;ig++)
    {
        Porter[pw.ig2fftc[ig]] = vh_g[ig];
    }

    //==========================================
    //transform hartree potential to real space
    //==========================================
    pw.FFT_chg.FFT3D(Porter, 1);
    //==========================================
    //Add hartree potential to the xc potential
    //==========================================

	if(NSPIN==4)
		for (int ir = 0;ir < pw.nrxx;ir++)
		{
			v(0, ir) += Porter[ir].real();
		}
	else
    	for (int is = 0;is < NSPIN;is++)
    	{
        	for (int ir = 0;ir < pw.nrxx;ir++)
        	{
            		v(is, ir) += Porter[ir].real();
        	}
    	}

	//-------------------------------------------
	// output the Hartree potential into a file.
	//-------------------------------------------
	if(out_potential==-2)
	{
		cout << " output VH" << endl;
		int is = 0;
		int iter = 0;
		int precision = 3;
		string fn = "VH.dat";
		stringstream ss;
		ss << global_out_dir << fn;
		matrix v;
		v.create(1,pw.nrxx);
		for(int ir=0; ir<pw.nrxx; ++ir)
		{
			v(0,ir) = Porter[ir].real();
		}
		this->write_potential( is, iter, ss.str(), v, precision, 1 );
	}

    timer::tick("potential","v_hartree");
    delete[] vh_g;
    return;
} // end subroutine v_h

