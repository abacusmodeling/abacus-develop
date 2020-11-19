#include "exx_pw.h"
#include "global.h"

Exx_pw::Exx_pw()
{
	ik_now = 0;	
	alpha = 0.0;
	start = false;
}

Exx_pw::~Exx_pw()
{}

void Exx_pw::init(bool start_in)
{
	TITLE("Exx_pw","init");
	this->start = start_in;
	this->exxdiv.init();
	if(DFT_FUNCTIONAL == "PBE0")
	{
		this->alpha = 0.25;
	}
	return;
}

void Exx_pw::get_exx(void)
{
	TITLE("Exx_pw","get_exx");
	timer::tick("Exx_pw","get_exx");

	en.exx = 0.0;

	for(int ik=0; ik<kv.nks; ik++)
	{
		this->ik_now = ik;
		wf.npw = kv.ngk[ik];
		if(NSPIN==2) CURRENT_SPIN = kv.isk[ik];
			
		complex<double> *vpsi = new complex<double>[wf.npw];
		complex<double> *psi = new complex<double>[wf.npw];
		
		for(int ib=0; ib<NBANDS; ib++)
		{
			for(int ig=0; ig<wf.npw; ig++)
			{
				psi[ig] = wf.evc[ik](ib, ig);
			}	
		
			ZEROS(vpsi, wf.npw);
			this->vxx_psi(psi, vpsi);
		
			for(int ig=0; ig<wf.npw; ig++)
			{
				en.exx = en.exx - 0.5 * wf.wg(ik, ib) *  (conj(psi[ig]) * vpsi[ig]).real(); 
			}
		}

		delete[] psi;
		delete[] vpsi;
	}
	
	cout << " Exx = " << en.exx << " (Rydberg)" << endl;
	timer::tick("Exx_pw","get_exx");
	return;
}

void Exx_pw::get_exx2(void)
{
	if(!start) return;
	timer::tick("Exx_pw","get_exx2");
	complex<double> *psi = new complex<double>[pw.nrxx];
	complex<double> *phi = new complex<double>[pw.nrxx];
	complex<double> *prho = new complex<double>[pw.nrxx];
	en.exx = 0.0;

	for(int ik=0; ik<kv.nks; ik++)
	{
		this->ik_now = ik;
		wf.npw = kv.ngk[ik];
		for(int ib=0; ib<NBANDS; ib++)
		{
			// wave function transfer
			ZEROS(psi, pw.nrxx);
			for(int ig=0; ig<wf.npw; ig++)
			{
				psi[ pw.ig2fftw[ wf.igk(ik,ig) ] ] = wf.evc[ik](ib,ig);
			}
			
			// to real space
			pw.FFT_wfc.FFT3D(psi, 1);

			for(int iq=0; iq<kv.nkstot; iq++)
			{
				this->exxdiv.get_factor(ik, iq);
				for(int m=0; m<NBANDS; m++)
				{
					// wave function transfer
					/*
					ZEROS(phi, pw.nrxx);
					for(int ig=0; ig<kv.ngk[iq]; ig++)
					{
						phi[ pw.ig2fftw[ wf.igk(iq, ig) ] ] = wf.evc[iq](m,ig);
					}

					// to real space
					pw.FFT_wfc.FFT3D(phi, 1);
					*/
					for(int ir=0; ir<pw.nrxx; ir++)
					{
						phi[ir] = ONE;
					}
					
					// calculate the pair density in real space. <phi_nq|psi_mk>
					for(int ir=0; ir<pw.nrxx; ir++)
					{
						prho[ir] = conj(psi[ir]) * phi[ir] / ucell.omega;
					}

					// to G space
					pw.FFT_wfc.FFT3D(prho, -1);
				
					double vc = 0.0;
					for(int ig=0; ig<wf.npw; ig++)
					{
						const complex<double> rho = prho[ pw.ig2fftw[ wf.igk(ik,ig) ] ];
						vc += exxdiv.factor[ig] * ( rho * conj(rho) ).real();	
					}
					cout << " ib=" << ib << " m=" << m <<  " vc=" << vc << endl;;

					en.exx = en.exx - this->alpha * vc * ucell.omega * wf.wg(ik,ib) / (double) kv.nkstot; 
				} // end m
			} // end iq
		}// end ib
	} // end ik
	delete[] psi;
	delete[] phi;
	delete[] prho;

	cout << " Exx = " << en.exx << " (Rydberg)" << endl;
	en.exx *= -0.5;
	timer::tick("Exx_pw","get_exx2");
	return;
}

void Exx_pw::vxx_psi(complex<double>* psi_in, complex<double> *hpsi)
{
	return;
	if(!start) return;
	timer::tick("Exx_pw","vxx_psi");

	// SUMMARAY : 	
	// Vxx|psi(kn)> = \int dr sum_{mq} psi_mq(r') * psi_mq(r) * phi_kn(r) / |r-r'|
	// where the last term, <psi(mq)|phi(kn)> is pair density(pdensity).

	//(2) FFT(psi)
	complex<double> *psi = new complex<double>[pw.nrxx]; 
	ZEROS(psi, pw.nrxx);
	for(int ig=0; ig<wf.npw; ig++)
	{
		psi[ pw.ig2fftw[ wf.igk(ik_now,ig) ] ] = psi_in[ig];
	}
	
	// This part is :
	// Vpair(r') = \int dr ( psi_mq(r) * phi_kn(r) ) / |r-r'|
	// which is done in G space instead, then FFT back to real space.

	complex<double> *phi = new complex<double>[pw.nrxx];
	complex<double> *prho = new complex<double>[pw.nrxx];	
	complex<double> *denG = new complex<double>[wf.npw];
	const double e2 = 2.0;
	
	for (int iq=0; iq<kv.nkstot; iq++)
	{
		// (3.0) set the factor |G + k - q|
		this->exxdiv.get_factor(this->ik_now, iq);
		for(int n=0; n<NBANDS; n++)
		{	
			// (3.1) doing fft, move the phi in Vxx to real space.
			UFFT.ToRealSpace_psi(iq, n, wf.evc[iq], phi);

			// (3.2) calculate the pair density in real space. <phi_nq|psi_mk>
			for(int ir=0; ir<pw.nrxx; ir++)
			{
				prho[ir] = conj(psi[ir]) * phi[ir] / ucell.omega;
			}
				
			// (3.3) bring the pair density back to G space.
			pw.FFT_wfc.FFT3D(prho, -1);
			for(int ig=0; ig<wf.npw; ig++)
			{
				denG[ig] = prho[ pw.ig2fftw[ wf.igk(ik_now, ig) ] ];
			}	

			// (3.4) multiply the pair density in G space by 1/|k-q+G|^2
			for (int ig=0; ig<wf.npw; ig++)
			{
				denG[ig] = this->exxdiv.factor[ig] * denG[ig] * wf.wg(iq,n) / (double)kv.nkstot;
			}

			// (3.5) FFT prho/|r-r'| back to real space
			ZEROS(prho, pw.nrxx);
			for(int ig=0; ig<wf.npw; ig++)
			{
				prho[ pw.ig2fftw[ wf.igk(ik_now, ig) ] ] = denG[ig];
			}
			pw.FFT_wfc.FFT3D(prho, 1);
			
			// (3.6) phi(r') * V(r')	
			for(int ir=0; ir<pw.nrxx; ir++)
			{
				prho[ir] *= phi[ir];
			}
	
			// (3.7)
			pw.FFT_wfc.FFT3D(prho, -1);
			for(int ig=0; ig<wf.npw; ig++)
			{
				hpsi[ig] -= this->alpha * prho[ pw.ig2fftw[ wf.igk(ik_now, ig) ] ];
			}
	
			// Vxx_psi(r') = Vpari(r') * psi_{mq}(r')
			// note: Vpair is diagnolized in real space, but not in G space,
			// this is the same as Vlocal.
			// get Vxx_psi(G) according to FFT( Vxx_psi(r') ).	 
		}
	}

	delete[] denG;
	delete[] prho;
	delete[] psi;
	delete[] phi;
	timer::tick("Exx_pw","vxx_psi");
	return;
}


