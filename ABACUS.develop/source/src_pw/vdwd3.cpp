//==========================================================
// AUTHOR : Yuyang Ji
// DATE : 2019-04-22
// UPDATE : 2021-4-19
//==========================================================

#include "vdwd3.h"
#include "src_global/global_function.h"
#include "src_global/constants.h"
#include "src_global/element_name.h"

Vdwd3::Vdwd3(const UnitCell_pseudo &unit_in, Vdwd3_Parameters &para_in):
	ucell(unit_in),
	para(para_in){}

void Vdwd3::set_criteria(double &rthr, vector<Vector3<double>> &lat, vector<double> &tau_max)
{
    tau_max.resize(3);
	double r_cutoff = sqrt(rthr);
	Vector3<double> norm1 = (lat[1] ^ lat[2]).normalize();
	Vector3<double> norm2 = (lat[2] ^ lat[0]).normalize();
	Vector3<double> norm3 = (lat[0] ^ lat[1]).normalize();
	double cos10 = norm1 * lat[0];
	double cos21 = norm2 * lat[1];
	double cos32 = norm3 * lat[2];
	tau_max[0] = abs(r_cutoff/cos10);
	tau_max[1] = abs(r_cutoff/cos21);
	tau_max[2] = abs(r_cutoff/cos32);
}

void Vdwd3::init(const UnitCell_pseudo &ucell)
{
    lat.resize(3);
    lat[0] = ucell.a1*ucell.lat0;
    lat[1] = ucell.a2*ucell.lat0;
    lat[2] = ucell.a3*ucell.lat0;

    vector<double>atom_kind = atomkind(ucell);
    iz.reserve(ucell.nat);
	xyz.reserve(ucell.nat);
    for(size_t it=0; it!=ucell.ntype; it++)
        for(size_t ia=0; ia!=ucell.atoms[it].na; ia++)
		{
            iz.push_back(atom_kind[it]);
			xyz.push_back(ucell.atoms[it].tau[ia] * ucell.lat0);
		}

    vector<double> tau_max(3);
    if(para.model=="radius")
	{
        rep_vdw.resize(3);
		set_criteria(para.rthr2,lat,tau_max);
		for(size_t i=0 ; i<3 ; i++)
			rep_vdw[i] = ceil(tau_max[i]);
	}
	else if(para.model=="period")
		rep_vdw = {para.period.x, para.period.y, para.period.z};
	rep_cn.resize(3);
    set_criteria(para.cn_thr2,lat,tau_max);
    for(size_t i=0 ; i<3 ; i++)
        rep_cn[i] = ceil(tau_max[i]);
}

vector<double> Vdwd3::atomkind (const UnitCell_pseudo &ucell)
{
	vector<double>atom_kind(ucell.ntype);
	for(size_t i=0 ; i!=ucell.ntype ; i++)
		for(int j=0; j!=element_name.size(); j++)
			if (ucell.atoms[i].psd == element_name[j])
			{
				atom_kind[i] = j;
				break;
			}
	return atom_kind;
}

void Vdwd3::getc6(int &iat, int &jat, double &nci, double &ncj, double &c6)
{
	double c6mem = -1e99, rsum = 0.0, csum = 0.0, r_save = 1e99;
	double cn1 = 0.0, cn2 = 0.0, r = 0.0;
	for(size_t i=0; i!=para.mxc[iat]; i++)
		for(size_t j=0; j!=para.mxc[jat]; j++)
		{
			c6 = para.c6ab[0][j][i][jat][iat];
			if(c6 > 0)
			{
				cn1 = para.c6ab[1][j][i][jat][iat];
				cn2 = para.c6ab[2][j][i][jat][iat];
				r = pow((cn1-nci), 2)+pow((cn2-ncj), 2);
				if(r < r_save)
				{
					r_save = r;
					c6mem = c6;
				}
				double tmp1 = exp(para.k3*r);
				rsum += tmp1;
				csum += tmp1*c6;
			}
		}
	c6 = (rsum > 1e-99) ? csum/rsum : c6mem;
}

void Vdwd3::pbcncoord(vector<double> &cn)
{
	for(size_t i=0; i!=ucell.nat; i++)
	{
		double xn = 0.0;
		Vector3<double> tau;
		double r2 = 0.0, rr = 0.0;
		for(size_t iat=0; iat!=ucell.nat; iat++)
			for(int taux=-rep_cn[0]; taux<=rep_cn[0]; taux++)
				for(int tauy=-rep_cn[1]; tauy<=rep_cn[1]; tauy++)
					for(int tauz=-rep_cn[2]; tauz<=rep_cn[2]; tauz++)
					{
						if(iat==i && taux==0 && tauy==0 && tauz==0) continue;
						tau = static_cast<double>(taux)*lat[0]+static_cast<double>(tauy)*lat[1]+static_cast<double>(tauz)*lat[2];
						r2 = (xyz[iat]-xyz[i]+tau).norm2();
						if(r2 > para.cn_thr2) continue;
						rr = (para.rcov[iz[i]]+para.rcov[iz[iat]])/sqrt(r2);
						xn += 1.0/(1.0+exp(-para.k1*(rr-1.0)));
					}
		cn[i] = xn;
	}
}

//void Vdwd3::pbcthreebody(vector<double> &cc6ab, double &e63)
//{

//}

void Vdwd3::cal_energy()
{
	TITLE("Vdwd3","cal_energy");
	init(ucell);

	int ij;	
	double c6 = 0.0, c8 = 0.0, r2 = 0.0, r6 = 0.0, r8 = 0.0, rr = 0.0, damp6 = 0.0, damp8 = 0.0;
	double e6 = 0.0, e8 = 0.0, e63 = 0.0;
	vector<double> cc6ab(ucell.nat*ucell.nat), cn(ucell.nat);
	pbcncoord(cn);
	Vector3<double> tau;
	if(para.version == "d3_0") // DFT-D3(zero-damping)
	{
		double tmp;
		for(int iat=0; iat!=ucell.nat-1; iat++)
			for(int jat=iat+1; jat!=ucell.nat; jat++)
			{
				getc6(iz[iat], iz[jat], cn[iat], cn[jat], c6);
				if(para.abc)
				{
					ij = lin(iat, jat);
					cc6ab[ij] = sqrt(c6);
				}
				for(int taux=-rep_vdw[0]; taux<=rep_vdw[0]; taux++)
					for(int tauy=-rep_vdw[1]; tauy<=rep_vdw[1]; tauy++)
						for(int tauz=-rep_vdw[2]; tauz<=rep_vdw[2]; tauz++)
						{
							tau = static_cast<double>(taux)*lat[0]+static_cast<double>(tauy)*lat[1]+static_cast<double>(tauz)*lat[2];
							r2 = (xyz[iat]-xyz[jat]+tau).norm2();
							if(r2 > para.rthr2) continue;
							rr = para.r0ab[iz[iat]][iz[jat]]/sqrt(r2);
							// zero-damping function
							tmp = para.rs6*rr;
							damp6 = 1.0/(1.0+6.0*pow(tmp, para.alp6));
							tmp = para.rs18*rr;
							damp8 = 1.0/(1.0+6.0*pow(tmp, para.alp8));

							r6 = pow(r2, 3);
							e6 += damp6/r6*c6;

							c8 = 3.0*para.r2r4[iz[iat]]*para.r2r4[iz[jat]]*c6;
							r8 = r6*r2;
							e8 += c8*damp8/r8;
						} // end tau
			} // end jat

		for(int iat=0; iat!=ucell.nat; iat++)
		{
			int jat = iat;
			getc6(iz[iat], iz[jat], cn[iat], cn[jat], c6);
			if(para.abc)
			{
				ij = lin(iat, jat);
				cc6ab[ij] = sqrt(c6);
			}
			for(int taux=-rep_vdw[0]; taux<=rep_vdw[0]; taux++)
				for(int tauy=-rep_vdw[1]; tauy<=rep_vdw[1]; tauy++)
					for(int tauz=-rep_vdw[2]; tauz<=rep_vdw[2]; tauz++)
					{
						if(taux==0 && tauy==0 && tauz==0) continue;
						tau = static_cast<double>(taux)*lat[0]+static_cast<double>(tauy)*lat[1]+static_cast<double>(tauz)*lat[2];
						r2 = tau.norm2();
						if(r2 > para.rthr2) continue;
						rr = para.r0ab[iz[iat]][iz[jat]]/sqrt(r2);

						// zero-damping function
						tmp = para.rs6*rr;
						damp6 = 1.0/(1.0+6.0*pow(tmp,para.alp6));
						tmp = para.rs18*rr;
						damp8 = 1.0/(1.0+6.0*pow(tmp, para.alp8));

						r6 = pow(r2, 3);
						e6 += damp6/r6*c6*0.5;

						c8 = 3.0*para.r2r4[iz[iat]]*para.r2r4[iz[jat]]*c6;
						r8 = r6*r2;
						e8 += c8*damp8/r8*0.5;
					} // end tau
		} // end iat
	} // end d3_0
	else if(para.version == "d3_bj") // DFT-D3(BJ-damping)
	{
		double r42;
		for(int iat=0; iat!=ucell.nat; iat++)
		{
			for(int jat=iat+1; jat!=ucell.nat; jat++)
			{
				getc6(iz[iat], iz[jat], cn[iat], cn[jat], c6);
				if(para.abc)
				{
					ij = lin(iat, jat);
					cc6ab[ij] = sqrt(c6);
				}
				for(int taux=-rep_vdw[0]; taux<=rep_vdw[0]; taux++)
					for(int tauy=-rep_vdw[1]; tauy<=rep_vdw[1]; tauy++)
						for(int tauz=-rep_vdw[2]; tauz<=rep_vdw[2]; tauz++)
						{
							tau = static_cast<double>(taux)*lat[0]+static_cast<double>(tauy)*lat[1]+static_cast<double>(tauz)*lat[2];
							r2 = (xyz[iat]-xyz[jat]+tau).norm2();
							if(r2 > para.rthr2) continue;
							rr = para.r0ab[iz[iat]][iz[jat]]/sqrt(r2);

							// BJ-damping function
							r42 = para.r2r4[iz[iat]]*para.r2r4[iz[jat]];
							damp6 = pow((para.rs6*sqrt(3.0*r42)+para.rs18), 6);
							damp8 = pow((para.rs6*sqrt(3.0*r42)+para.rs18), 8);

							r6 = pow(r2, 3);
							e6 += c6/(r6+damp6);

							c8 = 3.0*c6*r42;
							r8 = r6*r2;
							e8 += c8/(r8+damp8);
						} // end tau
			} // end jat
			int jat = iat;
			getc6(iz[iat], iz[jat], cn[iat], cn[jat], c6);
			r42 = para.r2r4[iz[iat]]*para.r2r4[iz[jat]];
			damp6 = pow((para.rs6*sqrt(3.0*r42)+para.rs18), 6);
			damp8 = pow((para.rs6*sqrt(3.0*r42)+para.rs18), 8);
			if(para.abc)
			{
				ij = lin(iat, jat);
				cc6ab[ij] = sqrt(c6);
			}
			for(int taux=-rep_vdw[0]; taux<=rep_vdw[0]; taux++)
				for(int tauy=-rep_vdw[1]; tauy<=rep_vdw[1]; tauy++)
					for(int tauz=-rep_vdw[2]; tauz<=rep_vdw[2]; tauz++)
					{
						if(taux==0 && tauy==0 && tauz==0) continue;
						tau = static_cast<double>(taux)*lat[0]+static_cast<double>(tauy)*lat[1]+static_cast<double>(tauz)*lat[2];
						r2 = tau.norm2();
						if(r2 > para.rthr2) continue;
						rr = para.r0ab[iz[iat]][iz[jat]]/sqrt(r2);

						r6 = pow(r2, 3);
						e6 += c6/(r6+damp6)*0.5;

						c8 = 3.0*c6*r42;
						r8 = r6*r2;
						e8 += c8/(r8+damp8)*0.5;
					} // end tau
		} // end iat
	} // end d3_bj

	//if(para.abc)
	//{
	//	pbcthreebody(cc6ab, e63)
	//}
	energy = (-para.s6*e6-para.s18*e8-e63)*2;
}

void Vdwd3::get_dc6_dcnij(int &mxci, int &mxcj, double &cni, double &cnj, int &izi, int &izj, int &iat, int &jat, double &c6check, double &dc6i, double &dc6j)
{
	double r_save = 9999.0, c6mem = -1e99, zaehler = 0.0, nenner = 0.0;
	double dzaehler_i = 0.0, dnenner_i = 0.0, dzaehler_j = 0.0, dnenner_j = 0.0;
	double c6ref = 0.0, cn_refi = 0.0, cn_refj = 0.0, r = 0.0, expterm = 0.0, term = 0.0;
	for(size_t a=0; a!=mxci; a++)
		for(size_t b=0; b!=mxcj; b++)
		{
			c6ref = para.c6ab[0][b][a][izj][izi];
			if(c6ref > 0)
			{
				cn_refi = para.c6ab[1][b][a][izj][izi];
				cn_refj = para.c6ab[2][b][a][izj][izi];
				r = (cn_refi-cni)*(cn_refi-cni)+(cn_refj-cnj)*(cn_refj-cnj);
				if(r < r_save)
				{
					r_save = r;
					c6mem = c6ref;
				}
				expterm = exp(para.k3*r);
				zaehler += c6ref*expterm;
				nenner += expterm;
				expterm *= 2.0*para.k3;
				term = expterm*(cni-cn_refi);
				dzaehler_i += c6ref*term;
				dnenner_i += term;

				term = expterm*(cnj-cn_refj);
				dzaehler_j += c6ref*term;
				dnenner_j += term;
			}
		}

	if(nenner > 1e-99)
	{
		c6check = zaehler/nenner;
		dc6i = ((dzaehler_i*nenner)-(dnenner_i*zaehler))/(nenner*nenner);
		dc6j = ((dzaehler_j*nenner)-(dnenner_j*zaehler))/(nenner*nenner);
	}
	else
	{
		c6check = c6mem;
		dc6i = 0.0;
		dc6j = 0.0;
	}
}

void Vdwd3::pbcgdisp(vector<Vector3<double>> &g, matrix &sigma)
{
	init(ucell);
	vector<double> c6save(ucell.nat*(ucell.nat+1)), dc6_rest_sum(ucell.nat*(ucell.nat+1)/2), dc6i(ucell.nat), cn(ucell.nat);
	pbcncoord(cn);
	vector<vector<double>> dc6ij(ucell.nat, vector<double>(ucell.nat));
	double c6 = 0.0, dc6iji = 0.0, dc6ijj = 0.0;
	double r = 0.0, r0 = 0.0, r2 = 0.0, r6 = 0.0, r7 = 0.0, r8 = 0.0, r9 = 0.0;
	double r42 = 0.0, rcovij = 0.0, t6 = 0.0, t8 = 0.0, dc6_rest = 0.0;
	int linii = 0, linij = 0;
	Vector3<double> tau;
	vector<vector<vector<vector<double>>>> drij(ucell.nat*(ucell.nat+1)/2, vector<vector<vector<double>>>(2*rep_vdw[0]+1, vector<vector<double>>(2*rep_vdw[1]+1, vector<double>(2*rep_vdw[2]+1))));
	if(para.version == "d3_0")
	{
		double damp6 = 0.0, damp8 = 0.0;
		for(int iat=0; iat!=ucell.nat; iat++)
		{
			get_dc6_dcnij(para.mxc[iz[iat]], para.mxc[iz[iat]], cn[iat], cn[iat], iz[iat], iz[iat], iat, iat, c6, dc6iji, dc6ijj);

			linii = lin(iat, iat);
			c6save[linii] = c6;
			dc6ij[iat][iat] = dc6iji;
			r0 = para.r0ab[iz[iat]][iz[iat]];
			r42 = para.r2r4[iz[iat]]*para.r2r4[iz[iat]];
			rcovij = para.rcov[iz[iat]]+para.rcov[iz[iat]];

			for(int taux=-rep_vdw[0]; taux<=rep_vdw[0]; taux++)
				for(int tauy=-rep_vdw[1]; tauy<=rep_vdw[1]; tauy++)
					for(int tauz=-rep_vdw[2]; tauz<=rep_vdw[2]; tauz++)
					{
						tau = static_cast<double>(taux)*lat[0]+static_cast<double>(tauy)*lat[1]+static_cast<double>(tauz)*lat[2];

						// first dE/d(tau) 
						r2 = tau.norm2();
						if(r2 > 0.1 && r2 < para.rthr2)
						{
							r = sqrt(r2);
							r6 = pow(r2, 3);
							r7 = r6*r;
							r8 = r6*r2;
							r9 = r8*r;

							t6 = pow(r/(para.rs6*r0), -para.alp6);
							damp6 = 1.0/(1.0+6.0*t6);
							t8 = pow(r/(para.rs18*r0), -para.alp8);
							damp8 = 1.0/(1.0+6.0*t8);

							// d(r^(-6))/d(tau)
							drij[linii][taux+rep_vdw[0]][tauy+rep_vdw[1]][tauz+rep_vdw[2]] += (-para.s6*(6.0/(r7)*c6*damp6)-para.s18*(24.0/(r9)*c6*r42*damp8))*0.5;
							// d(f_dmp)/d(tau)
							drij[linii][taux+rep_vdw[0]][tauy+rep_vdw[1]][tauz+rep_vdw[2]] += (para.s6*c6/r7*6.0*para.alp6*t6*damp6*damp6+para.s18*c6*r42/r9*18.0*para.alp8*t8*damp8*damp8)*0.5;

							dc6_rest = (para.s6/r6*damp6+3.0*para.s18*r42/r8*damp8)*0.5;
							dc6i[iat] += dc6_rest*(dc6iji+dc6ijj);
							dc6_rest_sum[linii] += dc6_rest;
						}
					} // end tau
			for(int jat=0; jat!=iat; jat++)
			{
				get_dc6_dcnij(para.mxc[iz[iat]], para.mxc[iz[jat]], cn[iat], cn[jat], iz[iat], iz[jat], iat, jat, c6, dc6iji, dc6ijj);

				linij = lin(iat, jat);
				c6save[linij] = c6;
				r0 = para.r0ab[iz[iat]][iz[jat]];
				r42 = para.r2r4[iz[iat]]*para.r2r4[iz[jat]];
				rcovij = para.rcov[iz[iat]]+para.rcov[iz[jat]];
				dc6ij[jat][iat] = dc6iji;
				dc6ij[iat][jat] = dc6ijj;
				for(int taux=-rep_vdw[0]; taux<=rep_vdw[0]; taux++)
					for(int tauy=-rep_vdw[1]; tauy<=rep_vdw[1]; tauy++)
						for(int tauz=-rep_vdw[2]; tauz<=rep_vdw[2]; tauz++)
						{
							tau = static_cast<double>(taux)*lat[0]+static_cast<double>(tauy)*lat[1]+static_cast<double>(tauz)*lat[2];
							r2 = (xyz[jat]-xyz[iat]+tau).norm2();
							if(r2 > para.rthr2) continue;

							r = sqrt(r2);
							r6 = pow(r2, 3);
							r7 = r6*r;
							r8 = r6*r2;
							r9 = r8*r;

							t6 = pow(r/(para.rs6*r0), -para.alp6);
							damp6 = 1.0/(1.0+6.0*t6);
							t8 = pow(r/(para.rs18*r0), -para.alp8);
							damp8 = 1.0/(1.0+6.0*t8);

							// d(r^(-6))/d(r_ij)
							drij[linij][taux+rep_vdw[0]][tauy+rep_vdw[1]][tauz+rep_vdw[2]] += -para.s6*(6.0/(r7)*c6*damp6)-para.s18*(24.0/(r9)*c6*r42*damp8);
							// d(f_dmp)/d(r_ij)
							drij[linij][taux+rep_vdw[0]][tauy+rep_vdw[1]][tauz+rep_vdw[2]] += para.s6*c6/r7*6.0*para.alp6*t6*damp6*damp6+para.s18*c6*r42/r9*18.0*para.alp8*t8*damp8*damp8;

							dc6_rest = para.s6/r6*damp6+3.0*para.s18*r42/r8*damp8;
							dc6i[iat] += dc6_rest*dc6iji;
							dc6i[jat] += dc6_rest*dc6ijj;
							dc6_rest_sum[linij] += dc6_rest;
						} // end tau
			} // end jat
		} // end iat
	} // end d3_0
	else if(para.version == "d3_bj")
	{
		double r4;
		for(int iat=0; iat!=ucell.nat; iat++)
		{
			get_dc6_dcnij(para.mxc[iz[iat]], para.mxc[iz[iat]], cn[iat], cn[iat], iz[iat], iz[iat], iat, iat, c6, dc6iji, dc6ijj);

			linii = lin(iat, iat);
			c6save[linii] = c6;
			dc6ij[iat][iat] = dc6iji;
			r42 = para.r2r4[iz[iat]]*para.r2r4[iz[iat]];
			r0 = para.rs6*sqrt(3.0*r42)+para.rs18;
			rcovij = para.rcov[iz[iat]]+para.rcov[iz[iat]];

			for(int taux=-rep_vdw[0]; taux<=rep_vdw[0]; taux++)
				for(int tauy=-rep_vdw[1]; tauy<=rep_vdw[1]; tauy++)
					for(int tauz=-rep_vdw[2]; tauz<=rep_vdw[2]; tauz++)
					{
						tau = static_cast<double>(taux)*lat[0]+static_cast<double>(tauy)*lat[1]+static_cast<double>(tauz)*lat[2];

						// first dE/d(tau) 
						r2 = tau.norm2();
						if(r2 > 0.1 && r2 < para.rthr2)
						{
							r = sqrt(r2);
							r4 = r2*r2;
							r6 = pow(r2, 3);
							r7 = r6*r;
							r8 = r6*r2;
							r9 = r8*r;

							t6 = r6+pow(r0, 6);
							t8 = r8+pow(r0, 8);

							// d(1/r^(-6)+r0^6)/d(r)
							drij[linii][taux+rep_vdw[0]][tauy+rep_vdw[1]][tauz+rep_vdw[2]] += -para.s6*c6*6.0*r4*r/(t6*t6)*0.5-para.s18*c6*24.0*r42*r7/(t8*t8)*0.5;

							dc6_rest = (para.s6/t6+3.0*para.s18*r42/t8)*0.5;
							dc6i[iat] += dc6_rest*(dc6iji+dc6ijj);
							dc6_rest_sum[linii] += dc6_rest;
						}
					} // end tau
			for(int jat=0; jat!=iat; jat++)
			{
				get_dc6_dcnij(para.mxc[iz[iat]], para.mxc[iz[jat]], cn[iat], cn[jat], iz[iat], iz[jat], iat, jat, c6, dc6iji, dc6ijj);

				linij = lin(iat, jat);
				c6save[linij] = c6;
				r42 = para.r2r4[iz[iat]]*para.r2r4[iz[jat]];
				r0 = para.rs6*sqrt(3.0*r42)+para.rs18;
				rcovij = para.rcov[iz[iat]]+para.rcov[iz[jat]];
				dc6ij[jat][iat] = dc6iji;
				dc6ij[iat][jat] = dc6ijj;
				for(int taux=-rep_vdw[0]; taux<=rep_vdw[0]; taux++)
					for(int tauy=-rep_vdw[1]; tauy<=rep_vdw[1]; tauy++)
						for(int tauz=-rep_vdw[2]; tauz<=rep_vdw[2]; tauz++)
						{
							tau = static_cast<double>(taux)*lat[0]+static_cast<double>(tauy)*lat[1]+static_cast<double>(tauz)*lat[2];
							r2 = (xyz[jat]-xyz[iat]+tau).norm2();
							if(r2 > para.rthr2) continue;

							r = sqrt(r2);
							r4 = r2*r2;
							r6 = pow(r2, 3);
							r7 = r6*r;
							r8 = r6*r2;
							r9 = r8*r;

							t6 = r6+pow(r0, 6);
							t8 = r8+pow(r0, 8);

							drij[linij][taux+rep_vdw[0]][tauy+rep_vdw[1]][tauz+rep_vdw[2]] += -para.s6*c6*6.0*r4*r/(t6*t6)-para.s18*c6*24.0*r42*r7/(t8*t8);

							dc6_rest = para.s6/t6+3.0*para.s18*r42/t8;
							dc6i[iat] += dc6_rest*dc6iji;
							dc6i[jat] += dc6_rest*dc6ijj;
							dc6_rest_sum[linij] += dc6_rest;
						} // end tau
			} // end jat
		} // end iat
	} // end d3_bj 

	//if(para.abc)
	//{
	//	
	//}

	// dE/dr_ij * dr_ij/dxyz_i
	double expterm, dcnn, x1;
	Vector3<double> rij, vec3;
	for(int iat=1; iat!=ucell.nat; iat++)
		for(int jat=0; jat!=iat; jat++)
		{
			linij = lin(iat, jat);
			rcovij = para.rcov[iz[iat]]+para.rcov[iz[jat]];
			for(int taux=-rep_vdw[0]; taux<=rep_vdw[0]; taux++)
				for(int tauy=-rep_vdw[1]; tauy<=rep_vdw[1]; tauy++)
					for(int tauz=-rep_vdw[2]; tauz<=rep_vdw[2]; tauz++)
					{
						tau = static_cast<double>(taux)*lat[0]+static_cast<double>(tauy)*lat[1]+static_cast<double>(tauz)*lat[2];
						rij = xyz[jat]-xyz[iat]+tau;
						r2 = rij.norm2();
						if(r2 > para.rthr2 || r2 < 0.5) continue;
						r = sqrt(r2);
						if(r2 < para.cn_thr2)
						{
							expterm = exp(-para.k1*(rcovij/r-1.0));
							dcnn = -para.k1*rcovij*expterm/(r2*(expterm+1.0)*(expterm+1.0));
						}
						else
							dcnn = 0.0;
						x1 = drij[linij][taux+rep_vdw[0]][tauy+rep_vdw[1]][tauz+rep_vdw[2]]+dcnn*(dc6i[iat]+dc6i[jat]);
						vec3 = x1*rij/r;
						g[iat] += vec3;
						g[jat] -= vec3;

						vector<double> vec = {vec3.x, vec3.y, vec3.z};
						vector<double> rij_vec = {rij.x, rij.y, rij.z};
						for(size_t i=0; i!=3; i++)
							for(size_t j=0; j!=3; j++)
							{
								sigma(i, j) += vec[j]*rij_vec[i];
							}
					} // end tau
		} // end iat, jat
	for(int iat=0; iat!=ucell.nat; iat++)
	{
		linii = lin(iat, iat);
		rcovij = para.rcov[iz[iat]]+para.rcov[iz[iat]];
		for(int taux=-rep_vdw[0]; taux<=rep_vdw[0]; taux++)
			for(int tauy=-rep_vdw[1]; tauy<=rep_vdw[1]; tauy++)
				for(int tauz=-rep_vdw[2]; tauz<=rep_vdw[2]; tauz++)
				{
					if(taux == 0 && tauy == 0 && tauz == 0) continue;
					tau = static_cast<double>(taux)*lat[0]+static_cast<double>(tauy)*lat[1]+static_cast<double>(tauz)*lat[2];
					r2 = tau.norm2();
					r = sqrt(r2);
					if(r2 < para.cn_thr2)
					{
						expterm = exp(-para.k1*(rcovij/r-1.0));
						dcnn = -para.k1*rcovij*expterm/(r2*(expterm+1.0)*(expterm+1.0));
					}
					else
						dcnn = 0.0;
					x1 = drij[linii][taux+rep_vdw[0]][tauy+rep_vdw[1]][tauz+rep_vdw[2]]+dcnn*dc6i[iat];

					vec3 = x1*tau/r;
					vector<double> vec = {vec3.x, vec3.y, vec3.z};
					vector<double> tau_vec = {tau.x, tau.y, tau.z};
					for(size_t i=0; i!=3; i++)
						for(size_t j=0; j!=3; j++)
						{
							sigma(i, j) += vec[j]*tau_vec[i];
						}
				} // end tau
	} // end iat
}

void Vdwd3::cal_force()
{
	TITLE("Vdwd3","cal_force");

	force.clear();
	force.resize(ucell.nat);

	vector<Vector3<double>> g;
	g.clear();
	g.resize(ucell.nat);
	matrix sigma(3, 3);

	pbcgdisp(g, sigma);

	for(size_t iat=0; iat!=ucell.nat; iat++)
		force[iat] = -2.0*g[iat];

}

void Vdwd3::cal_stress()
{
	TITLE("Vdwd3","cal_stress");

	vector<Vector3<double>> g;
	g.clear();
	g.resize(ucell.nat);
	matrix sigma(3, 3);

	pbcgdisp(g, sigma);
	
	stress = Matrix3(
		2.0*sigma(0, 0), 2.0*sigma(0, 1), 2.0*sigma(0, 2),
		2.0*sigma(1, 0), 2.0*sigma(1, 1), 2.0*sigma(1, 2),
		2.0*sigma(2, 0), 2.0*sigma(2, 1), 2.0*sigma(2, 2)
	)/ucell.omega;

}