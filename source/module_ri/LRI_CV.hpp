//=======================
// AUTHOR : Peize Lin
// DATE :   2022-08-17
//=======================

#ifndef LRI_CV_HPP
#define LRI_CV_HPP

#include "LRI_CV.h"
#include "src_ri/exx_abfs-abfs_index.h"
#include "Inverse_Matrix.h"
#include "RI_Util.h"
#include "module_base/tool_title.h"
#include "module_base/timer.h"
#include <omp.h>

template<typename Tdata>
LRI_CV<Tdata>::LRI_CV()
{
	pthread_rwlock_init(&rwlock_Vw,NULL);
}

template<typename Tdata>
LRI_CV<Tdata>::~LRI_CV()
{
	pthread_rwlock_destroy(&rwlock_Vw);
}


template<typename Tdata>
void LRI_CV<Tdata>::set_orbitals(
	const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &lcaos_in,
	const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &abfs_in,
	const std::vector<std::vector<std::vector<Numerical_Orbital_Lm>>> &abfs_ccp_in,
	const double &kmesh_times,
	const double &ccp_rmesh_times_in)
{
	ModuleBase::TITLE("LRI_CV", "set_orbitals");
	ModuleBase::timer::tick("LRI_CV", "set_orbitals");

	this->lcaos = lcaos_in;
	this->abfs = abfs_in;
	this->abfs_ccp = abfs_ccp_in;
	this->ccp_rmesh_times = ccp_rmesh_times_in;

	const ModuleBase::Element_Basis_Index::Range
		range_lcaos = Exx_Abfs::Abfs_Index::construct_range( lcaos );
	this->index_lcaos = ModuleBase::Element_Basis_Index::construct_index( range_lcaos );

	const ModuleBase::Element_Basis_Index::Range
		range_abfs = Exx_Abfs::Abfs_Index::construct_range( abfs );
	this->index_abfs = ModuleBase::Element_Basis_Index::construct_index( range_abfs );		

	this->m_abfs_abfs.init( 2, kmesh_times, (1+this->ccp_rmesh_times)/2.0 );
	this->m_abfs_abfs.init_radial( this->abfs_ccp, this->abfs );
	this->m_abfs_abfs.init_radial_table();
	
	this->m_abfslcaos_lcaos.init( 1, kmesh_times, 1 );
	this->m_abfslcaos_lcaos.init_radial( this->abfs_ccp, this->lcaos, this->lcaos );	
	this->m_abfslcaos_lcaos.init_radial_table();

	ModuleBase::timer::tick("LRI_CV", "set_orbitals");
}



template<typename Tdata>
auto LRI_CV<Tdata>::cal_datas(
	const std::vector<TA> &list_A0,
	const std::vector<TAC> &list_A1,
	const double threshold,
	const bool flag_writable,
	const T_func_DPcal_data &func_DPcal_data)
-> std::map<TA,std::map<TAC,Tensor<Tdata>>>
{
	ModuleBase::TITLE("LRI_CV","cal_datas");
	ModuleBase::timer::tick("LRI_CV", "cal_datas");

	std::map<TA,std::map<TAC,Tensor<Tdata>>> Datas;
	#pragma omp parallel for
	for(size_t i0=0; i0<list_A0.size(); ++i0)
	{
		#pragma omp parallel for
		for(size_t i1=0; i1<list_A1.size(); ++i1)
		{
			const TA iat0 = list_A0[i0];
			const TA iat1 = list_A1[i1].first;
			const TC &cell1 = list_A1[i1].second;
			const int it0 = GlobalC::ucell.iat2it[iat0];
			const int ia0 = GlobalC::ucell.iat2ia[iat0];
			const int it1 = GlobalC::ucell.iat2it[iat1];
			const int ia1 = GlobalC::ucell.iat2ia[iat1];
			const ModuleBase::Vector3<double> tau1 = GlobalC::ucell.atoms[it0].tau[ia0];
			const ModuleBase::Vector3<double> tau2 = GlobalC::ucell.atoms[it1].tau[ia1];
			const double Rcut = std::min(
				GlobalC::ORB.Phi[it0].getRcut() * this->ccp_rmesh_times + GlobalC::ORB.Phi[it1].getRcut(),
				GlobalC::ORB.Phi[it1].getRcut() * this->ccp_rmesh_times + GlobalC::ORB.Phi[it0].getRcut());
			const Abfs::Vector3_Order<double> R_delta = -tau1+tau2+(RI_Util::array3_to_Vector3(cell1)*GlobalC::ucell.latvec);
			if( R_delta.norm()*GlobalC::ucell.lat0 < Rcut )
			{
				const Tensor<Tdata> Data = func_DPcal_data(it0, it1, R_delta, flag_writable);
				if(Data.norm(std::numeric_limits<double>::max()) > threshold)
				{
					#pragma omp critical(LRI_CV_cal_datas)
					Datas[list_A0[i0]][list_A1[i1]] = Data;
				}
			}
		}
	}
	ModuleBase::timer::tick("LRI_CV", "cal_datas");
	return std::move(Datas);
}


template<typename Tdata>
auto LRI_CV<Tdata>::cal_Vs(
	const std::vector<TA> &list_A0,
	const std::vector<TAC> &list_A1,
	const double threshold_V,
	const bool flag_writable)
-> std::map<TA,std::map<TAC,Tensor<Tdata>>>
{
	ModuleBase::TITLE("LRI_CV","cal_Vs");
	const T_func_DPcal_data func_DPcal_V = std::bind( &LRI_CV<Tdata>::DPcal_V, this,
		std::placeholders::_1,
		std::placeholders::_2,
		std::placeholders::_3,
		std::placeholders::_4);
	return this->cal_datas(list_A0, list_A1, threshold_V, flag_writable, func_DPcal_V);
}

template<typename Tdata>
auto LRI_CV<Tdata>::cal_Cs(
	const std::vector<TA> &list_A0,
	const std::vector<TAC> &list_A1,
	const double threshold_C,
	const bool flag_writable)
-> std::map<TA,std::map<TAC,Tensor<Tdata>>>
{
	ModuleBase::TITLE("LRI_CV","cal_Cs");
	const T_func_DPcal_data func_DPcal_C = std::bind( &LRI_CV<Tdata>::DPcal_C, this,
		std::placeholders::_1,
		std::placeholders::_2,
		std::placeholders::_3,
		std::placeholders::_4);
	return this->cal_datas(list_A0, list_A1, threshold_C, flag_writable, func_DPcal_C);
}



template<typename Tdata>
Tensor<Tdata>
LRI_CV<Tdata>::DPcal_V(
	const int it0,
	const int it1,
	const Abfs::Vector3_Order<double> &R,
	const bool flag_writable)
{
	pthread_rwlock_rdlock(&this->rwlock_Vw);
	Tensor<Tdata> V = this->Vws[it0][it1][R];
	Tensor<Tdata> VT = this->Vws[it1][it0][-R];
	pthread_rwlock_unlock(&this->rwlock_Vw);

	if(!V && !VT)
	{
		V = this->m_abfs_abfs.cal_overlap_matrix<Tdata>(
			it0, it1, {0,0,0}, R,
			this->index_abfs, this->index_abfs,
			Matrix_Orbs11::Matrix_Order::AB);
		if(flag_writable)
		{
			VT = V.transpose();
			pthread_rwlock_wrlock(&this->rwlock_Vw);
			this->Vws[it0][it1][R] = V;
			this->Vws[it1][it0][-R] = VT;
			pthread_rwlock_unlock(&this->rwlock_Vw);
		}
	}
	else if(V && !VT)
	{	
		if(flag_writable)
		{
			VT = V.transpose();
			pthread_rwlock_wrlock(&this->rwlock_Vw);
			this->Vws[it1][it0][-R] = VT;
			pthread_rwlock_unlock(&this->rwlock_Vw);
		}
	}
	else if(!V && VT)
	{
		V = VT.transpose();
		if(flag_writable)
		{		
			pthread_rwlock_wrlock(&this->rwlock_Vw);
			this->Vws[it0][it1][R] = V;
			pthread_rwlock_unlock(&this->rwlock_Vw);
		}
	}
	return V;
}



template<typename Tdata>
Tensor<Tdata>
LRI_CV<Tdata>::DPcal_C(
	const int it0,
	const int it1,
	const Abfs::Vector3_Order<double> &R,
	const bool flag_writable)
{
	auto transpose12 = [](const Tensor<Tdata> &c_in) -> Tensor<Tdata>
	{
		Tensor<Tdata> c_out({c_in.shape[0], c_in.shape[2], c_in.shape[1]});
		for(size_t i0=0; i0<c_in.shape[0]; ++i0)
			for(size_t i1=0; i1<c_in.shape[1]; ++i1)
				for(size_t i2=0; i2<c_in.shape[2]; ++i2)
					c_out(i0,i2,i1) = c_in(i0,i1,i2);
		return c_out;
	};

	pthread_rwlock_rdlock(&this->rwlock_Cw);
	const std::array<Tensor<Tdata>, 2>
		Cr = {this->Cws[it0][it1][R],
		      this->Cws[it1][it0][-R]};
	pthread_rwlock_unlock(&this->rwlock_Cw);

	if(Cr[0])
	{
		return Cr[0].data->size() ? Cr[0] : Tensor<Tdata>({});
	}
	else
	{
		if( (ModuleBase::Vector3<double>(0,0,0)==R) && (it0==it1) )
		{
			const Tensor<Tdata> A = 
				this->m_abfslcaos_lcaos.cal_overlap_matrix<Tdata>(
					it0, it1, {0,0,0}, {0,0,0},
					this->index_abfs, this->index_lcaos, this->index_lcaos,
					Matrix_Orbs21::Matrix_Order::A1A2B);
			const size_t sa=A.shape[0], sl0=A.shape[1], sl1=A.shape[2];
			const Tensor<Tdata> V = DPcal_V( it0, it0, {0,0,0}, true);
			const Tensor<Tdata> L = cal_I(V);
			Tensor<Tdata> C = Global_Func::convert<Tdata>(0.5) * (L * A.reshape({sa, sl0*sl1})).reshape({sa,sl0,sl1});					// Attention 0.5!
			if(flag_writable)
			{
				pthread_rwlock_wrlock(&this->rwlock_Cw);
				this->Cws[it0][it1][{0,0,0}] = C;
				pthread_rwlock_unlock(&this->rwlock_Cw);
			}
			return C;
		}
		else
		{
			const std::array<Tensor<Tdata>, 2>
				 A = {this->m_abfslcaos_lcaos.cal_overlap_matrix<Tdata>(
							it0, it1, {0,0,0}, R,
							this->index_abfs, this->index_lcaos, this->index_lcaos,
							Matrix_Orbs21::Matrix_Order::A1A2B),
				      this->m_abfslcaos_lcaos.cal_overlap_matrix<Tdata>(
							it1, it0, {0,0,0}, -R,
							this->index_abfs, this->index_lcaos, this->index_lcaos,
							Matrix_Orbs21::Matrix_Order::A1BA2)};
			const size_t sa0=A[0].shape[0], sa1=A[1].shape[0], sl0=A[0].shape[1], sl1=A[0].shape[2];

			const std::vector<std::vector<Tensor<Tdata>>>
				V = {{DPcal_V(it0, it0, {0,0,0}, true),
				      DPcal_V(it0, it1, R,       false)},
				     {DPcal_V(it1, it0, -R,      false),
				      DPcal_V(it1, it1, {0,0,0}, true)}};

			const std::vector<std::vector<Tensor<Tdata>>>
				L = cal_I(V);

			std::array<Tensor<Tdata>, 2>
				C = {( L[0][0]*A[0].reshape({sa0,sl0*sl1}) + L[0][1]*A[1].reshape({sa1,sl0*sl1}) ).reshape({sa0,sl0,sl1}),
				     ( L[1][0]*A[0].reshape({sa0,sl0*sl1}) + L[1][1]*A[1].reshape({sa1,sl0*sl1}) ).reshape({sa1,sl0,sl1})};
			if(flag_writable)
			{
				pthread_rwlock_wrlock(&this->rwlock_Cw);
				this->Cws[it0][it1][R] = C[0];
				this->Cws[it1][it0][-R] = transpose12(C[1]);
				pthread_rwlock_unlock(&this->rwlock_Cw);
			}
			return C[0];
		}
	} // end if(!Cr[0])
}




template<typename Tdata>
Tensor<Tdata>
LRI_CV<Tdata>::cal_I( const Tensor<Tdata> &m )
{
	Inverse_Matrix<Tdata> I;
	I.input(m);
	I.cal_inverse( Inverse_Matrix<Tdata>::Method::potrf );
	return I.output();
}

template<typename Tdata>
std::vector<std::vector<Tensor<Tdata>>>
LRI_CV<Tdata>::cal_I( const std::vector<std::vector<Tensor<Tdata>>> &ms )
{
	Inverse_Matrix<Tdata> I;
	I.input(ms);
	I.cal_inverse( Inverse_Matrix<Tdata>::Method::potrf );
	return I.output({ms[0][0].shape[0], ms[1][0].shape[0]}, {ms[0][0].shape[1], ms[0][1].shape[1]});
}

#endif