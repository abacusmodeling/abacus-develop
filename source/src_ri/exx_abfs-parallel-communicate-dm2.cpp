#include "exx_abfs-parallel-communicate-dm2.h"
#include "../src_pw/global.h"


void Exx_Abfs::Parallel::Communicate::DM2::init( 
	const set<pair<size_t,size_t>> &H_atom_pairs_core,
	const double threshold_in)
{
	threshold = threshold_in;
	
	atom_in_exx.row.resize(ucell.nat,false);
	atom_in_exx.col.resize(ucell.nat,false);
	atom_in_exx.row_col.resize(ucell.nat,vector<bool>(ucell.nat,false));
	
	for( const auto &pair : H_atom_pairs_core )
	{
		atom_in_exx.row[pair.first] = true;
		atom_in_exx.col[pair.second] = true;
		atom_in_exx.row_col[pair.first][pair.second] = true;
	}
}

void Exx_Abfs::Parallel::Communicate::DM2::clear_DMr()
{
	DMr.clear();
	DMr.resize(GlobalV::NSPIN);
}

void Exx_Abfs::Parallel::Communicate::DM2::set_DM_gamma( const matrix &DM_2D, const int is, const pair<int,int> &index_begin )
{
	auto get_iats = []( const int n, const int i_begin, const vector<bool> &in_exx ) -> map<int,pair<int,int>>
	{
		map<int,pair<int,int>> iats;
		for( int i=0; i<n; ++i )
		{
			const int iat = ucell.iwt2iat[i+i_begin];
			if( in_exx[iat] )
				iats[iat].second=i;
		}
		for( int i=n-1; i>=0; --i )
		{
			const int iat = ucell.iwt2iat[i+i_begin];
			if( in_exx[iat] )
				iats[iat].first=i;
		}
		return iats;
	};
	const map<int,pair<int,int>> iat1s = get_iats( DM_2D.nr, index_begin.first, atom_in_exx.row );
	const map<int,pair<int,int>> iat2s = get_iats( DM_2D.nc, index_begin.second, atom_in_exx.col );
	
	const double SPIN_multiple = 0.5*GlobalV::NSPIN;	
	for( const auto &iat1sA : iat1s )
	{
		const int iat1 = iat1sA.first;
		const auto &iat1_range = iat1sA.second;
		for( const auto &iat2sA : iat2s )
		{
			const int iat2 = iat2sA.first;
			const auto &iat2_range = iat2sA.second;
			if( !atom_in_exx.row_col[iat1][iat2] )	continue;
			
			matrix DM_local( ucell.atoms[ucell.iat2it[iat1]].nw, ucell.atoms[ucell.iat2it[iat2]].nw );
			const int iw2_begin = ucell.iwt2iw[ iat2_range.first + index_begin.second ];
			for( int iw1_2D=iat1_range.first; iw1_2D<=iat1_range.second; ++iw1_2D )
			{
				const int iw1 = ucell.iwt2iw[ iw1_2D + index_begin.first ];
				memcpy(
					&DM_local( iw1, iw2_begin ), 
					&DM_2D( iw1_2D, iat2_range.first ), 
					( iat2_range.second - iat2_range.first + 1 )*sizeof(double) );
			}
			DM_local *= SPIN_multiple;
			
			if( DM_local.absmax() < threshold )	continue;
			
			if( matrix*const DM_ptr = static_cast<matrix*const>(MAP_EXIST( DMr[is], iat1, iat2, Abfs::Vector3_Order<int>{0,0,0} )) )
				*DM_ptr += DM_local;
			else
				DMr[is][iat1][iat2][{0,0,0}] = std::move(DM_local);
		}
	}
}




void Exx_Abfs::Parallel::Communicate::DM2::cal_DM_k( 
	const Abfs::Vector3_Order<int> &Born_von_Karman_period,
	const set<pair<size_t,size_t>> &H_atom_pairs_core,
	const double threshold )
{
	TITLE("Exx_Abfs::Parallel::Communicate::DM::cal_DM");
	
	vector<Abfs::Vector3_Order<int>> Born_von_Karman_boxes;
	for( int ix=0; ix<Born_von_Karman_period.x; ++ix )
		for( int iy=0; iy<Born_von_Karman_period.y; ++iy )
			for( int iz=0; iz<Born_von_Karman_period.z; ++iz )
				Born_von_Karman_boxes.push_back( Abfs::Vector3_Order<int>{ix,iy,iz} % Born_von_Karman_period );
			
	Exx_Abfs::DM dm_my;
	dm_my.flag_mix = false;
	dm_my.cal_DM( H_atom_pairs_core, Born_von_Karman_boxes );
	
	this->DMr.clear();
	this->DMr.resize(GlobalV::NSPIN);
	for( const auto &DMrA : dm_my.DMr )
	{
		const size_t iat1 = DMrA.first;
		for( const auto &DMrB : DMrA.second )
		{
			const size_t iat2 = DMrB.first;
			for( const auto &DMrC : DMrB.second )
			{
				const auto box2 = DMrC.first;
				for( size_t is=0; is!=DMrC.second.size(); ++is )
					if( DMrC.second[is].absmax() >= threshold )
						DMr[is][iat1][iat2][box2] = std::move(DMrC.second[is]);
			}
		}
	}
		
}
