#include "exx_abfs-parallel-distribute-htime.h"

#include "src_pw/global.h"
#include "src_lcao/abfs.h"
#include <algorithm>

vector<pair<size_t,size_t>> Exx_Abfs::Parallel::Distribute::Htime::distribute( 
	const Abfs::Vector3_Order<int> & Born_von_Karman_period, const double rmesh_times )
{
	TITLE("Exx_Abfs::Parallel::distribute");
	const vector<size_t> Nadj = cal_Nadj(Born_von_Karman_period);

//for( const size_t & i : Nadj )
//	cout<<i<<endl;

	const vector<pair<size_t,pair<size_t,size_t>>> pair_costs = cal_pair_costs(Nadj,rmesh_times);
	
//for( const auto & i : pair_costs )
//	cout<<i.first<<"\t"<<i.second.first<<"\t"<<i.second.second<<endl;
	
	const vector<vector<pair<size_t,size_t>>> rank_work = cal_rank_work(pair_costs);
	
//for( size_t irank=0; irank!=rank_work.size(); ++irank )
//{
//	cout<<irank<<endl;
//	for( const auto & work : rank_work[irank] )
//		cout<<"\t"<<work.first<<"\t"<<work.second<<endl;
//}	
	
	return rank_work[MY_RANK];
}

// Nadj[iat]
vector<size_t> Exx_Abfs::Parallel::Distribute::Htime::cal_Nadj( const Abfs::Vector3_Order<int> & Born_von_Karman_period )
{
	TITLE("Exx_Abfs::Parallel::cal_Nadj");
	vector<size_t> Nadj(ucell.nat);
	for( size_t iat=0; iat!=ucell.nat; ++iat )
	{
		const map<size_t,vector<Abfs::Vector3_Order<int>>> 
			adjs = Abfs::get_adjs(iat);
		for( const auto & adj_i : adjs )
		{
			set<Abfs::Vector3_Order<int>> boxp;
			for( const auto & box : adj_i.second )
				boxp.insert( box % Born_von_Karman_period );
			Nadj[iat] += boxp.size();
		}
	}
	return Nadj;
}

// { Ni*Nj, {i,j} }
vector<pair<size_t,pair<size_t,size_t>>> Exx_Abfs::Parallel::Distribute::Htime::cal_pair_costs( 
	const vector<size_t> &Nadj, const double rmesh_times )
{
	TITLE("Exx_Abfs::Parallel::cal_pair_costs");
	
	vector<Vector3<int>> boxes;
	for(const int ix:{-1,0,1})
		for(const int iy:{-1,0,1})
			for(const int iz:{-1,0,1})
				boxes.push_back({ix,iy,iz});
			
	auto neighbour = [&](const size_t iat1, const size_t iat2) -> bool
	{
		const int it1 = ucell.iat2it[iat1];
		const int it2 = ucell.iat2it[iat2];
		const Vector3<double> tau1 = ucell.atoms[it1].tau[ucell.iat2ia[iat1]];
		const Vector3<double> tau2 = ucell.atoms[it2].tau[ucell.iat2ia[iat2]];
		double R_min = std::numeric_limits<double>::max();
		for(const Vector3<int> box2 : boxes)
		{
			const double R = (-tau1 + tau2 + box2 * ucell.latvec).norm();
			if(R<R_min)
				R_min = R;
		}
		return R_min < ORB.Phi[it1].getRcut()*rmesh_times+ORB.Phi[it2].getRcut() &&
			   R_min < ORB.Phi[it1].getRcut()+ORB.Phi[it2].getRcut()*rmesh_times;
	};
			
	vector<pair<size_t,pair<size_t,size_t>>> pair_costs;
	for( size_t iat1=0; iat1<ucell.nat; ++iat1 )
		for( size_t iat2=iat1; iat2<ucell.nat; ++iat2 )
			if(neighbour(iat1,iat2))
				pair_costs.push_back( {Nadj[iat1]*Nadj[iat2], {iat1,iat2} } );
		
	auto comp = []( 
		const pair<size_t,pair<size_t,size_t>> & pair_cost1, 
		const pair<size_t,pair<size_t,size_t>> & pair_cost2 ) -> bool
	{	return pair_cost1.first > pair_cost2.first; };
	std::sort( pair_costs.begin(), pair_costs.end(), comp );
	return pair_costs;
}

vector<vector<pair<size_t,size_t>>> Exx_Abfs::Parallel::Distribute::Htime::cal_rank_work( const vector<pair<size_t,pair<size_t,size_t>>> & pair_costs )
{
	TITLE("Exx_Abfs::Parallel::cal_rank_work");
	vector<pair<size_t,size_t>> rank_cost(NPROC);				// rank_cost[i] = { irank, cost }
	for( size_t irank=0; irank!=NPROC; ++irank )
		rank_cost[irank] = { irank, 0 };
	
	auto comp = [](
		const pair<size_t,size_t> & rank_cost1, 
		const pair<size_t,size_t> & rank_cost2 )
	{	return rank_cost1.second > rank_cost2.second;	};
	
	vector<vector<pair<size_t,size_t>>> rank_work(NPROC);		// rank_work[irank] = { {iat1,iat2}, {iat1,iat2}, ... }
	for( const auto pair_cost : pair_costs )
	{
		pop_heap( rank_cost.begin(), rank_cost.end(), comp );
//cout<<rank_cost.back().first<<"\t"<<rank_cost.back().second<<"\t"<<pair_cost.first<<"\t"<<pair_cost.second.first<<"\t"<<pair_cost.second.second<<endl;
		rank_cost.back().second += pair_cost.first;
		rank_work[ rank_cost.back().first ].push_back( pair_cost.second );
		push_heap( rank_cost.begin(), rank_cost.end(), comp );
	}
	return rank_work;
}