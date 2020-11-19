#include "exx_abfs-parallel-distribute-order.h"

#include "src_pw/global.h"
#include "src_lcao/abfs.h"
#include <cmath>

vector<pair<size_t,size_t>> Exx_Abfs::Parallel::Distribute::Order::distribute(
	const double rmesh_times )
{
	const vector<Abfs::Vector3_Order<int>> Coulomb_potential_boxes = Abfs::get_Coulomb_potential_boxes(rmesh_times);
	auto neighbour = [&](const size_t iat1, const size_t iat2) -> bool
	{
		const int it1 = ucell.iat2it[iat1];
		const int it2 = ucell.iat2it[iat2];
		const Vector3<double> tau1 = ucell.atoms[it1].tau[ucell.iat2ia[iat1]];
		const Vector3<double> tau2 = ucell.atoms[it2].tau[ucell.iat2ia[iat2]];
		const double Rcut = std::min( ORB.Phi[it1].getRcut()*rmesh_times+ORB.Phi[it2].getRcut(), ORB.Phi[it1].getRcut()+ORB.Phi[it2].getRcut()*rmesh_times );
		for(const Vector3<int> box2 : Coulomb_potential_boxes)
		{
			const double R = (-tau1 + tau2 + box2 * ucell.latvec).norm();
			if(R*ucell.lat0 < Rcut)
				return true;
		}
		return false;
	};
	
	vector<pair<size_t,size_t>> pairs;
	for( size_t iat1=0; iat1<ucell.nat; ++iat1 )
		for( size_t iat2=iat1; iat2<ucell.nat; ++iat2 )
			if(neighbour(iat1,iat2))
				pairs.push_back( {iat1,iat2} );
			
	const size_t work = std::ceil(static_cast<double>(pairs.size())/NPROC);
	return vector<pair<size_t,size_t>>{ 
		std::min( pairs.begin()+ MY_RANK   *work, pairs.end() ),
		std::min( pairs.begin()+(MY_RANK+1)*work, pairs.end() )};
}