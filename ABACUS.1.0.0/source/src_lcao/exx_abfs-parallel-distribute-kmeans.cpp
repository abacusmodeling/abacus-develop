#include "exx_abfs-parallel-distribute-kmeans.h"
#include "src_pw/global.h"
#include <random>

vector<Exx_Abfs::Parallel::Distribute::Kmeans::Atom> 
Exx_Abfs::Parallel::Distribute::Kmeans::cluster( const int Nc )
{
	TITLE("Exx_Abfs::Parallel::Distribute::Kmeans::cluster");
	
	vector<Cluster> clusters(Nc+1);						// clusters[Nc] just for atoms init
	vector<Atom> atoms(ucell.nat);
	
ofstream ofs_mpi(exx_lcao.test_dir.process+"kmeans_"+TO_STRING(Nc)+"_"+TO_STRING(MY_RANK),ofstream::app);
	
	auto init = [&]() -> void
	{
ofs_mpi<<Nc<<endl;
		const double volumn = abs(ucell.a1.norm()*ucell.a2.norm()*ucell.a3.norm());
ofs_mpi<<volumn<<endl;
		const double rate = pow(Nc/volumn,1.0/3.0);
ofs_mpi<<rate<<endl;
ofs_mpi<<ucell.a1<<"\t"<<ucell.a2<<"\t"<<ucell.a3<<endl;
ofs_mpi<<ucell.a1.norm()<<"\t"<<ucell.a2.norm()<<"\t"<<ucell.a3.norm()<<endl;
		const int Nx = ceil(ucell.a1.norm()*rate);
		const int Ny = ceil(ucell.a2.norm()*rate);
		const int Nz = ceil(ucell.a3.norm()*rate);
ofs_mpi<<Nx<<"\t"<<Ny<<"\t"<<Nz<<endl;

		vector<bool> flag_is_center(Nx*Ny*Nz,true);
		for( int ic_big=Nc/2; ic_big<Nc/2+(Nx*Ny*Nz-Nc); ++ic_big )
			flag_is_center[ic_big] = false;
		
//		// flag_is_center has Nc trues and Nx*Ny*Nz-Nc falses randomly
//		vector<int> flag_is_center(Nc,true);
//		flag_is_center.resize(Nx*Ny*Nz,false);
//		std::shuffle( flag_is_center.begin(), flag_is_center.end(), std::default_random_engine(Nc*Nx*Ny*Nz) );	// Nc*Nx*Ny*Nz is just a seed, for all thread be the same

		
for( int ic_big=0; ic_big<flag_is_center.size(); ++ic_big )
	ofs_mpi<<flag_is_center[ic_big]<<"\t";
ofs_mpi<<endl;
		
		Vector3<double> taud_max = {0,0,0},
		                taud_min = {1,1,1};
		for( int iat=0; iat<ucell.nat; ++iat )
		{
			const Vector3<double> & taud = ucell.atoms[ ucell.iat2it[iat] ].taud[ ucell.iat2ia[iat] ];
			taud_max.x = max( taud.x, taud_max.x );
			taud_max.y = max( taud.y, taud_max.y );
			taud_max.z = max( taud.z, taud_max.z );
			taud_min.x = min( taud.x, taud_min.x );
			taud_min.y = min( taud.y, taud_min.y );
			taud_min.z = min( taud.z, taud_min.z );
		}
		const Vector3<double> taud_delta(
			(taud_max.x-taud_min.x)/Nx,
			(taud_max.y-taud_min.y)/Ny,
			(taud_max.z-taud_min.z)/Nz);
ofs_mpi<<taud_min<<"\t"<<taud_max<<"\t"<<taud_delta<<"\t"<<endl;			
		
		for( int ix=0,ic=0,ic_big=0; ix<Nx; ++ix )
			for( int iy=0; iy<Ny; ++iy )
				for( int iz=0; iz<Nz; ++iz,++ic_big )
					if(flag_is_center[ic_big])
					{
						clusters[ic].tau = ((ix+0.5)*taud_delta.x+taud_min.x)*ucell.a1 
						                 + ((iy+0.5)*taud_delta.y+taud_min.y)*ucell.a2 
						                 + ((iz+0.5)*taud_delta.z+taud_min.z)*ucell.a3;
						++ic;
					}
/*
		std::default_random_engine random_engine;
		random_engine.seed(time(0));
		std::uniform_real_distribution<double> uniform(0,1);
		for( int ic=0; ic<Nc; ++ic )
			clusters[ic].tau = uniform(random_engine)*ucell.a1
			                 + uniform(random_engine)*ucell.a2
			                 + uniform(random_engine)*ucell.a3;
*/							 
		for( int ic=0; ic<Nc; ++ic )
		{
			clusters[ic].tau_sum = {0,0,0};
			clusters[ic].size = 0;
		}
		clusters[Nc].tau_sum = {0,0,0};
		clusters[Nc].size = ucell.nat;
		
		for( int iat=0; iat<ucell.nat; ++iat )
		{
			atoms[iat].tau = ucell.atoms[ ucell.iat2it[iat] ].tau[ ucell.iat2ia[iat] ];
			atoms[iat].center = Nc;
			atoms[iat].distance = std::numeric_limits<double>::max();
		}
for( int ic=0; ic<Nc; ++ic )
	ofs_mpi<<clusters[ic].tau<<"\t";
ofs_mpi<<endl;
for( int iat=0; iat<ucell.nat; ++iat )
	ofs_mpi<<atoms[iat].center<<"\t";
ofs_mpi<<endl;		
	};
							
							
	auto update = [&]() -> bool
	{
//ofs_mpi<<__FILE__<<__LINE__<<endl;
		bool flag_atom_move = false;
		for( int iat=0; iat<ucell.nat; ++iat )
		{
			for( int ic=0; ic<Nc; ++ic )
			{
				const double distance_new = (atoms[iat].tau-clusters[ic].tau).norm2();
				if( distance_new < atoms[iat].distance )
				{
					const int ic_old = atoms[iat].center;
					clusters[ic_old].tau_sum -= atoms[iat].tau;
					--clusters[ic_old].size;
					
					atoms[iat].center = ic;
					atoms[iat].distance = distance_new;
					clusters[ic].tau_sum += atoms[iat].tau;
					++clusters[ic].size;
					
					flag_atom_move = true;
				}
			}
		}
		for( Cluster & cluster : clusters )
			cluster.tau = cluster.tau_sum / static_cast<double>(cluster.size);
		
for( int ic=0; ic<Nc; ++ic )
	ofs_mpi<<clusters[ic].tau<<"\t";
ofs_mpi<<endl;
for( int iat=0; iat<ucell.nat; ++iat )
	ofs_mpi<<atoms[iat].center<<"\t";
ofs_mpi<<endl;

		return flag_atom_move;	
	};

	init();
	while( update() );
	return atoms;
	
ofs_mpi.close();
}


vector<pair<size_t,size_t>> 
Exx_Abfs::Parallel::Distribute::Kmeans::distribute( const MPI_Comm & mpi_comm, const int multiple_core )
{
	TITLE("Exx_Abfs::Parallel::Distribute::Kmeans::distribute");
	
	assert(multiple_core>=1);
	int comm_size;	MPI_Comm_size( mpi_comm, &comm_size );
	int my_rank;	MPI_Comm_rank( mpi_comm, &my_rank );
	const int comm_size_nominal = comm_size * multiple_core;
	
//for( int N=2; N<=20; ++N )
//	cluster(N);
	
	{
		auto classify_atom = []( const int Ng, const vector<Atom> &atoms ) -> vector<vector<size_t>>
		{
			vector<vector<size_t>> clusters_atoms(Ng);
			for( size_t iat=0; iat<atoms.size(); ++iat )
				clusters_atoms[atoms[iat].center].push_back(iat);
			return clusters_atoms;
		};
		
		auto g_same_cluster = [](const vector<size_t> &cluster_atoms) -> vector<pair<size_t,size_t>>
		{
			vector<pair<size_t,size_t>> rank_work;
			for( size_t i1=0; i1<cluster_atoms.size(); ++i1 )
				for( size_t i2=i1; i2<cluster_atoms.size(); ++i2 )
					rank_work.push_back({cluster_atoms[i1],cluster_atoms[i2]});
			return rank_work;
		};
		auto g_different_cluster = [&]( const vector<size_t> &cluster_atoms1, const vector<size_t> &cluster_atoms2 ) -> vector<pair<size_t,size_t>>
		{
			vector<pair<size_t,size_t>> rank_work;
			for( const size_t iat1 : cluster_atoms1 )
				for( const size_t iat2 : cluster_atoms2 )
					rank_work.push_back({iat1,iat2});
			return rank_work;
		};
		
		const int Ng = static_cast<int>(sqrt(2*comm_size_nominal));
		if( Ng*(Ng+1)/2 == comm_size_nominal )
		{
											
			vector<pair<size_t,size_t>> rank_work;
			const vector<Atom> atoms = cluster(Ng);
			const vector<vector<size_t>> clusters_atoms = classify_atom(Ng,atoms);
			for( size_t ig1=0, rank_tmp=0; ig1<Ng; ++ig1 )
				for( size_t ig2=ig1; ig2<Ng; ++ig2, ++rank_tmp )
					if( rank_tmp%comm_size == my_rank )
					{
																						
						const vector<pair<size_t,size_t>> rank_work_tmp = 
							(ig1==ig2) ?
							g_same_cluster(clusters_atoms[ig1]) :
							g_different_cluster(clusters_atoms[ig1],clusters_atoms[ig2]);
						rank_work.insert( rank_work.end(), rank_work_tmp.begin(), rank_work_tmp.end() );
					}
			return rank_work;
		}
	}
	{
		int N1 = static_cast<int>(sqrt(comm_size_nominal));
		while(comm_size_nominal%N1)	--N1;
		const int N2=comm_size_nominal/N1;
		
		auto f = [&]( const vector<Atom> &atoms1, const vector<Atom> &atoms2 ) -> vector<pair<size_t,size_t>>
		{
			auto index = [&](const int ic1, const int ic2){ return (ic1*N2+ic2)%comm_size; };
			
			vector<vector<pair<size_t,size_t>>> rank_work(comm_size);
			for( size_t iatA=0; iatA<ucell.nat; ++iatA )
			{
				rank_work[index(atoms1[iatA].center, atoms2[iatA].center)].push_back({iatA,iatA});
				for( size_t iatB=iatA+1; iatB<ucell.nat; ++iatB )
					if( atoms1[iatA].center==atoms1[iatB].center && atoms2[iatA].center==atoms2[iatB].center )
						rank_work[index(atoms1[iatA].center, atoms2[iatB].center)].push_back({iatA,iatB});
			}
			for( size_t iatA=0; iatA<ucell.nat; ++iatA )
				for( size_t iatB=iatA+1; iatB<ucell.nat; ++iatB )
					if( atoms1[iatA].center!=atoms1[iatB].center || atoms2[iatA].center!=atoms2[iatB].center )
					{
						vector<pair<size_t,size_t>> & rank_work1 = rank_work[index(atoms1[iatA].center, atoms2[iatB].center)];
						vector<pair<size_t,size_t>> & rank_work2 = rank_work[index(atoms1[iatB].center, atoms2[iatA].center)];
						if( rank_work1.size() < rank_work2.size() )
							rank_work1.push_back({iatA,iatB});
						else
							rank_work2.push_back({iatB,iatA});
					}
			return rank_work[my_rank];
		};
		
		if(N1==N2)
		{
			const vector<Atom> atoms = cluster(N1);
			return f(atoms,atoms);
		}
		else
		{
			const vector<Atom> atoms1 = cluster(N1);
			const vector<Atom> atoms2 = cluster(N2);
			return f(atoms1,atoms2);
		}
	}
 
				 
}
