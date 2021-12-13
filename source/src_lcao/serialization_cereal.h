#ifndef SERIALIZATION_CEREAL_H
#define SERIALIZATION_CEREAL_H

#include <cereal/archives/binary.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/functional.hpp>

#include "../module_base/vector3.h"
#include "../src_ri/abfs-vector3_order.h"
#include "../module_base/matrix.h"
#include "src_io/read_txt_input_value.h"
#include "src_io/read_txt_input_item.h"



template<class Archive, typename T> void serialize( Archive & ar, Abfs::Vector3_Order<T> & v ){ ar(v.x); ar(v.y); ar(v.z); }

namespace ModuleBase
{

template<class Archive, typename T> void serialize( Archive & ar, Vector3<T> & v ){ ar(v.x); ar(v.y); ar(v.z); }

template<class Archive> void save( Archive & ar, const matrix & m )
{
	ar(m.nr);	ar(m.nc);
	ar(cereal::binary_data(m.c, m.nr*m.nc*sizeof(double)));
}
template<class Archive> void load( Archive & ar, matrix & m )
{
	int nr, nc;
	ar(nr);	ar(nc);
	m.create(nr,nc);
	ar(cereal::binary_data(m.c, m.nr*m.nc*sizeof(double)));
}
}

namespace Read_Txt_Input
{
template<class Archive> void Input_Value::serialize( Archive & ar )
{
	ar( this->b );
	ar( this->i );
	ar( this->d );
	ar( this->s );
}
template<class Archive> void Input_Item::serialize( Archive & ar )
{
	ar( this->annotation );
	ar( this->value );
	ar( this->label );
	ar( this->value_read_size );
}
}

#endif
