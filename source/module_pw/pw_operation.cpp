#include "pw_basis_k.h"
#include "pw_basis.h"
#include <assert.h>

ModuleBase::Vector3<double> PW_Basis_K:: get_GPlusK_cartesian(const int ik, const int ig) const {
    assert(ig>=0 && ig<this->npw && ik>=0 && ik<this->nks);
    ModuleBase::Vector3<double> g_temp_ = this->kvec_c[ik] + this->gcar[ig];
    return g_temp_;
}


double PW_Basis_K::get_GPlusK_cartesian_projection(const int ik, const int ig, const int axis) const
{
    assert(ig >= 0 && ig < this->npw && ik >= 0 && ik<this->nks && axis >= 0 && axis <= 2);
    ModuleBase::Vector3<double> g_temp_ = this->kvec_c[ik] + this->gcar[ig];
    if (axis == 0)
    {
        return g_temp_.x;
    }
    else if (axis == 1)
    {
        return g_temp_.y;
    }
    else if (axis == 2)
    {
        return g_temp_.z;
    }
    return 0.0;
}

double PW_Basis_K::get_SquareGPlusK_cartesian(const int ik, const int ig) const 
{
    assert(ig>=0 && ig<this->npw && ik>=0 && ik<this->nks);
    ModuleBase::Vector3<double> g_temp_ = this->kvec_c[ik] + this->gcar[ig];
    return (g_temp_ * g_temp_);
}