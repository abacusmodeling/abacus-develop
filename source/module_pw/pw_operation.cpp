ModuleBase::Vector3<double> get_GPlusK_cartesian(const int ik, const int ig) const {
    assert(ig>=0 && ig<this->ngmc && ik>=0 && ik<Klist->nks);
    ModuleBase::Vector3<double> g_temp_ = Klist->kvec_c[ik] + this->gcar[ig];
    return g_temp_;
};
double get_GPlusK_cartesian_projection(const int ik, const int ig, const int axis) const
{
    assert(ig >= 0 && ig < this->ngmc && ik >= 0 && ik < Klist->nks && axis >= 0 && axis <= 2);
    ModuleBase::Vector3<double> g_temp_ = Klist->kvec_c[ik] + this->gcar[ig];
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
double get_SquareGPlusK_cartesian(const int ik, const int ig) const 
{
    assert(ig >= 0 && ig < this->ngmc && ik >= 0 && ik < Klist->nks);
    ModuleBase::Vector3<double> g_temp_ = Klist->kvec_c[ik] + this->gcar[ig];
    return (g_temp_ * g_temp_);
};
ModuleBase::Vector3<double> get_G_cartesian(const int ig) const 
{
    assert(ig>=0 && ig<this->ngmc);
    return this->gcar[ig];
};
double get_G_cartesian_projection(const int ig, const int axis) const 
{
    assert(ig>=0 && ig<this->ngmc && axis>=0 && axis<=2);
    if(axis == 0) 
    {
        return this->gcar[ig].x;
    }
    else if(axis == 1)
    {
        return this->gcar[ig].y;
    }
    else if(axis == 2)
    {
        return this->gcar[ig].z;
    }
    return 0.0;
}
double get_NormG_cartesian(const int ig) const
{
    assert(ig >= 0 && ig < this->ngmc);
    return (this->gcar[ig].x * this->gcar[ig].x + this->gcar[ig].y * this->gcar[ig].y + this->gcar[ig].z * this->gcar[ig].z);
}