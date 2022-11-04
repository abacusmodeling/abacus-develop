#ifndef VDW_H
#define VDW_H

#include <vector>
#include "module_cell/unitcell_pseudo.h"
#include "module_vdw/vdw_parameters.h"
#include "module_vdw/vdwd2_parameters.h"
#include "module_vdw/vdwd3_parameters.h"

namespace vdw
{

template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args &&... args) {
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

class Vdw
{
  public:
    Vdw(const UnitCell_pseudo &unit_in) : ucell_(unit_in) {};

    virtual ~Vdw(){};

    inline double get_energy(bool cal=true) {
        if (cal) { cal_energy(); }
        return energy_;
    }
    inline const std::vector<ModuleBase::Vector3<double>> &get_force(bool cal=true) {
        if (cal) { cal_force(); }
        return force_;
    }
    inline const ModuleBase::Matrix3 &get_stress(bool cal=true) {
        if (cal) { cal_stress(); }
        return stress_;
    }

  protected:
    const UnitCell_pseudo &ucell_;

    double energy_ = 0;
    std::vector<ModuleBase::Vector3<double>> force_;
    ModuleBase::Matrix3 stress_;

    virtual void cal_energy() { throw std::runtime_error("No cal_energy method in base Vdw class"); }
    virtual void cal_force() { throw std::runtime_error("No cal_energy method in base Vdw class"); }
    virtual void cal_stress() { throw std::runtime_error("No cal_energy method in base Vdw class"); }
};

std::unique_ptr<Vdw> make_vdw(const UnitCell_pseudo &ucell, const Input &input);

} // namespace vdw

#endif // VDW_H
