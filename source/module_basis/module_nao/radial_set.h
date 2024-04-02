#ifndef RADIAL_SET_H_
#define RADIAL_SET_H_

#include <map>
#include <memory>
#include <string>
#include <utility>

#include "module_base/spherical_bessel_transformer.h"
#include "module_basis/module_nao/numerical_radial.h"
#include "module_basis/module_ao/ORB_nonlocal.h"
#include "module_basis/module_ao/ORB_atomic.h"

/**
 * @brief An abstract class representing the set of all numerical radial
 * functions (of the same kind) from a single element.
 *
 * This class is supposed to be the base class for concrete classes like
 * AtomicRadials or BetaRadials, which represent the set of all numerical
 * atomic orbitals or Kleinman-Bylander beta functions of a single element
 * respectively.
 *
 * @see AtomicRadials BetaRadials
 */
class RadialSet
{
  public:
    RadialSet() = default;
    RadialSet(const RadialSet&);            ///< deep copy
    RadialSet& operator=(const RadialSet&); ///< deep copy
    virtual RadialSet* clone() const = 0;   ///< for polymorphic copy

    virtual ~RadialSet();

    /**
     * @name Builders.
     *
     * The derived classes should implement these functions so as to
     * set symbol_, itype_, nchi_ and chi_.
     *
     */
    ///@{
    /**
     * @brief Builds the object from a file.
     *
     * Currently only AtomicRadials objects are supposed to used this
     * interface.
     */
    virtual void build(const std::string&,             ///< file name
                       const int = 0,                  ///< the element index in calculation
                       std::ofstream* const = nullptr, ///< output file stream for logging
                       const int = 0                   ///< MPI rank
                       ) {}

    /**
     * @brief Builds from a Numerical_Nonlocal object.
     *
     * Currently nonlocal projectors are read in module_cell and passed
     * to Numerical_Nonlocal objects.
     */
    virtual void build(const Numerical_Nonlocal&,     ///< Numerical_Nonlocal object
                       const int = 0,                 ///< the element index in calculation
                       std::ofstream* const = nullptr ///< output file stream for logging
                       ) {}

    /**
     * @brief Builds from hydrogen-like radial functions.
     *
     * Currently only HydrogenRadials objects are supposed to used this
     * interface.
     */
    virtual void build(const int = 0,                         ///< the element index in calculation
                       const double = 1.0,                    ///< nuclear charge
                       const bool = false,                    ///< whether to include Slater screening
                       const int = 0,                         ///< maximal principal quantum number or electrons
                       const double = 10.0,                   ///< maximal radius
                       const double = 0.01,                   ///< radial grid step
                       const double = 1e-6,                   ///< convergence threshold for norm of pseudowavefunction
                       const int = 0,                         ///< MPI rank
                       const std::string = "",                ///< the name of the element
                       const std::string = "minimal-valence", ///< the strategy to generate whole set of radial functions
                       std::ofstream* const = nullptr         ///< output file stream for logging
                       ) {}

    /**
     * @brief Builds from pseudopotential file
     *
     * Currently only PswfcRadials objects are supposed to used this
     * interface.
     */
    virtual void build(const std::string&,             ///< file name
                       const int,                      ///< the element index in calculation
                       const double,                   ///< radial screening coefficient, sometimes to be 0.1
                       const double,                   ///< convergence threshold of norm of pseudowavefunction
                       std::ofstream* const = nullptr, ///< output file stream for logging
                       const int = 0                   ///< MPI rank
                       ) {}
    /**
     * @brief write any RadialSet object to a file in ABACUS numerical atomic orbital format.
     * 
     * This function will write any RadialSet object to file in ABACUS numerical atomic orbital format.
     * However its counterparts,  
     * `read_abacus_orb()` is in AtomicRadials,  
     * `read_beta_upf100()` and `read_beta_upf201()` are in BetaRadials,
     * `read_upf_pswfc()` is in PswfcRadials,  
     * `read_coeff()` is in SphbesRadials,  
     * due to read-in procedures always vary from one to another.
     */
    void write_abacus_orb(const std::string&,      ///< file name
                          const int = 0) const;    ///< MPI rank

    ///@}

    /**
     * @brief Overwrites the content of a Numerical_Orbital object with the current object.
     *
     * This function provides an interface to the corresponding object in the old module_ao.
     */
    virtual void to_numerical_orbital(Numerical_Orbital&,
                                      const int nk_legacy = 4005, // equivalent to lcao_ecut = 1600
                                      const double lcao_dk = 0.01
                                      ) const;

    /**
     * @name Getters
     */
    ///@{
    const std::string& symbol() const { return symbol_; }
    int itype() const { return itype_; }
    int lmax() const { return lmax_; }
    double rcut_max() const { return rcut_max_; }

    int nzeta(const int l) const { return (l >= 0 && l <= lmax_) ? nzeta_[l] : 0; }
    int nzeta_max() const { return nzeta_max_; }
    int nchi() const { return nchi_; }

    const NumericalRadial& chi(const int l, const int izeta);
    const NumericalRadial* cbegin() const { return chi_; }
    const NumericalRadial* cend() const { return chi_ + nchi_; }
    ///@}

    /**
     * @name Property setters for all NumericalRadial objects
     */
    ///@{
    /// Sets a spherical Bessel transformers for all NumericalRadial objects.
    void set_transformer(ModuleBase::SphericalBesselTransformer sbt, const int update = 0);

    /// Sets a common grid for all NumericalRadial objects.
    void set_grid(const bool for_r_space, const int ngrid, const double* grid, const char mode = 'i');

    /// Sets a common uniform grid for all NumericalRadial objects.
    void set_uniform_grid(const bool for_r_space,
                          const int ngrid,
                          const double cutoff,
                          const char mode = 'i',
                          const bool enable_fft = false);
    ///@}

  protected:
    std::string symbol_ = "";   ///< usually the chemical symbol
    int itype_ = 0;             ///< usually the index for element in calculation
    int lmax_ = -1;             ///< maximum angular momentum among all NumericalRadial objects
    double rcut_max_ = 0.0;     ///< maximum rcut (NOT rmax!) among all NumericalRadial objects

    int* nzeta_ = nullptr;      ///< number of NumericalRadial objects for each angular momentum
    int nzeta_max_ = 0;         ///< maximum number of NumericalRadial objects among each angular momentum
    int nchi_ = 0;              ///< total number of NumericalRadial objects

    NumericalRadial* chi_ = nullptr; ///< array of NumericalRadial objects

    /**
     * @brief A map from (l,izeta) to an index in chi_ array.
     *
     * An array of nzeta_max_ * (lmax_ + 1) elements which stores the index
     * of chi_ array for each (l,izeta) pair.
     */
    int* index_map_ = nullptr;

    /// Deallocates memory and reset all class members to default values.
    void cleanup();

    /// Gets the index in chi_ array from (l,izeta).
    int index(const int l, const int izeta) const;

    /// Builds index_map_ from nzeta_, nzeta_max_ and lmax_.
    void indexing();

    /// Sets rcut_max_ to be the maximum rcut of all NumericalRadial objects.
    void set_rcut_max();
};

#endif
