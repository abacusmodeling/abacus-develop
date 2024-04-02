#include "module_basis/module_nao/hydrogen_radials.h"
#include "module_base/global_variable.h"
#include "module_base/math_integral.h"
#include "module_base/atom_in.h" // for calculating slater screening constant
#include <map>
#include <iostream>
#include <algorithm>

HydrogenRadials& HydrogenRadials::operator=(const HydrogenRadials& rhs)
{
    RadialSet::operator=(rhs);
    return *this;
}

void HydrogenRadials::build(const int itype,
                            const double charge,
                            const bool with_slater_screening,
                            const int nmax,
                            const double rcut,
                            const double dr,
                            const double conv_thr,
                            const int rank,
                            const std::string symbol,
                            const std::string strategy,
                            std::ofstream* ptr_log)
{
    if(with_slater_screening) {printf("Build hydrogen_radials with Slater screening coefficients.\n");}
    cleanup();
    itype_ = itype;
    symbol_ = symbol;
    // rcut should be determined as soon as possible...
    //generate_hydrogen_radials(charge, nmax, 10.0, dr, rank, ptr_log);
    hydrogen(charge, with_slater_screening, nmax, dr, conv_thr, rank, strategy, ptr_log);
    set_rcut_max();
}

std::vector<double> HydrogenRadials::generate_hydrogen_radial_segment(const double charge,
                                                                      const bool with_slater_screening,
                                                                      const int n,
                                                                      const int l,
                                                                      const double rmin,
                                                                      const double rmax,
                                                                      const double dr,
                                                                      const int rank,
                                                                      std::ofstream* ptr_log)
{
    double a0 = 1.0; // Bohr radius
    int ngrid = static_cast<int>((rmax - rmin) / dr) + 1;
    std::vector<double> rgrid(ngrid);
    std::vector<double> rvalue(ngrid);

    // initialize value for rgrid
    for(int ir = 0; ir != ngrid; ++ir)
    {
        rgrid[ir] = rmin + ir * dr;
    }

    double screened_charge = charge;
    if(with_slater_screening)
    {
        double sigma = slater_screening(symbol_, n, l);
        screened_charge = charge - sigma;
    }

    double norm_factor = sqrt(
        4.0*std::pow(screened_charge, 3)*
        static_cast<double>(this->assoc_laguerre_.factorial(n - l - 1)) /
        std::pow(double(n), 4) / 
        static_cast<double>(this->assoc_laguerre_.factorial(n + l)) /
        std::pow(a0, 3)
    );

    for(int ir = 0; ir != ngrid; ++ir)
    {
        // Bohr radius is 1.0
        double rho = 2.0 * rgrid[ir] * screened_charge / n / a0;
        rvalue[ir] = norm_factor * std::pow(rho, l) * exp(-rho/2.0) * this->assoc_laguerre_.value(
            n, l, rho
        );
    }
    return rvalue;
}

double HydrogenRadials::radial_norm(const std::vector<double> rgrid,
                                    const std::vector<double> rvalue)
{
    std::vector<double> integrand(rvalue.size());
    for(int ir = 0; ir != rvalue.size(); ++ir)
    {
        integrand[ir] = rvalue[ir] * rvalue[ir] * rgrid[ir] * rgrid[ir];
    }
    double dr = rgrid[1] - rgrid[0];
    double norm = ModuleBase::Integral::simpson(rvalue.size(), integrand.data(), dr);
    norm = sqrt(norm);
    return norm;
}

double HydrogenRadials::generate_hydrogen_radial_toconv(const double charge,
                                                        const bool with_slater_screening,
                                                        const int n,
                                                        const int l,
                                                        const double conv_thr,
                                                        const int rank,
                                                        std::vector<double>& rgrid,
                                                        std::vector<double>& rvalue,
                                                        std::ofstream* ptr_log)
{
    double norm = 0.0;
    double rmax_ = 0.0; // in Bohr
    double rmin_ = 0.0; // always to be 0, in Bohr
    double delta_r = 0.5; // stepsize for radius cutoff searching, in Bohr

    // clear the input vectors
    rgrid.clear(); rgrid.shrink_to_fit();
    rvalue.clear(); rvalue.shrink_to_fit();
    double dr = 0.01; // radial function realspace grid stepsize, in Bohr
    if(delta_r < dr)
    {
        dr = delta_r;
    }
    printf("Searching for the cutoff radius for n = %d, l = %d, conv_thr = %6.4e\n", n, l, conv_thr);
    printf("%10s%12s%14s%18s", "Step Nr.", "Rmax (a.u.)", "Norm", "Delta Norm\n");
    int istep = 1;
    double delta_norm = 1.0;
    while((std::fabs(delta_norm) > conv_thr))
    {
        rmin_ = rmax_;
        rmax_ += delta_r;
        int ngrid = static_cast<int>((rmax_ - rmin_) / dr) + 1; // [rmin, rmax]
        std::vector<double> rgrid_segment(ngrid);
        for(int ir = 0; ir != ngrid; ++ir)
        {
            rgrid_segment[ir] = rmin_ + ir * dr;
        }
        std::vector<double> rvalue_segment = generate_hydrogen_radial_segment(
            charge, with_slater_screening, n, l, rmin_, rmax_, dr, rank, ptr_log);
        // before push back, pop back the last element
        if(rgrid.size() != 0)
        {
            rgrid.pop_back();
            rvalue.pop_back();
        }
        rgrid.insert(rgrid.end(), rgrid_segment.begin(), rgrid_segment.end());
        rvalue.insert(rvalue.end(), rvalue_segment.begin(), rvalue_segment.end());
        delta_norm = norm;
        norm = radial_norm(rgrid, rvalue);
        delta_norm = norm - delta_norm;
        if(istep == 0) printf("%10d%12.2f%14.10f%18.10e\n", istep, rmax_, norm, delta_norm);
        ++istep;
    }
    printf("...\n");
    printf("%10d%12.2f%14.10f%18.10e\n", istep, rmax_, norm, delta_norm);
    return rmax_;
}

std::vector<std::pair<int, int>> HydrogenRadials::unzip_strategy(const int nmax,
                                                                 const std::string strategy)
{
    if(strategy.substr(0, 6) != "energy")
    {
        // because for "energy", the nmax is used as the number of electrons
        assert(nmax < 8);
    }
    std::vector<std::pair<int, int>> nl_pairs;
    if(strategy.substr(0, 7) == "minimal")
    {
        if(strategy == "minimal-nodeless")
        {
            for(int n = 1; n <= nmax; n++)
            {
                std::pair<int, int> nl_pair = std::make_pair(n, n - 1);
                nl_pairs.push_back(nl_pair);
            }
        }
        else// if(strategy == "minimal-valence")
        {
            for(int l = 0; l < nmax; l++)
            {
                std::pair<int, int> nl_pair = std::make_pair(nmax, l);
                nl_pairs.push_back(nl_pair);
            }
        }
    }
    else if(strategy.substr(0, 6) == "energy")
    {
        // 1s, -(n+1)-> 2s, -(l+1)-> 2p, 3s, -(n+1)-> 3p, 4s, -(l+1)-> 3d, 4p, 5s, -(n+1)-> 4d, 5p, 6s
        int starting_n = 1;
        int starting_l = 0;
        int nelec = 0;
        int nl_switch = 0;
        while(nelec < nmax) // now nmax is the full number of electrons of certain atom
        {
            int n = starting_n;
            int l = starting_l;
            do
            {
                nl_pairs.push_back(std::make_pair(n, l));
                nelec += (2*l + 1)*2;
                n++;
                l--;
            } while ((l >= 0)&&(nelec < nmax));
            
            if(nl_switch%2 == 1)
            {
                starting_l++;
            }
            else
            {
                starting_n++;
            }
            nl_switch++;
        }
        if(strategy == "energy-valence")
        {
            std::vector<int> nmax_ls;
            std::vector<int> nmax_minus1_ls;
            int real_nmax = 0;
            for(auto nl_pair : nl_pairs)
            {
                if(nl_pair.first > real_nmax) real_nmax = nl_pair.first;
            }
            for(auto it = nl_pairs.begin(); it != nl_pairs.end();)
            {
                if(it->first == real_nmax) nmax_ls.push_back(it->second);
                else if(it->first == real_nmax - 1) nmax_minus1_ls.push_back(it->second);
                else
                {
                    it = nl_pairs.erase(it);
                    continue;
                }
                ++it;
            }
            for(auto it = nl_pairs.begin(); it != nl_pairs.end();)
            {
                if(it->first == real_nmax - 1)
                {
                    if(std::find(nmax_ls.begin(), nmax_ls.end(), it->second) != nmax_ls.end())
                    {
                        it = nl_pairs.erase(it);
                        continue;
                    }
                }
                ++it;
            }
        }
    }
    else
    {
        for(int n = 1; n <= nmax; n++)
        {
            for(int l = 0; l < n; l++)
            {
                std::pair<int, int> nl_pair = std::make_pair(n, l);
                nl_pairs.push_back(nl_pair);
            }
        }
    }
    return nl_pairs;
}

void HydrogenRadials::smooth(std::vector<double>& rgrid,
                             std::vector<double>& rvalue,
                             const double sigma)
{
    double prefactor = 1.0 / sqrt(2.0 * M_PI) / sigma;
    double rmax = rgrid.back();
    for(int ir = 0; ir != rgrid.size(); ++ir)
    {
        double delta_r = rgrid[ir] - rmax;
        double smooth = prefactor * exp(-delta_r * delta_r / 2.0 / sigma / sigma);
        rvalue[ir] *= (1 - smooth);
    }
}

std::map<std::pair<int, int>, std::pair<std::vector<double>, std::vector<double>>> 
HydrogenRadials::generate_orb(const double charge,
                              const bool with_slater_screening,
                              const int nmax,
                              const double dr,
                              const double conv_thr,
                              const int rank,
                              const std::string strategy,
                              std::ofstream* ptr_log)
{
    // create space for storing all generated orbitals
    // (n, l) to (rgrid, rvalue)
    std::map<std::pair<int, int>, std::pair<std::vector<double>, std::vector<double>>> radials; 
    std::map<std::pair<int, int>, double> rmaxs;

    std::vector<std::pair<int, int>> nl_pairs = unzip_strategy(nmax, strategy);
    double rmax = 0.0;
    for(auto nl_pair : nl_pairs)
    {
        int n = nl_pair.first;
        int l = nl_pair.second;
        std::vector<double> rgrid;
        std::vector<double> rvalue;
        double rmax_nl = generate_hydrogen_radial_toconv(charge,
                                                         with_slater_screening,
                                                         n,
                                                         l,
                                                         conv_thr,
                                                         rank,
                                                         rgrid,
                                                         rvalue,
                                                         ptr_log);
        radials[nl_pair] = std::make_pair(rgrid, rvalue);
        rmaxs[nl_pair] = rmax_nl;
        if(rmax < rmax_nl)
        {
            rmax = rmax_nl;
        }
    }
    // zero padding towards rmax
    for(auto& radial : radials)
    {
        int n = radial.first.first;
        int l = radial.first.second;
        std::pair<int, int> nl_pair = std::make_pair(n, l);

        std::vector<double>& rgrid = radial.second.first;
        std::vector<double>& rvalue = radial.second.second;

        if(rmaxs[nl_pair] < rmax)
        {
            int ngrid = static_cast<int>((rmax - rmaxs[nl_pair]) / dr); // (r, rmax]
            for(int ir = 1; ir <= ngrid; ++ir)
            {
                rgrid.push_back(rmaxs[nl_pair] + ir * dr);
                rvalue.push_back(0.0);
            }
        }
        // smooth the tail
        smooth(rgrid, rvalue, 0.1);
    }
    return radials;
}

std::map<std::pair<int, int>, std::pair<int, int>>
HydrogenRadials::mapping_nl_lzeta(const int nmax,
                                  const std::string strategy)
{
    std::map<std::pair<int, int>, std::pair<int, int>> nl_lzeta;
    std::vector<std::pair<int, int>> nl_pairs = unzip_strategy(nmax, strategy);
    // initialize nzetas by all zeros
    // get lmax
    int lmax = 0;
    for(auto nl_pair: nl_pairs)
    {
        lmax = lmax>nl_pair.second?lmax:nl_pair.second;
    }
    std::vector<int> nzetas(lmax + 1, 0);
    for(auto nl_pair : nl_pairs)
    {
        int n = nl_pair.first;
        int l = nl_pair.second;
        nl_lzeta[nl_pair] = std::make_pair(l, nzetas[l]);
        nzetas[l] += 1;
    }
    // import to standard RadialSet data structure
    lmax_ = lmax;
    nzeta_ = new int[lmax_ + 1];
    for(int l = 0; l <= lmax_; ++l)
    {
        nzeta_[l] = nzetas[l];
        nzeta_max_ = nzeta_max_>nzeta_[l]?nzeta_max_:nzeta_[l];
    }
    indexing();

    return nl_lzeta;
}

void HydrogenRadials::hydrogen(const double charge,
                               const bool with_slater_screening,
                               const int nmax,
                               const double dr,
                               const double conv_thr,
                               const int rank,
                               const std::string strategy,
                               std::ofstream* ptr_log)
{
    std::map<std::pair<int, int>, std::pair<std::vector<double>, std::vector<double>>> orbitals = 
        generate_orb(charge, with_slater_screening, nmax, dr, conv_thr, rank, strategy, ptr_log);
    std::map<std::pair<int, int>, std::pair<int, int>> nl_lzeta = mapping_nl_lzeta(nmax, strategy);

    nchi_ = orbitals.size();
    chi_ = new NumericalRadial[nchi_];

    int ichi = 0;
    for(auto orbital : orbitals)
    {
        int n = orbital.first.first;
        int l = orbital.first.second;
        std::pair<int, int> nl_pair = std::make_pair(n, l);
        std::vector<double>& rgrid = orbital.second.first;
        std::vector<double>& rvalue = orbital.second.second;
        int lzeta = nl_lzeta[nl_pair].second;
        chi_[index(l, lzeta)].build(l, true, rgrid.size(), rgrid.data(), rvalue.data(), 0, lzeta, symbol_, itype_, false);
        chi_[index(l, lzeta)].normalize();
        //++ichi;
    }
}

double HydrogenRadials::slater_screening(const std::string symbol,
                                         const int n,
                                         const int l)
{
    atom_in atom_db;
    std::vector<int> electron_config = atom_db.groundstate_electronconfiguration[symbol];
    int isubshell = 0;
    double sigma = 0.0;
    int _len = 0;
    for(int n_ = 1; n_ <= n; ++n_)
    {
        if(n_ == n) _len += l + 1;
        else _len += n_;
    }
    if(_len > electron_config.size())
    {
        printf("Error: electron configuration is not enough for %s\n", symbol.c_str());
        printf("n = %d, l = %d\n", n, l);
        exit(1);
    }
    // special case for 1s: for H and He, use 0.30 for 1s screening constant
    if(symbol == "H") return 0.0; // only one 1s electron, no screening by "other electrons"
    else if(symbol == "He") return 0.30; // only two 1s electrons, one screening the other
    else if(n == 1) return 0.30;
    for(int n_ = 1; n_ <= n; ++n_)
    {
        int lmax = (n_ == n) ? l : n_ - 1;
        for(int l_ = 0; l_ <= lmax; ++l_)
        {
            int nelec = electron_config[isubshell];
            if(n - n_ >= 2) sigma += nelec * 1.0;
            else if(n - n_ == 1)
            {
                double screening = (l > 1)? 1.00 : 0.85;
                sigma += nelec * screening;
            }
            else
            {
                if(l_ == l) sigma += (nelec - 1) * 0.35;
                else
                {
                    double screening = (l > 1)? 1.00 : 0.35;
                    sigma += nelec * screening;
                }
            }
            ++isubshell;
        }
    }
    return sigma;
}