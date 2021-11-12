#include "pw_basis.h"

//
//distribute plane waves to different cores
//Known: G, GT, GGT, nx, ny, nz, poolnproc, poolrank, ggecut
//output: ig2isz[ig], istot2ixy[is], ixy2istot[ixy], is2ixy[is], ixy2ip[ixy], startnsz_per[ip], nstnz_per[ip], gg[ig], gcar[ig], gdirect[ig], nst, nstot
//
void PW_Basis::distribute_g()
{
    if(this->distribution_type == 1)
    {
        this->distribution_method1();
    }
    else
    {
        ModuleBase::WARNING_QUIT("divide", "No such division type.");
    }
    return;
}

//
// Distribute planewaves in reciprocal space to coreors.
// Firstly, divide the sphere in reciprocal space into sticks, which are vertical to x-y plane.
// Secondly, distribute these sticks to coreors.
//Known: G, GT, GGT, nx, ny, nz, poolnproc, poolrank, ggecut
//output: ig2isz[ig], istot2ixy[is], ixy2istot[nxy], is2ixy[is], ixy2ip[ixy], startnsz_per[ip], nst_per[ip], gg[ig], gcar[ig], gdirect[ig], nst
//
void PW_Basis::distribution_method1()
{
    // if use gamma point only, when convert real function f(r) to F(k) = FFT(f),
    // we have F(-k) = F(k)*, so that only half of planewaves are needed.
    if (this->gamma_only) this->nx = int(this->nx / 2) + 1;

    // initial the variables needed by all proc.
    const int nxy = this->nx * this->ny; // number of points in x-y plane.
    int tot_npw = 0;                     // total number of planewaves.
    int tot_nst = 0;                     // total number of sticks.
    int st_start = 0;                    // index of the first stick on current proc.
    int *st_i = NULL;                    // x or x + nx (if x < 0) of stick.
    int *st_j = NULL;                    // y or y + ny (if y < 0) of stick.
    int *st_bottom = NULL;               // minimum z of stick.
    int *st_length = NULL;               // number of planewaves in stick.
    int *is2ip = NULL;                   // ip of core that contains is^th stick, map is to ip.

    if (poolrank == 0)
    {
        // (1) Count the total number of planewaves (tot_npw) and sticks (tot_nst).                  
        
        // Actually we will scan [(2 * ibox[0] + 1) * (2 * ibox[1] + 1)] points on x-y plane,
        // but we define st_length2D with (nx * ny) points here, because we assume that the diameter
        // of the sphere is shorter than the sides of the cube.
        int *st_length2D = new int[nxy];                    // the number of planewaves that belong to the stick located on (x, y).
        int *st_bottom2D = new int[nxy];                    // the z-coordinate of the bottom of stick on (x, y).
        ModuleBase::GlobalFunc::ZEROS(st_length2D, nxy);
        ModuleBase::GlobalFunc::ZEROS(st_bottom2D, nxy);

        this->count_pw_st(tot_npw, tot_nst, st_length2D, st_bottom2D);
        this->nstot = tot_nst;

        // (2) Collect all planewaves and the x, y indexs, length, bottom of the sticks.
        double *gg_global = new double[tot_npw];                // the modulus of all planewaves.
        ModuleBase::Vector3<double>* gdirect_global = new ModuleBase::Vector3<double>[tot_npw]; // direct coordinates of planewaves.
        st_i = new int[tot_nst];                           // x or x + nx (if x < 0) of stick.
        st_j = new int[tot_nst];                           // y or y + ny (if y < 0) of stick.
        st_bottom = new int[tot_nst];                      // minimum z of stick.
        st_length = new int[tot_nst];                      // number of planewaves in stick.
        
        this->collect_pw_st(&tot_npw, &tot_nst, st_length2D, st_bottom2D, gg_global, gdirect_global, st_i, st_j, st_length, st_bottom);

        delete[] st_length2D;
        delete[] st_bottom2D;

#ifdef __MPI
        // Parallel line
        // (3) Distribute sticks to cores.
        int *npw_per = new int[this->poolnproc];  // number of planewaves on each core.
        int *nst_per = new int[this->poolnproc]; // number of sticks on each core.
        this->nstnz_per = new int[this->poolnproc]; // nz * nst(number of sticks) on each core.
        ModuleBase::GlobalFunc::ZEROS(npw_per, poolnproc);
        ModuleBase::GlobalFunc::ZEROS(nst_per, poolnproc);
        ModuleBase::GlobalFunc::ZEROS(this->nstnz_per, poolnproc);
        is2ip = new int[tot_nst];                 // ip of core that contains is^th stick, map is to ip.
        for (int i = 0; i < tot_nst; i++)
        {
            is2ip[i] = -1;                        // meaning this stick has not been distributed.
        }
        this->ixy2ip = new int[nxy];              // ip of core which contains stick on (x, y).
        for (int i = 0; i < nx * ny; i++)
        {
            this->ixy2ip[i] = -1;                 // meaning this stick has not been distributed or there is no stick on (x, y).
        }
        this->startnsz_per = new int[poolnproc];
        ModuleBase::GlobalFunc::ZEROS(startnsz_per, poolnproc);

        this->divide_sticks(tot_npw, tot_nst, st_i, st_j, st_length, npw_per, nst_per, is2ip);
        this->get_istot2ixy(tot_nst, st_i, st_j, is2ip);

        // (4) Divide planewaves to each core, construct gg2D and gdirect2D.
        double **gg2D = new double*[poolnproc];  // the i^th row contains the modulus of planewaves that belong to the i^th core.
        ModuleBase::Vector3<double> **gdirect2D = new ModuleBase::Vector3<double>*[poolnproc]; // the i^th row contains the direct coordinates of planewaves that belong to the i^th core.
        for (int ip = 0; ip < poolnproc; ip++)
        {
            gg2D[ip] = new double[npw_per[ip]];
            gdirect2D[ip] = new ModuleBase::Vector3<double>[npw_per[ip]];
        }

        this->divide_pw(tot_npw, gg_global, gdirect_global, gg2D, gdirect2D);

        delete[] gg_global;
        delete[] gdirect_global;

        // (5) Send gg, gdirect, npw_per, nst[poolrank], st_* to all cores.
        this->npw = npw_per[0];
        this->nst = nst_per[0];
        this->gg = gg2D[0];
        this->gdirect = gdirect2D[0];

        MPI_Bcast(&tot_npw, 1, MPI_INT, 0, POOL_WORLD);
        MPI_Bcast(&tot_nst, 1, MPI_INT, 0, POOL_WORLD);
        MPI_Bcast(&st_i, tot_nst, MPI_INT, 0, POOL_WORLD);
        MPI_Bcast(&st_j, tot_nst, MPI_INT, 0, POOL_WORLD);
        MPI_Bcast(&st_length, tot_nst, MPI_INT, 0, POOL_WORLD);
        MPI_Bcast(&st_bottom, tot_nst, MPI_INT, 0, POOL_WORLD);
        MPI_Bcast(&is2ip, tot_nst, MPI_INT, 0, POOL_WORLD);
        MPI_Bcast(&ixy2ip, nxy, MPI_INT, 0, POOL_WORLD);
        MPI_Bcast(&istot2ixy, tot_nst, MPI_INT, 0, POOL_WORLD);
        MPI_Bcast(&ixy2istot, nxy, MPI_INT, 0, POOL_WORLD);
        
        for (int ip = 1; ip < this->poolnproc; ip++)
        {
            // MPI_Send(&tot_npw, 1, MPI_INT, ip, 0, POOL_WORLD);
            // MPI_Send(&tot_nst, 1, MPI_INT, ip, 0, POOL_WORLD);
            MPI_Send(&npw_per[ip], 1, MPI_INT, ip, 0, POOL_WORLD);
            MPI_Send(&nst_per[ip], 1, MPI_INT, ip, 0, POOL_WORLD);
            MPI_Send(&gg2D[ip], npw_per[ip], MPI_DOUBLE, ip, 0, POOL_WORLD);
            MPI_Send(&gdirect2D[ip], npw_per[ip]*3, MPI_DOUBLE, ip, 0, POOL_WORLD); // I'm not sure about the send size and type here.
            // MPI_Send(&st_i, tot_nst, MPI_INT, ip, 0, POOL_WORLD);
            // MPI_Send(&st_j, tot_nst, MPI_INT, ip, 0, POOL_WORLD);
            // MPI_Send(&st_length, tot_nst, MPI_INT, ip, 0, POOL_WORLD);
            // MPI_Send(&st_bottom, tot_nst, MPI_INT, ip, 0, POOL_WORLD);
            // MPI_Send(&is2ip, tot_nst, MPI_INT, ip, 0, POOL_WORLD);
            // MPI_Send(&ixy2ip, nxy, MPI_INT, ip, 0, POOL_WORLD);
            // MPI_Send(&istot2ixy, tot_nst, MPI_INT, ip, 0, POOL_WORLD);
            // MPI_Send(&ixy2istot, nxy, MPI_INT, ip, 0, POOL_WORLD);
        }
        for (int ip = 0; ip < this->poolnproc; ip++)
        {
            delete[] gg2D[ip];
            delete[] gdirect2D[ip];
        }
        delete[] gg2D;
        delete[] gdirect2D;
        delete[] npw_per;
        delete[] nst_per;
#else
        // Serial line
        this->npw = tot_npw;
        this->nst = tot_nst;
        this->gg = gg_global;
        this->gdirect = gdirect_global;
        delete[] gg_global;
        delete[] gdirect_global;

        this->nstnz_per = new int[1]{0};
        this->startnsz_per = new int[1]{0};

        this->ixy2istot = new int[nxy];
        this->istot2ixy = new int[tot_nst];
        this->ixy2ip = new int[nxy];              // ip of core which contains stick on (x, y).
        for (int i = 0; i < nxy; i++)
        {
            ixy2istot[i] = -1;
            ixy2ip[i] = -1;
        }
        for (int is = 0; is < tot_nst; is++)
        {
            int index = st_i[is] + st_j[is] * nx;
            ixy2istot[index] = is;
            istot2ixy[is] = index;
            ixy2ip[index] = 0;
        }

        is2ip = new int[tot_nst];
        ModuleBase::GlobalFunc::ZEROS(is2ip, tot_nst);
#endif
    }
    else
    {
#ifdef __MPI
        MPI_Status ierror;
        // MPI_Recv(&tot_npw, 1, MPI_INT, 0, 0, POOL_WORLD, &ierror);
        // MPI_Recv(&tot_nst, 1, MPI_INT, 0, 0, POOL_WORLD, &ierror);
        MPI_Recv(&npw, 1, MPI_INT, 0, 0, POOL_WORLD, &ierror);  // number of planewaves in current proc.
        MPI_Recv(&nst, 1, MPI_INT, 0, 0, POOL_WORLD, &ierror);
        MPI_Recv(&gg, npw, MPI_DOUBLE, 0, 0, POOL_WORLD, &ierror);
        MPI_Recv(&gdirect, npw*3, MPI_DOUBLE, 0, 0, POOL_WORLD, &ierror); // I'm not sure about the send size and type here.
        // st_i = new int[tot_nst];
        // st_j = new int[tot_nst];
        // st_length = new int[tot_nst];
        // st_bottom = new int[tot_nst];
        // is2ip = new int [tot_nst];
        // MPI_Recv(&st_i, tot_nst, MPI_INT, 0, 0, POOL_WORLD, &ierror);
        // MPI_Recv(&st_j, tot_nst, MPI_INT, 0, 0, POOL_WORLD, &ierror);
        // MPI_Recv(&st_length, tot_nst, MPI_INT, 0, 0, POOL_WORLD, &ierror);
        // MPI_Recv(&st_bottom, tot_nst, MPI_INT, 0, 0, POOL_WORLD, &ierror);
        // MPI_Recv(&is2ip, tot_nst, MPI_INT, 0, 0, POOL_WORLD, &ierror);
        // MPI_Recv(&ixy2ip, nxy, MPI_INT, 0, 0, POOL_WORLD, &ierror);
        // MPI_Recv(&istot2ixy, tot_nst, MPI_INT, 0, 0, POOL_WORLD, &ierror);
        // MPI_Recv(&ixy2istot, nxy, MPI_INT, 0, 0, POOL_WORLD, &ierror);
#endif
    }
    
    this->nstnz = this->nst * this->nz;

    this->gcar = new ModuleBase::Vector3<double>[this->npw];
    for (int i = 0; i < this->npw; i++)
    {
        gcar[i] = gdirect[i] * this->G;
    }

    // (6) Construct ig2isz and is2ixy. 
    this->get_ig2isz_is2ixy(tot_nst, st_i, st_j, st_bottom, st_length, is2ip);

    if (st_i != NULL) delete[] st_i;
    if (st_j != NULL) delete[] st_j;
    if (st_bottom != NULL) delete[] st_bottom;
    if (st_length != NULL) delete[] st_length;
    if (is2ip != NULL) delete[] is2ip;

    return;
}

//
// (1) We count the total number of planewaves (tot_npw) and sticks (tot_nst) here.
// Meanwhile, we record the number of planewaves on (x, y) in st_length2D, and store the smallest z-coordinate of each stick in st_bottom2D,
// so that we can scan a much smaller area in step(2).
// known: nx, ny, nz, ggecut, GGT
// output: tot_npw, tot_nst, st_length2D, st_bottom2D
//
void PW_Basis::count_pw_st(
        int &tot_npw,     // total number of planewaves.
        int &tot_nst,     // total number of sticks.
        int* st_length2D, // the number of planewaves that belong to the stick located on (x, y).
        int* st_bottom2D  // the z-coordinate of the bottom of stick on (x, y).
)
{
    int ibox[3] = {0, 0, 0};                            // an auxiliary vector, determine the boundary of the scanning area.
    ibox[0] = int(this->nx / 2) + 1;                    // scan x from -ibox[0] to ibox[0].
    ibox[1] = int(this->ny / 2) + 1;                    // scan y from -ibox[1] to ibox[1].
    ibox[2] = int(this->nz / 2) + 1;                    // scan z from -ibox[2] to ibox[2].

    ModuleBase::Vector3<double> f;

    int ix_start = -ibox[0]; // determine the scaning area along x-direct, if gamma-only, only positive axis is used.
    int ix_end = ibos[0];
    if (this->gamma_only)
    {
        ix_start = 0;
        ix_end = this->nx;
    }
    for (int ix = ix_start; ix <= ix_end; ix++)
    {
        for (int iy = -ibox[1]; iy <= ibox[1]; iy++)
        {
            // we shift all sticks to the first quadrant in x-y plane here.
            // (ix, iy, iz) is the direct coordinates of planewaves.
            // x and y is the coordinates of shifted sticks in x-y plane.
            // for example, if nx = ny = 10, we will shift the stick on (-1, 2) to (9, 2),
            // so that its index in st_length and st_bottom is 9 * 10 + 2 = 20.
            int x = ix;
            int y = iy;
            if (x < 0) x += this->nx;
            if (y < 0) y += this->ny;
            int index = x * this->ny + y;

            int length = 0; // number of planewave in stick (x, y).
            for (int iz = -ibox[2]; iz <= ibox[2]; iz++)
            {
                f.x = ix;
                f.y = iy;
                f.z = iz;
                double modulus = f * (this->GGT * f);
                if (modulus <= this->ggecut)
                {
                    if (length == 0) st_bottom2D[index] = iz; // length == 0 means this point is the bottom of stick (x, y).
                    ++tot_npw;
                    ++length;
                }
            }
            if (length > 0)
            {
                st_length2D[index] = length;
                ++tot_nst;
            }
        }
    }
    return;
}

//        
// (2) Collect all planewaves and the x, y indexs, length, z_start of the sticks.
// Firstly, we collect the modulus of planewaves to gg_global, direct coordinates to gdirect_global,
// and construct temp_st_* simultaneously.
// Then, we will sort gg_global from lagest to smallest, and rearrange gdirect_global in the same order,
// that's why we define "temp"_gdirect.
// We define "temp"_st_* because we will rearrange them in the order of length decreasing, too.
// known: tot_npw, tot_nst, st_length2D, st_bottom2D
// output: gg_global, gdirect_global, st_i, st_j, st_length, st_bottom
//
void PW_Basis::collect_pw_st(
    const int tot_npw,                              // total number of planewaves.
    const int tot_nst,                              // total number of sticks.
    int* st_length2D,                               // the number of planewaves that belong to the stick located on (x, y), stored in 2d x-y plane.
    int* st_bottom2D,                               // the z-coordinate of the bottom of stick on (x, y), stored in 2d x-y plane.
    double* gg_global,                              // the modulus of all planewaves.
    ModuleBase::Vector3<double> *gdirect_global,    // direct coordinates of planewaves.
    int* st_i,                                      // x or x + nx (if x < 0) of stick.
    int* st_j,                                      // y or y + ny (if y < 0) of stick.
    int* st_length,                                 // number of planewaves in stick, stored in 1d array with tot_nst elements.
    int* st_bottom                                  // minimum z of stick, stored in 1d array with tot_nst elements.
)
{
    ModuleBase::Vector3<double> *temp_gdirect = new ModuleBase::Vector3<double>[tot_npw]; // direct coordinates of all planewaves, in the order of (x * ny * nz + y * nx + z).
    int *temp_st_i = new int[tot_nst];                      // x or x + nx (if x < 0) of stick.
    int *temp_st_j = new int[tot_nst];                      // y or y + ny (if y < 0) of stick.
    int *temp_st_bottom = new int[tot_nst];                 // minimum z of stick.
    ModuleBase::GlobalFunc::ZEROS(st_length, tot_nst);

    int ibox[3] = {0, 0, 0};                            // an auxiliary vector, determine the boundary of the scanning area.
    ibox[0] = int(this->nx / 2) + 1;                    // scan x from -ibox[0] to ibox[0].
    ibox[1] = int(this->ny / 2) + 1;                    // scan y from -ibox[1] to ibox[1].
    ibox[2] = int(this->nz / 2) + 1;                    // scan z from -ibox[2] to ibox[2].

    ModuleBase::Vector3<double> f;
    int ig = 0; // index of planewave.
    int is = 0; // index of stick.

    int ix_start = -ibox[0]; // determine the scaning area along x-direct, if gamma-only, only positive axis is used.
    int ix_end = ibos[0];
    if (this->gamma_only)
    {
        ix_start = 0;
        ix_end = this->nx;
    }
    for (int ix = ix_start; ix <= ix_end; ix++)
    {
        for (int iy = -ibox[1]; iy <= ibox[1]; iy++)
        {
            // we have shifted all sticks to the first quadrant in x-y plane before.
            // (ix, iy, iz) is the direct coordinates of planewaves.
            // x and y is the coordinates of shifted sticks in x-y plane.
            // for example, if nx = ny = 10, we will shift the stick on (-1, 2) to (9, 2),
            // so that its index in st_length and st_bottom is 9 * 10 + 2 = 20.
            int x = ix;
            int y = iy;
            if (x < 0) x += nx;
            if (y < 0) y += ny;
            int index = x * ny + y;
            if (st_length[index] > 0) // meaning there is a stick on (x, y) point.
            {
                for (int iz = st_bottom2D[index]; iz < st_bottom2D[index] + st_length2D[index]; iz++)
                {
                    f.x = ix;
                    f.y = iy;
                    f.z = iz;
                    double modulus = f * (GGT * f);
                    assert (modulus <= ggecut);
                    gg_global[ig] = modulus;
                    temp_gdirect[ig] = f;
                    ig++;                    
                }

                temp_st_i[is] = x;
                temp_st_j[is] = y;
                st_length[is] = st_length2D[index];
                temp_st_bottom[is] = st_bottom2D[index];
                is++;
            }
        }
    }
    assert(ig == tot_npw);
    assert(is == tot_nst);

    // Rearrange gg_global and gdirect in the order of modulus decreasing, and sort st_* from longest to shortest.
    // Firstly, we sort gg_global from lagest to smallest here, and get gdirect by rearranging gdirect_global in the same order.
    // Next, as we will distribute the longest sticks preferentially in Step(3), we rearrange st_* in the order of length decreasing.

    int *gg_sorted_index = new int[tot_npw];
    gg_sorted_index[0] = 0;
    ModuleBase::heapsort(tot_npw, gg_global, gg_sorted_index); // sort gg_global in the order of modulus decreasing.
    for (int i = 0; i < tot_npw; i++)
    {
        gdirect_global[i] = temp_gdirect[gg_sorted_index[i]]; // rearrange gdirect_global in the same order of gg_global.
    }

    int *st_sorted_index = new int[tot_nst];
    st_sorted_index[0] = 0;
    ModuleBase::heapsort(tot_nst, st_length, st_sorted_index); // sort st_* in the order of length decreasing.

    for (int i = 0; i < tot_nst; i++)
    {
        st_i[i] = temp_st_i[st_sorted_index[i]];
        st_j[i] = temp_st_j[st_sorted_index[i]];
        st_bottom[i] = temp_st_bottom[st_sorted_index[i]];
    }

    delete[] temp_gdirect;
    delete[] temp_st_i;
    delete[] temp_st_j;
    delete[] temp_st_bottom;
    delete[] gg_sorted_index;
    delete[] st_sorted_index;
    return;
}

//
// (3-1) Distribute sticks to cores.
// We have rearranged sticks in the order of length decreasing, so that we will distribute the longest stick preferentially here.
// For each stick, we find the core that contains the least planewaves firstly, and distribute the stick to it,
// then update npw_per, this->nstnz_per, is2ip this->ixy2ip and this->startnsz_per.
// known: tot_npw, tot_nst, st_i, st_j, st_length
// output: npw_per, nst_per, this->nstnz_per, is2ip, this->ixy2ip, this->startnsz_per
//
void PW_Basis::divide_sticks(
    const int tot_npw,  // total number of planewaves.
    const int tot_nst,  // total number of sticks.
    int* st_i,          // x or x + nx (if x < 0) of stick.
    int* st_j,          // y or y + ny (if y < 0) of stick.
    int* st_length,     // the stick on (x, y) consists of st_length[x*ny+y] planewaves.
    int* npw_per,       // number of planewaves on each core.
    int* nst_per,       // number of sticks on each core.
    int* is2ip          // ip of core containing is^th stick, map is to ip.
)
{
    int ipmin = 0; // The ip of core containing least number of planewaves.
    for (int is = 0; is < tot_nst; is++)
    {
        // find the ip of core containing the least planewaves.
        for (int ip = 0; ip < this->poolnproc; ip++)
        {
            const int ngrid = this->nx * this->ny * this->nrxx;     // number of real space planewaves on this core.
            const int non_zero_grid = nst_per[ip] * this->nz;       // number of reciprocal planewaves on this core.
            const int npwmin = npw_per[ipmin];
            const int npw_ip = npw_per[ip];
            const int nstmin = nst_per[ipmin];
            const int nst_ip = nst_per[ip];

            if (npw_ip == 0 || npw_ip < npwmin)
            {
                if (non_zero_grid + nz < ngrid) // assert reciprocal planewaves is less than real space planewaves.
                {
                    ipmin = ip;
                }
            }
            else if (npw_ip == npwmin && nst_ip < nstmin)
            {
                if (non_zero_grid + nz < ngrid)
                {
                    ipmin = ip;
                }
            }
        }
        nst_per[ipmin]++;
        this->nstnz_per[ipmin] += this->nz;
        npw_per[ipmin] += st_length[is];
        is2ip[is] = ipmin;
        this->ixy2ip[st_j[is] * this->nx + st_i[is]] = ipmin;
    }
    for (int ip = 1; ip < poolnproc; ip++)
    {
        this->startnsz_per[ip] = this->startnsz_per[ip - 1] + this->nstnz_per[ip - 1];
    }
    return;
}

//
// (3-2) Rearrange sticks in the order of the ip of core increasing, (st_start + st_move) is the new index of sticks.
// Then get istot2ixy (istot2ixy[is]: ix + iy * nx of is^th stick among all sticks) on the first core
// and ixy2istot (ixy2istot[ix + iy * nx]: is of stick on (ix, iy) among all sticks).
// known: tot_nst, st_i, st_j, is2ip, this->startnsz_per
// output: istot2ixy, ixy2istot
//
void PW_Basis::get_istot2ixy(
    const int tot_nst,  // total number of sticks.
    int* st_i,          // x or x + nx (if x < 0) of stick.
    int* st_j,          // y or y + ny (if y < 0) of stick.
    int* is2ip         // ip of core containing is^th stick, map is to ip.
)
{
    assert(this->poolrank == 0);
    this->istot2ixy = new int[tot_nst];
    ModuleBase::GlobalFunc::ZEROS(this->istot2ixy, poolnproc);
    this->ixy2istot = new int[this->nx * this->ny];
    for (int i = 0; i < this->nx * this->ny; i++)
    {
        this->ixy2istot[i] = -1;
    }

    int* st_move = new int[this->poolnproc]; // st_move[ip]: this is the st_move^th stick on ip^th core.
    ModuleBase::GlobalFunc::ZEROS(st_move, poolnproc);

    for (int is = 0; is < tot_nst; is++)
    {
        int ip = is2ip[is];
        this->istot2ixy[this->startnsz_per[ip] / this->nz + st_move[ip]] = st_j[is] * nx + st_i[is];
        this->ixy2istot[st_j[is] * nx + st_i[is]] = is;
        st_move[ip]++;
    }
    delete[] st_move;
    return;
}
//
// (4) Divide planewaves to each core, construct gg2D and gdirect2D.
// gg2D is a 2D array, the i^th row contains the planewaves that belong to the i^th core, and so is gdirect2D.
// Hence, we can send i^th row of gg2D and gdirect2D to i^th core instead of sending whole gg_global and gdirect_global.
// For each planewave, we map its x, y coordinates to the ip of core that contains it with this->ixy2ip.
// known: tot_npw, gg_global, gdirect_global, this->ixy2ip
// output: gg2D, gdirect2D
//
void PW_Basis::divide_pw(
    const int tot_npw,                          // total number of planewaves.
    double* gg_global,                          // the modulus of all planewaves.
    ModuleBase::Vector3<double>*gdirect_global, // direct coordinates of planewaves.
    double** gg2D,                               // the i^th row contains the modulus of planewaves that belong to the i^th core.
    ModuleBase::Vector3<double>**gdirect2D       // the i^th row contains the direct coordinates of planewaves that belong to the i^th core.
)
{
    int* pw_found = new int[this->poolnproc]; // how many planewaves have been found.
    ModuleBase::GlobalFunc::ZEROS(pw_found, this->poolnproc);
    for (int ig = 0; ig < tot_npw; ig++)
    {
        int x = gdirect_global[ig].x;
        int y = gdirect_global[ig].y;
        if (x < 0) x += nx;
        if (y < 0) y += ny;
        int ip = this->ixy2ip[y * nx + x];
        gg2D[ip][pw_found[ip]] = gg_global[ig];
        gdirect2D[ip][pw_found[ip]] = gdirect_global[ig];
        pw_found[ip]++;
    }
    delete[] pw_found;
    return;
}
//
// (6) Construct ig2isz and is2ixy.
// is2ixy contains the x-coordinate and y-coordinate of sticks on current core.
// ig2isz contains the z-coordinate of planewaves on current core.  
//
void PW_Basis::get_ig2isz_is2ixy(
    const int tot_nst,  // total number of sticks.
    int* st_i,          // x or x + nx (if x < 0) of stick.
    int* st_j,          // y or y + ny (if y < 0) of stick.
    int* st_bottom,     // minimum z of stick, stored in 1d array with tot_nst elements.
    int* st_length,     // the stick on (x, y) consists of st_length[x*ny+y] planewaves.
    int* is2ip          // ip of core containing is^th stick, map is to ip.
)
{
    this->ig2isz = new int[this->npw]; // map ig to the z coordinate of this planewave.
    ModuleBase::GlobalFunc::ZEROS(ig2isz, npw);
    this->is2ixy = new int[this->nst]; // map is (index of sticks) to ixy (ix + iy * nx).
    for (int i = 0; i < nst; i++) 
    {
        is2ixy[i] = -1;
    }

    int st_move = 0; // this is the st_move^th stick on current core.
    int pw_filled = 0; // how many current core's planewaves have been found.
    for (int is = 0; is < tot_nst; is++)
    {
        if (is2ip[is] == poolrank)
        {
            int zstart = st_bottom[is];
            for (int ig = 0; ig < st_length[is]; ig++)
            {
                int z = ig + zstart; // z-coordinate of this planewave.
                if (z < 0) z += nz;
                this->ig2isz[pw_filled] = st_move * nz + z;
                pw_filled++;
            }
            this->is2ixy[st_move] = st_j[is] * nx + st_i[is];
            st_move++;
        }
        if (st_move == this->nst && pw_filled == this->npw) break;
    }
    return;
}
