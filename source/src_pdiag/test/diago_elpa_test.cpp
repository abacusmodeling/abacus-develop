
#include "diago_elpa_utils.h"
#include "mpi.h"

#include "gtest/gtest.h"
#include <type_traits>

#define DETAILINFO       true
#define PASSTHRESHOLD    1e-4
#define PRINT_NEW_MATRIX false // print out the H matrix after doing block-cycle distribution

/************************************************
 *  unit test of ELPA
 ***********************************************/

/**
 * Test the diagonalization by ELPA with parallel
 * H/S matrix is complex double or double.
 */

template <class T> class ElpaPrepare
{
  public:
    ElpaPrepare(int n, int nblk, int sparsity) : n(n), nblk(nblk), sparsity(sparsity)
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &mypnum);
    }

    ElpaPrepare()
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &mypnum);
    }

    ~ElpaPrepare()
    {
        delete[] hmatrix;
        delete[] smatrix;
        delete[] e_lapack;
        delete[] e_elpa;
    }

    int n = 10; // dimension of H/S matrix
    int nblk = 1; // dimension of sub-matrix
    int sparsity = 100; // sparsity of H matrix;
    T *hmatrix = 0;
    T *smatrix = 0;
    double *e_lapack = 0; // eigenvalues solvered by LAPACK
    double *e_elpa = 0;
    double lapack_time = 0.0, elpa_time = 0.0;
    int mypnum = 0;

    void random_HS()
    {
        LCAO_DIAGO_TEST::random_hs<T>(hmatrix, smatrix, n, sparsity);
    }

    inline void readmatrix(std::ifstream &inf, T *&matrix, int &npw)
    {
        inf >> npw;
        matrix = new T[npw * npw];
        for (int i = 0; i < npw; i++)
        {
            for (int j = i; j < npw; j++)
            {
                inf >> matrix[i * npw + j];
                if (i != j)
                {
                    matrix[j * npw + i] = matrix[i * npw + j];
                }
            }
        }
    }

    void read_HS()
    {
        if (mypnum == 0)
        {
            int hdim = 1;
            int sdim = 1;

            std::ifstream ifs;
            if (std::is_same<T, double>::value)
            {
                ifs.open("H-GammaOnly.dat");
            }
            else if ((std::is_same<T, std::complex<double>>::value))
                ifs.open("H-KPonts.dat");
            readmatrix(ifs, hmatrix, hdim);
            ifs.close();

            if (std::is_same<T, double>::value)
            {
                ifs.open("S-GammaOnly.dat");
            }
            else if ((std::is_same<T, std::complex<double>>::value))
                ifs.open("S-KPoints.dat");
            readmatrix(ifs, smatrix, sdim);
            ifs.close();

            if (hdim != sdim)
            {
                printf("Error: dimensions of H and S are not equal, %d, %d", hdim, sdim);
                exit(1);
            }
            n = hdim;
            nblk = 1;
        }
        MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&nblk, 1, MPI_INT, 0, MPI_COMM_WORLD);
        e_lapack = new double[n];
        e_elpa = new double[n];
        for (int i = 0; i < n; i++)
        {
            e_lapack[i] = 0.0;
            e_elpa[i] = 1.0;
        }
    }

    void new_matrix()
    {
        hmatrix = new T[n * n];
        smatrix = new T[n * n];
        e_lapack = new double[n];
        e_elpa = new double[n];
        for (int i = 0; i < n; i++)
        {
            e_lapack[i] = 0.0;
            e_elpa[i] = 1.0;
        }
    }

    void print_hs()
    {
        if (mypnum == 0)
        {
            std::ofstream fp("hmatrix.dat");
            LCAO_DIAGO_TEST::print_matrix(fp, hmatrix, n, n, true);
            fp.close();
            fp.open("smatrix.dat");
            LCAO_DIAGO_TEST::print_matrix(fp, smatrix, n, n, true);
            fp.close();
        }
    }

    void lapack_ev()
    {
        double starttime, endtime;
        starttime = MPI_Wtime();
        LCAO_DIAGO_TEST::lapack_diago(hmatrix, smatrix, e_lapack, n);
        endtime = MPI_Wtime();
        lapack_time = endtime - starttime;
    }

    void elpa_ev()
    {
        double starttime, endtime;
        starttime = MPI_Wtime();
        LCAO_DIAGO_TEST::elpa_diago<T>(hmatrix, smatrix, e_elpa, n, nblk, PRINT_NEW_MATRIX, DETAILINFO);
        endtime = MPI_Wtime();
        elpa_time = endtime - starttime;
    }

    void run_diago()
    {
        if (this->mypnum == 0)
        {
            this->lapack_ev();
        }
        this->elpa_ev();
        if (this->mypnum == 0)
        {
            double maxerror = 0.0;
            int iindex = 0;
            for (int i = 0; i < this->n; i++)
            {
                EXPECT_NEAR(this->e_lapack[i], this->e_elpa[i], PASSTHRESHOLD);
                if (abs(this->e_lapack[i] - this->e_elpa[i]) > maxerror)
                {
                    maxerror = abs(this->e_lapack[i] - this->e_elpa[i]);
                    iindex = i;
                }
            }
            if (DETAILINFO)
            {
                std::cout << "n=" << this->n << ", nblk=" << this->nblk;
                if (sparsity != 100)
                {
                    std::cout << ", H sparsity=" << this->sparsity;
                }
                std::cout << std::endl;
                std::cout << "ELPA time(+data distribution): " << this->elpa_time
                          << "s, LAPACK time(1 core):" << this->lapack_time << "s," << std::endl;
                std::cout << "Maximum difference between ELPA and LAPACK is " << maxerror << " (" << iindex
                          << "-th eigenvalue)" << std::endl;
            }
        }
    }
};

class ElpaDoubleTest : public ::testing::TestWithParam<ElpaPrepare<double>>
{
};
class ElpaComplexDoubleTest : public ::testing::TestWithParam<ElpaPrepare<std::complex<double>>>
{
};

TEST(VerifyElpaDiaoDouble, ReadHS)
{
    ElpaPrepare<double> edp;
    edp.read_HS();
    // edp.print_hs();
    edp.run_diago();
}

TEST(VerifyElpaDiaoComplexDouble, ReadHS)
{
    ElpaPrepare<std::complex<double>> edp;
    edp.read_HS();
    // edp.print_hs();
    edp.run_diago();
}

TEST_P(ElpaDoubleTest, RandomHS)
{
    ElpaPrepare<double> edp = GetParam();
    edp.new_matrix();
    edp.random_HS();
    edp.print_hs();
    edp.run_diago();
}

INSTANTIATE_TEST_SUITE_P(VerifyElpaB2DDiag,
                         ElpaDoubleTest,
                         ::testing::Values(ElpaPrepare<double>(5, 1, 0),
                                           // ElpaPrepare<double>(100,1,7),
                                           // ElpaPrepare<double>(500,1,0),
                                           ElpaPrepare<double>(500, 1, 7),
                                           ElpaPrepare<double>(500, 32, 7)
                                           // ElpaPrepare<double>(1000,64,10)
                                           // ElpaDoublePrepare<double>(2000,128,7)
                                           ));

TEST_P(ElpaComplexDoubleTest, RandomHS)
{
    ElpaPrepare<std::complex<double>> edp = GetParam();
    edp.new_matrix();
    edp.random_HS();
    edp.print_hs();
    edp.run_diago();
}

INSTANTIATE_TEST_SUITE_P(VerifyElpaB2DDiag,
                         ElpaComplexDoubleTest,
                         ::testing::Values(ElpaPrepare<std::complex<double>>(5, 1, 0),
                                           // ElpaPrepare<std::complex<double>>(100,1,7),
                                           // ElpaPrepare<std::complex<double>>(500,1,0),
                                           ElpaPrepare<std::complex<double>>(500, 1, 7),
                                           ElpaPrepare<std::complex<double>>(500, 32, 7)
                                           // ElpaPrepare<std::complex<double>>(800,64,9)
                                           // ElpaDoublePrepare<double>(2000,128,7)
                                           ));

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    testing::InitGoogleTest(&argc, argv);

    int mypnum;
    MPI_Comm_rank(MPI_COMM_WORLD, &mypnum);
    ::testing::TestEventListeners &listeners = ::testing::UnitTest::GetInstance()->listeners();
    if (mypnum != 0)
    {
        delete listeners.Release(listeners.default_result_printer());
    }

    int result = RUN_ALL_TESTS();
    if (mypnum == 0 && result != 0)
    {
        std::cout << "ERROR:some tests are not passed" << std::endl;
        return result;
    }
    else
    {
        MPI_Finalize();
        return 0;
    }
}
