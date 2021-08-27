#include "blas_connector.h"
#include "gtest/gtest.h"

#include <algorithm>
#include <array>
#include <complex>
#include <cstdlib>
TEST(blas_connector, sscal_)
{
	typedef float T;
	const int size = 8;
	const T scale = 2;
	const int incx = 1;
	std::array<T, size> result, answer;
	std::generate(result.begin(), result.end(), std::rand);
	for (int i = 0; i < size; i++)
		answer[i] = result[i] * scale;
	sscal_(&size, &scale, result.data(), &incx);
	for (int i = 0; i < size; i++)
		EXPECT_FLOAT_EQ(answer[i], result[i]);
}

TEST(blas_connector, dscal_)
{
	typedef double T;
	const int size = 8;
	const T scale = 2;
	const int incx = 1;
	std::array<T, size> result, answer;
	std::generate(result.begin(), result.end(), std::rand);
	for (int i = 0; i < size; i++)
		answer[i] = result[i] * scale;
	dscal_(&size, &scale, result.data(), &incx);
	for (int i = 0; i < size; i++)
		EXPECT_FLOAT_EQ(answer[i], result[i]);
}

TEST(blas_connector, cscal_)
{
	typedef std::complex<float> T;
	const int size = 8;
	const T scale = {2, 3};
	const int incx = 1;
	std::array<T, size> result, answer;
	std::generate(result.begin(), result.end(), []()
				  { return T{static_cast<float>(std::rand()),
							 static_cast<float>(std::rand())}; });
	for (int i = 0; i < size; i++)
		answer[i] = result[i] * scale;
	cscal_(&size, &scale, result.data(), &incx);
	for (int i = 0; i < size; i++)
	{
		EXPECT_FLOAT_EQ(answer[i].real(), result[i].real());
		EXPECT_FLOAT_EQ(answer[i].imag(), result[i].imag());
	}
}

TEST(blas_connector, zscal_)
{
	typedef std::complex<double> T;
	const int size = 8;
	const T scale = {2, 3};
	const int incx = 1;
	std::array<T, size> result, answer;
	std::generate(result.begin(), result.end(), []()
				  { return T{static_cast<double>(std::rand()),
							 static_cast<double>(std::rand())}; });
	for (int i = 0; i < size; i++)
		answer[i] = result[i] * scale;
	zscal_(&size, &scale, result.data(), &incx);
	for (int i = 0; i < size; i++)
	{
		EXPECT_FLOAT_EQ(answer[i].real(), result[i].real());
		EXPECT_FLOAT_EQ(answer[i].imag(), result[i].imag());
	}
}
