#include "gtest/gtest.h"
#include "module_base/blas_connector.h"

#include <vector>
#include <cstdlib>
#include <algorithm>

TEST(blas_connector, dscal_)
{
	const int size = 8, scale = 2;
	std::vector<double> a(size), result(size), answer(size);
	std::generate(a.begin(), a.end(), std::rand);
	const int incx = 1;
	dscal_(&size, a.data(), result.data(), &incx);
	for (auto i : a)
		answer.push_back(i * scale);
	for (int i = 0; i < size; i++)
		EXPECT_FLOAT_EQ(answer[i], result[i]);
}