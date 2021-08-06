#include "module_base/matrix3.h"
#include "gtest/gtest.h"
#include <random>
#include <vector>
class matrix3_test : public testing::Test
{
protected:
	Matrix3 matrix_a, matrix_a1, matrix_b;
	Matrix3 get_random_matrix3()
	{
		vector<double> v(9);
		for (auto &i : v)
		{
			i = std::rand();
		}
		auto matrix_a = Matrix3(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8]);
		return matrix_a;
	}
	void SetUp()
	{
		matrix_a = get_random_matrix3();
		// matrix_a.print();
		matrix_a1 = matrix_a;
		// matrix_a1.print();
		matrix_b = get_random_matrix3();
		// matrix_b.print();
		return;
	};
};

TEST_F(matrix3_test, op_equal)
{
	ASSERT_TRUE(matrix_a == matrix_a1);
	ASSERT_FALSE(matrix_a == matrix_b);
}

TEST_F(matrix3_test, op_inequal)
{
	ASSERT_FALSE(matrix_a != matrix_a1);
	ASSERT_TRUE(matrix_a != matrix_b);
}
