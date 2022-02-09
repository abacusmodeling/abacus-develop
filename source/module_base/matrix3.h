#ifndef MatriX3_H
#define MatriX3_H

#ifdef _MCD_CHECK
#include "../src_parallel/mcd.h"
#endif

#include "matrix.h"
#include "vector3.h"

namespace ModuleBase
{

class Matrix3
{

  public:
    /**
     * @brief element e_ij: i_row, j_column
     *
     */
    double e11, e12, e13, e21, e22, e23, e31, e32, e33;

    /**
     * @brief Construct a new Matrix 3 object
     * to Identity matrix
     *
     */
    Matrix3()
    {
        Identity();
    }

    /**
     * @brief Construct a new Matrix 3 object
     *
     * @param r11 element r_ij: i_row, j_column
     * @param r12
     * @param r13
     * @param r21
     * @param r22
     * @param r23
     * @param r31
     * @param r32
     * @param r33
     */
    Matrix3(const double &r11,
            const double &r12,
            const double &r13,
            const double &r21,
            const double &r22,
            const double &r23,
            const double &r31,
            const double &r32,
            const double &r33);

    /**
     * @brief Set a 3x3 matrix to identity matrix
     *
     */
    void Identity(void);

    /**
     * @brief Set all elements of a 3x3 matrix to zero
     *
     */
    void Zero(void);

    /**
     * @brief Calculate the determinant of a 3x3 matrix
     *
     * @return double
     */
    double Det(void) const;

    /**
     * @brief Transpose a 3x3 matrix
     *
     * @return Matrix3
     */
    Matrix3 Transpose(void) const;

    /**
     * @brief Inverse a 3x3 matrix
     *
     * @return Matrix3
     */
    Matrix3 Inverse(void) const;

    /**
     * @brief Overload operator "=" for 3x3 matrices
     * For example, assign mb = ma
     *
     * @param m
     * @return Matrix3&
     */
    Matrix3 &operator=(const Matrix3 &m);

    /**
     * @brief Overload operator "+=" for 3x3 matrices
     * For example, mb += ma
     *
     * @param m
     * @return Matrix3&
     */
    Matrix3 &operator+=(const Matrix3 &m);

    /**
     * @brief Overload operator "-=" for 3x3 matrices
     * For example, mb -= ma
     *
     * @param m
     * @return Matrix3&
     */
    Matrix3 &operator-=(const Matrix3 &m);

    /**
     * @brief Overload operator "*=" for 3x3 matrix and
     * a scalar
     * For example, mb *= 3.0
     *
     * @param s The scalar
     * @return Matrix3&
     */
    Matrix3 &operator*=(const double &s);

    /**
     * @brief Overload operator "/=" for 3x3 matrix and
     * a scalar
     * For example, mb /= 3.0
     *
     * @param s The scalar
     * @return Matrix3&
     */
    Matrix3 &operator/=(const double &s);

    /**
     * @brief Print a 3x3 matrix on screening
     *
     */
    void print(void) const;

    /**
     * @brief Change the form of a 3x3 matrix from that of
     * class Matrix3 to that of class matrix
     *
     * @return ModuleBase::matrix
     */
    ModuleBase::matrix to_matrix(void) const;
};

/**
 * @brief Overload operator "+" for two 3x3 matrices m1 and m2
 * i.e. m1+m2
 *
 * @param m1
 * @param m2
 * @return Matrix3
 */
Matrix3 operator+(const Matrix3 &m1, const Matrix3 &m2);

/**
 * @brief Overload operator "-" for two 3x3 matrices m1 and m2,
 * i.e. m1-m2
 *
 * @param m1
 * @param m2
 * @return Matrix3
 */
Matrix3 operator-(const Matrix3 &m1, const Matrix3 &m2);

/**
 * @brief Overload operator "/" for a (3x3 matrix)/(scalar)
 * i.e. m/s
 *
 * @param m The 3x3 matrix
 * @param s The scalar
 * @return Matrix3
 */
Matrix3 operator/(const Matrix3 &m, const double &s);

/**
 * @brief Overload operator "*" for two 3x3 matrices m1 and m2
 * i.e. m1*m2
 *
 * @param m1
 * @param m2
 * @return Matrix3
 */
Matrix3 operator*(const Matrix3 &m1, const Matrix3 &m2);

/**
 * @brief Overload operator "*" for (3x3 matrix)*(scalar)
 * i.e. m*s
 *
 * @param m The 3x3 matrix
 * @param s The scalar
 * @return Matrix3
 */
Matrix3 operator*(const Matrix3 &m, const double &s);

/**
 * @brief Overload operator "*" for (scalar)*(3x3 matrix)
 * i.e. s*m
 *
 * @param s The scalar
 * @param m The 3x3 matrix
 * @return Matrix3
 */
Matrix3 operator*(const double &s, const Matrix3 &m);

/**
 * @brief Over load operator "*" for (3x3 matrix)*(Vector3)
 *
 * @tparam T
 * @param m The 3x3 matrix
 * @param u The vector with 3 elements
 * @return ModuleBase::Vector3<double>
 * @author Peize Lin
 */
template <typename T> ModuleBase::Vector3<double> operator*(const Matrix3 &m, const ModuleBase::Vector3<T> &u);

/**
 * @brief Over load operator "*" for (Vector3)*(3x3 matrix)
 *
 * @tparam T
 * @param u The vector with 3 elements
 * @param m The 3x3 matrix
 * @return ModuleBase::Vector3<double>
 */
template <typename T> ModuleBase::Vector3<double> operator*(const ModuleBase::Vector3<T> &u, const Matrix3 &m);

/**
 * @brief Overload operator "==" to assert
 * the equality between two 3x3 matrices
 *
 * @param m1
 * @param m2
 * @return true
 * @return false
 */
bool operator==(const Matrix3 &m1, const Matrix3 &m2);

/**
 * @brief Overload operator "!=" to assert
 * the inequality between two 3x3 matrices
 *
 * @param m1
 * @param m2
 * @return true
 * @return false
 */
bool operator!=(const Matrix3 &m1, const Matrix3 &m2); // whethor m1 != m2

// m*u
template <typename T> ModuleBase::Vector3<double> operator*(const Matrix3 &m, const ModuleBase::Vector3<T> &u)
{
    return ModuleBase::Vector3<double>(m.e11 * u.x + m.e12 * u.y + m.e13 * u.z,
                                       m.e21 * u.x + m.e22 * u.y + m.e23 * u.z,
                                       m.e31 * u.x + m.e32 * u.y + m.e33 * u.z);
}

// u*m
template <typename T> ModuleBase::Vector3<double> operator*(const ModuleBase::Vector3<T> &u, const Matrix3 &m)
{
    return ModuleBase::Vector3<double>(u.x * m.e11 + u.y * m.e21 + u.z * m.e31,
                                       u.x * m.e12 + u.y * m.e22 + u.z * m.e32,
                                       u.x * m.e13 + u.y * m.e23 + u.z * m.e33);
}

} // namespace ModuleBase

#endif // MATRIX3_H
