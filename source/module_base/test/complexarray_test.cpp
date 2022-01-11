#include"../complexarray.h"
#include"gtest/gtest.h"

/************************************************
*  unit test of class ComplexArray
***********************************************/

/**
* - Tested Function:
*   - operator "=":
*       - assign a complex to all elements of a ComplexArray
*       - assign a ComplexArray to another ComplexArray
*
*   - operator "+":
*       - one ComplexArray plus another ComplexArray that has
*         the same dimension.
*       - throw error when one ComplexArray plus another ComplexArray 
*         that has different dimension.
*
*   - operator "+=":
*       - one ComplexArray plus another ComplexArray, and assign the 
*         value to first ComplexArray
*       - throw error when one ComplexArray plus another ComplexArray
*         that has different dimension.
*
*   - operator "-":
*       - one ComplexArray minus another ComplexArray that has
*         the same dimension.
*       - throw error when one ComplexArray minus another ComplexArray
*         that has different dimension.
*
*   - operator "-=":
*       - one ComplexArray minus another ComplexArray, and assign the
*         value to first ComplexArray
*       - throw error when one ComplexArray minus another ComplexArray
*         that has different dimension.
*
*   - operator "*":
*       - one ComplexArray is multiplied by a double
*       - one ComplexArray is multiplied by a complex
*       - one ComplexArray is multiplied by another ComplexArray that has same dimension
*       - throw error when one ComplexArray is miltiplied by another ComplexArray
*         that has different dimension.
*
*   - operator "*=":
*       similar as "*"
*
*   - oprator "()":
*       - access the element
*/


class ComplexArray_test : public testing::Test
{
    protected:
        ModuleBase::ComplexArray a2,a4,b2,b4,c2,c4,d2;
        ModuleBase::ComplexArray a2_plus_b2,a2_minus_b2,a2_mult_b2,a2_mult_3;
        std::complex<double> com1 {1.0,2.0};
        std::complex<double> com2 {3.0,4.0};

        void SetUp()
        {
            a2 = ModuleBase::ComplexArray(2,1,1,1);  // if define this class as a matrix
            b2 = ModuleBase::ComplexArray(2,1,1,1);  // of 4 dimesions,
            c2 = ModuleBase::ComplexArray(2,1,1,1);  // is a2 +/-/* d2 allowed or not???
            d2 = ModuleBase::ComplexArray(1,2,1,1);  // it does not matter, this situation will not appear in ABACUS
            a2_plus_b2 = ModuleBase::ComplexArray(2,1,1,1);
            a2_minus_b2 = ModuleBase::ComplexArray(2,1,1,1);
            a2_mult_b2 = ModuleBase::ComplexArray(2,1,1,1);
            a2_mult_3 = ModuleBase::ComplexArray(2,1,1,1);
            a4 = ModuleBase::ComplexArray(2,2,1,1);
            b4 = ModuleBase::ComplexArray(2,2,1,1);
            c4 = ModuleBase::ComplexArray(2,2,1,1);
            a2 = com1;
            b2 = com2;
            a4 = com1;
            b4 = com2;
            d2 = com1;
            a2_plus_b2 = com1 + com2;
            a2_minus_b2 = com1 - com2;
            a2_mult_b2 = com1 * com2;
            a2_mult_3 = com1 * 3.0;
        }

};

TEST_F(ComplexArray_test,operator_equal)
{   
    c2 = a2; //c2 is not assigned
    b2 = a2; //b2 is arrigned by a complex
    EXPECT_EQ(c2,a2);
    EXPECT_EQ(b2,a2);
    EXPECT_NE(c2,a4);
    EXPECT_DEATH(c2=a4,"");
}

TEST_F(ComplexArray_test,operator_plus)
{
    c2 = a2 + b2;
    EXPECT_EQ(c2,a2_plus_b2);
    EXPECT_DEATH(a2+a4,"");
    //EXPECT_DEATH(a2+d2,"");
}

TEST_F(ComplexArray_test,operator_plus_equal)
{
    a2 += b2;
    EXPECT_EQ(a2,a2_plus_b2);
    EXPECT_DEATH(a2+=a4,"");
    //EXPECT_DEATH(a2+=d2,"");
}

TEST_F(ComplexArray_test,operator_minus)
{
    c2 = a2 - b2;
    EXPECT_EQ(c2,a2_minus_b2);
    EXPECT_DEATH(a2-a4,"");
    //EXPECT_DEATH(a2-d2,"");
}

TEST_F(ComplexArray_test,operator_minus_equal)
{      
    a2 -= b2;
    EXPECT_EQ(a2,a2_minus_b2);
    EXPECT_DEATH(a2-=a4,"");
    //EXPECT_DEATH(a2-=d2,"");
}


TEST_F(ComplexArray_test,operator_multiply)
{
    c2 = a2 * com2;
    EXPECT_EQ(c2,a2_mult_b2);

    c2 = a2*3.0;
    EXPECT_EQ(c2,a2_mult_3);
   // EXPECT_ANY_THROW(a2*a4);
   // EXPECT_ANY_THROW(a2*d2);
}

TEST_F(ComplexArray_test,operator_multiply_equal)
{
    c2 = a2;
    c2 *= b2;
    EXPECT_EQ(c2,a2_mult_b2);

    c2 = a2;
    c2 *= 3.0;
    EXPECT_EQ(c2,a2_mult_3);

    c2 = a2;
    c2 *= com2;
    EXPECT_EQ(c2,a2_mult_b2);

    EXPECT_DEATH(a2*=a4,"");
 //   EXPECT_DEATH(a2*=d2,"");
}

TEST_F(ComplexArray_test,operator_parentheses)
{
    c2 = a2;
    EXPECT_EQ(c2(0,0,0,0), com1);

    c2(1,0,0,0) = com2;
    EXPECT_NE(c2,a2);

    EXPECT_DEATH(a2(1,1,1,0),"");
}

