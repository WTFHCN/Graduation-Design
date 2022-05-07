#include <Eigen/Dense>
#include <Eigen/Core>
#include "FDA.hpp"
#include <gtest/gtest.h>

using namespace std;
using namespace Eigen;

double error(MatrixXd B, MatrixXd L)
{
    return B.squaredNorm() * pow((L.cwiseAbs().colwise().sum().maxCoeff()), 2);
}
TEST(testFDA, FDA_m_3)
{
    int m = 3;
    MakeCoefficient test(m);
    int n = (1 << m) - 1;
    MatrixXd B = MatrixXd::Zero(n, n);
    MatrixXd L = B;
    MatrixXd P = B;
    for (int i = 1; i <= n; i++)
    {
        for (int j = i; j <= n; j += LOWBIT(j))
        {
            L(j - 1, i - 1) = 1;
        }
    }
    for (int i = 1; i <= n; i++)
    {
        for (int j = i; j; j -= LOWBIT(j))
        {
            B(i - 1, j - 1) = 1;
        }
    }
    for (int i = 1; i <= n; i++)
    {
        P(i - 1, i - 1) = test.getCof(i);
    }

    auto BB = B * P.inverse();
    auto LL = P * L;

    EXPECT_LE(error(BB, LL), error(B, L));
}
int main(int argc, char *argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
