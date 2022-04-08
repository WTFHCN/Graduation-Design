#include <iostream>
#include <cmath>
#include <vector>
#include "util.hpp"


class MakeCoefficient
{
private:
    std::vector<double> err;
    std::vector<double> lamda;
    int N, M;

public:
    MakeCoefficient(int m) : M(m), err(1 << m), lamda(1 << m)
    {

        int n = std::pow(2, m);

        err[1] = 1;

        for (int i = 2; i < err.size(); i++)
        {
            double a = pow(err[i - 1], 1 / 3.0);
            double b = pow(1.0 * (1 << (i - 1)), 1 / 3.0);
            err[i] = (a + b) * (a + b) * (a + b) + err[i - 1];
        }

        lamda[1] = 0;
        for (int i = 2; i < err.size(); i++)
        {
            lamda[i] = pow(err[i - 1], 1 / 3.0) / (pow(err[i - 1], 1 / 3.0) + pow(1.0 * (1 << (i - 1)), 1 / 3.0));
        }
    }

    double getCof(int k)
    {
        if (k < 1 || k >= N)
        {
            return 0;
        }
        double ans = 1;
        for (int i = M; i; i--)
        {
            if (k < N / 2)
            {
                ans *= lamda[i];
            }
            else if (k == N / 2)
            {
                ans *= (1 - lamda[i]);
                break;
            }
            else
            {
                k -= N / 2;
            }
            N /= 2;
        }
        return ans;
    }
    double getErr()
    {
        return err[M];
    }
};
class FDA
{
private:
    std::vector<double> a;
    std::vector<double> NoiCn;
    MakeCoefficient coefficient;
    double e;
    int N, M;

public:
    FDA(int m, int e) : M(m), a(1 << M), NoiCn(1 << M), N(1 << M), e(e), coefficient(m){};

    void addVal(int x, double val)
    {
        for (int i = x; i <= N; i += LOWBIT(i))
            a[i - 1] += val;
        NoiCn[x] = coefficient.getCof(x) * a[x] + util::getLap(0, e);
    }
    double getReal(int x)
    {
        double ans = 0;
        for (int i = x; i; i -= LOWBIT(i))
            ans += a[i - 1];
        return ans;
    }
    double query(int x)
    {
        double ans = 0;
        for (int i = x; i; i -= LOWBIT(i))
            ans += NoiCn[i] / coefficient.getCof(i);
        return ans;
    }
};

int main()
{
    int n;
    std::cout << "请输入数据容量参数n(容量为2^n-1): \n";
    std::cin >> n;
    double e;
    std::cout << "请输入隐私参数e:\n";
    std::cin >> e;

    FDA test(n, e);

    std::cout << "请输入要添加的数据:\n";
    for (int i = 1; i < (1 << n); i++)
    {
        double val;
        std::cin >> val;
        test.addVal(i, val);
        std::cout << "第" << i << "次添加数据后，结果是" << test.query(i) << "，真实值为" << test.getReal(i) << '\n';
    }
}
