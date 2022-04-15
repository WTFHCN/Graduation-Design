#include <iostream>
#include <cmath>
#include <vector>
#include "util.h"
#include "data.h"
class MakeCoefficient
{
private:
    std::vector<double> err;
    std::vector<double> lamda;
    int N, M;

public:
    MakeCoefficient(int m) : M(m), err(1 << m), lamda(1 << m), N(1 << m)
    {
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
        int K = N;
        for (int i = M; i; i--)
        {
            if (k < K / 2)
            {
                ans *= lamda[i];
            }
            else if (k == K / 2)
            {
                ans *= (1 - lamda[i]);
                break;
            }
            else
            {
                k -= K / 2;
            }
            K /= 2;
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
    int N, M, cnt;

public:
    FDA(int m, int e) : cnt(0), M(m), a((1 << m) + 1), NoiCn((1 << m) + 1), N(1 << m), e(e), coefficient(m){};

    void addVal(double val)
    {
        cnt++;
        for (int i = cnt; i <= N; i += LOWBIT(i))
            a[i] += val;
        NoiCn[cnt] = coefficient.getCof(cnt) * a[cnt] + util::getLap(0, e);
    }
    double getReal()
    {
        double ans = 0;
        for (int i = cnt; i; i -= LOWBIT(i))
            ans += a[i];
        return ans;
    }
    double query()
    {
        double ans = 0;
        for (int i = cnt; i; i -= LOWBIT(i))
            ans += NoiCn[i] / coefficient.getCof(i);
        return ans;
    }
};

int main()
{
    int m;
    std::cout << "请输入数据容量参数m(容量为2^m-1): \n";
    std::cin >> m;
    double e;
    std::cout << "请输入隐私参数e:\n";
    std::cin >> e;

    FDA test(m, e);

    std::cout << "请输入要添加的数据:\n";
    for (int i = 1; i < (1 << m); i++)
    {
        double val;
        std::cin >> val;
        test.addVal(val);
        std::cout << "第" << i << "次添加数据后，结果是" << test.query() << "，真实值为" << test.getReal() << '\n';
    }
}
