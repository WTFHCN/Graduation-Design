#include <iostream>
#include <cmath>
#include <vector>

class MakeCoefficient
{
private:
    std::vector<double> err;
    std::vector<double> lamda;
    int N, M;

public:
    MakeCoefficient(int m) : M(m), err(m + 1), lamda(m + 1), N(1 << m)
    {
        err[1] = 1;
        for (int i = 2; i < err.size(); i++)
        {
            double c = pow(err[i - 1], 1 / 3.0);
            double b = pow(1.0 * (1 << (i - 1)), 1 / 3.0);
            err[i] = (c + b) * (c + b) * (c + b) + err[i - 1];
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

// int main()
// {
//     int m = 3;
//     double e = 1;
//     FDA test(m, e);

//     std::cout << "请输入要添加的数据:\n";
//     for (int i = 1; i < (1 << m); i++)
//     {
//         double val;
//         std::cin >> val;
//         test.addVal(val);
//         std::cout << "第" << i << "次添加数据后，结果是" << test.query() << "，真实值为" << test.getReal() << '\n';
//     }
// }
