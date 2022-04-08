#include <iostream>
#include <cmath>
#include <vector>

#include <random>
#include <time.h>
#define LOWBIT(x) (x & -x)
namespace util
{
    std::default_random_engine random(time(NULL));
    std::uniform_real_distribution<double> randOut(0.0, 1.0);
    double getLap(double u, double b)
    {
        double x = randOut(random);
        if (x < 0.5)
            return b * std::log(2 * x) + u;
        else
            return u - b * std::log(2 - 2 * x);
    }
}
class PublishAlgorithm
{
protected:
    double e;
    double sum;
    static std::vector<int> a;

public:
    virtual double add(double val) {}
    PublishAlgorithm(double _e) : e(_e){};
};

class algorithm1 : public PublishAlgorithm
{
public:
    using PublishAlgorithm::PublishAlgorithm;
    virtual double add(double val)
    {
        a.push_back(val);
        sum += val;
        return sum + util::getLap(0, e);
    }
};

class algorithm2 : public PublishAlgorithm
{
public:
    using PublishAlgorithm::PublishAlgorithm;
    double add(double val)
    {
        a.push_back(val);
        sum += val;
        return sum + util::getLap(0, e);
    }
};

class algorithm3 : public PublishAlgorithm
{
private:
    int N, B, q, r, cnt;

    std::vector<int> block;

public:
    algorithm3(double _e, int n) : PublishAlgorithm(e), N(n), B(sqrt(n)), q(0), r(0), cnt(0), block(n / sqrt(n)) { a = std::vector<int>(n); }
    virtual double add(double val)
    {
        a[cnt++] = val;
        block[q] += a[cnt];
        r++;
        q += (r == B);
        sum = 0;
        for (int i = 0; i < q; i++)
            sum += block[i] + util::getLap(0, e);
        for (int i = q * r; i < cnt; i++)
            sum += a[i] + util::getLap(0, e);
        return sum;
    }
};

class algorithm4
{
private:
    double E, sum;
    int cnt, N;
    std::vector<int> a;

public:
    algorithm4(double e, int n) : E(e), N(n), a(n), cnt(0) {}

    double add(double val)
    {
        for (int i = cnt; i <= N; i += LOWBIT(i))
        {
            a[i - 1] += val;
        }
        sum = 0;
        for (int i = cnt; i; i -= LOWBIT(i))
        {
            sum += a[i - 1] + util::getLap(0, E);
        }
        return sum;
    }
};
int main()
{
    PublishAlgorithm *t;
    // t[0] = new algorithm1(1);
    // t[1] = new algorithm2(1);
    // t[2] = new algorithm3(1, 10);
    t = new algorithm1(1);
    std::cout << t->add(1) << '\n';
    // for (int i = 0; i < 10; i++)
    // {
    //     for (int j = 0; j < 3; j++)
    //     {
    //         std::cout << t[i]->add(i) << " ";
    //     }
    //     std::cout << '\n';
    // }
}