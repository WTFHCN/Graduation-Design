#include <iostream>
#include <cmath>
#include <vector>
#include <util.h>
class PublishAlgorithm
{
protected:
    double e;
    double sum;
    std::vector<double> a;

public:
    virtual double add(double val) { throw "Please select a mechanism"; }
    PublishAlgorithm(){};
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

/*
Using the Simple Counting Mechanism II as a building block,
we describe the Two-Level Counting Mechanism. The idea is
that when items from the stream come, we group them in
contiguous blocks of size B, Within a block, we run the Simple
Counting Mechanism II. On top of that, we run another Simple
Counting Mechanism II, treating each block as a single element.
*/
class algorithm3 : public PublishAlgorithm
{
private:
    int N, B, q, r, cnt;
    std::vector<double> block;

public:
    algorithm3(double _e, int n) : PublishAlgorithm(_e), N(n), B(sqrt(n)), q(0), r(0), cnt(0), block(n / sqrt(n)) { a = std::vector<double>(n); }
    virtual double add(double val)
    {
        a[cnt++] = val;
        block[q] += val;
        q += (r + 1 == B);
        r = (r + 1) % B;
        sum = 0;
        for (int i = 0; i < q; i++)
            sum += block[i] + util::getLap(0, e);
        for (int i = q * B; i < cnt; i++)
            sum += a[i] + util::getLap(0, e);
        return sum;
    }
};

/*
We could extend the idea of the Two-Level Counting
Mechanism to a Multi-level Counting Mechanism, and
compute the optimal number of levels given T, the
upper bound on time. However, we take a better approach
called the Binary Mechanism. The idea is that at any
time t, the counting mechanism internally groups the
items that have arrived to form p-sums of different
sizes. The precise grouping of the items depends on
the binary representation of the number t – hence
the name Binary Mechanism.
*/
class algorithm4 : public PublishAlgorithm
{
private:
    int cnt, N;

public:
    algorithm4(double _e, int n) : PublishAlgorithm(_e), cnt(0), N(n) { a = std::vector<double>(n); }
    double add(double val)
    {
        cnt++;
        for (int i = cnt; i <= N; i += LOWBIT(i))
        {
            a[i - 1] += val;
        }
        sum = 0;
        for (int i = cnt; i; i -= LOWBIT(i))
        {
            sum += a[i - 1] + util::getLap(0, e);
        }
        return sum;
    }
};
int main()
{
    PublishAlgorithm *t[4];
    t[0] = new algorithm1(1);
    t[1] = new algorithm2(1);
    t[2] = new algorithm3(1, 100);
    t[3] = new algorithm4(1, 100);
    // std::cout << t[2]->add(1) << '\n';
    for (int i = 0; i < 100; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            std::cout << t[j]->add(i) << " ";
        }
        std::cout << '\n';
    }
}