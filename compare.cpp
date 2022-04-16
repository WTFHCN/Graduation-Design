#include <iostream>
#include <cmath>
#include <vector>
#include <chrono>
#include <random>
#include <time.h>
#include <optional>
#include "util.h"
#include "data.h"
std::vector<std::string> algorithmName = {
    "reality",
    "Simple Counting Mechanism I",
    "Simple Counting Mechanism II",
    "Two-Level Counting Mechanism",
    "Binary Counting Mechanism",
    "Hybrid Mechanism",
};

class PublishAlgorithm
{
protected:
    double e;
    double sum;
    std::vector<double> a;

public:
    virtual double add(double val) { throw "Please select a mechanism"; }
    PublishAlgorithm(){};
    PublishAlgorithm(double _e) : e(_e), sum(0){};
    int getNum() { return a.size(); }
};
class algorithm0 : public PublishAlgorithm
{
public:
    using PublishAlgorithm::PublishAlgorithm;
    virtual double add(double val)
    {
        a.push_back(val);
        sum += val;
        return sum;
    }
};

class algorithm1 : public PublishAlgorithm
{
public:
    using PublishAlgorithm::PublishAlgorithm;
    virtual double add(double val)
    {
        a.push_back(val);
        sum += val;
        return sum + util::getLap(0, 1 / (e / a.size()));
    }
};

class algorithm2 : public PublishAlgorithm
{
public:
    using PublishAlgorithm::PublishAlgorithm;
    double add(double val)
    {
        a.push_back(val);
        sum += val + util::getLap(0, 1 / e);
        return sum;
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
            sum += block[i] + util::getLap(0, 1 / e);
        for (int i = q * B; i < cnt; i++)
            sum += a[i] + util::getLap(0, 1 / e);
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
    std::vector<double> psum;

public:
    algorithm4(double _e, int n) : PublishAlgorithm(_e), cnt(0), N(n), psum(n) {}
    double add(double val)
    {
        a.push_back(val);
        cnt++;
        for (int i = cnt; i <= N; i += LOWBIT(i))
        {
            psum[i - 1] += val;
        }
        sum = 0;
        for (int i = cnt; i; i -= LOWBIT(i))
        {
            sum += psum[i - 1] + util::getLap(0, 1 / (e / log2(N)));
        }
        return sum;
    }
};

class algorithm5 : public PublishAlgorithm
{
private:
    int cnt, N;
    std::vector<std::vector<double>> block;
    std::vector<double> psum;

    int len;
    void pushup(int pos)
    {
        psum[pos] = psum[pos << 1] + psum[pos << 1 | 1];
    }

    void addInTree(int k, int pos, int L, int R)
    {

        if (R - L < len)
        {
            int t = R - L + 1;
            int B = sqrt(t);
            if (!block[pos].size())
                block[pos].resize(t / B + (t % B != 0));
            int q = (k - L) / B;
            block[pos][q] += a[k];
            psum[pos] += a[k];
            return;
        }
        int mid = (L + R) / 2;
        if (k <= mid)
            addInTree(k, pos << 1, L, mid);
        else
            addInTree(k, pos << 1 | 1, mid + 1, R);
        pushup(pos);
    }
    double queryInTree(int ql, int qr, int pos, int L, int R)
    {

        if (R - L < len)
        {
            int t = R - L + 1;
            int B = sqrt(t);
            if (!block[pos].size())
                block[pos].resize(t / B + (t % B != 0));

            int fl = (ql - L) / B;
            int fr = (qr - L) / B;
            double res = 0;
            if (fr - fl <= 1)
            {
                for (int i = ql; i <= qr; i++)
                    res += a[i] + util::getLap(0, 1 / (e / log2(N / len)));
            }
            else
            {
                for (int i = fl + 1; i < fr; i++)
                    res += block[pos][i] + util::getLap(0, 1 / (e / log2(N / len)));
                for (int i = fr * B + L; i <= qr; i++)
                    res += a[i] + util::getLap(0, 1 / (e / log2(N / len)));
                for (int i = ql; i < (fl + 1) * B + L; i++)
                    res += a[i] + util::getLap(0, 1 / (e / log2(N / len)));
            }

            return res;
        }
        else if (ql <= L && R <= qr)
        {
            return psum[pos] + util::getLap(0, 1 / (e / log2(N / len)));
        }
        int mid = (L + R) / 2;
        if (qr <= mid)
            return queryInTree(ql, qr, pos << 1, L, mid);
        else if (ql > mid)
            return queryInTree(ql, qr, pos << 1 | 1, mid + 1, R);
        else
            return (queryInTree(ql, mid, pos << 1, L, mid) + queryInTree(mid + 1, qr, pos << 1 | 1, mid + 1, R));
    }

public:
    algorithm5(double _e, int n, int _len) : PublishAlgorithm(_e), len(_len), cnt(0), N(n), block(n << 2), psum(n << 2) {}

    double add(double val)
    {
        a.push_back(val);
        addInTree(a.size() - 1, 1, 0, N - 1);
        return queryInTree(0, a.size() - 1, 1, 0, N - 1);
    }
};

int main()
{
    int T = 256;
    int NUM = 6;
    double e = 1;
    std::string csvName = "test.csv";
    PublishAlgorithm *t[6];
    t[0] = new algorithm0(e);
    t[1] = new algorithm1(e);
    t[2] = new algorithm2(e);
    t[3] = new algorithm3(e, T);
    t[4] = new algorithm4(e, T);
    t[5] = new algorithm5(e, T, 32);

    std::default_random_engine random(time(NULL));
    std::uniform_int_distribution<int> dis1(0, 100);
    std::vector<std::vector<double>> testData(NUM, std::vector<double>(T));
    for (int i = 0; i < T; i++)
    {
        double c = dis1(random);
        for (int j = 0; j < NUM; j++)
        {
            testData[j][i] = t[j]->add(c);
        }
    }
    wirteCsv(csvName, algorithmName, testData);
}