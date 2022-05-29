#include <iostream>
#include <cmath>
#include <vector>
#include <chrono>
#include <random>
#include <time.h>
#include <optional>
#include "data.h"
#include "FDA.hpp"
const int NUM = 7;
std::string const csvName = "census.csv";
std::vector<std::string> algorithmName = {
    "reality",
    "Simple Counting Mechanism I",
    "Simple Counting Mechanism II",
    "TM",
    "BM",
    "Hybrid Mechanism",
    "FDA",
};
std::vector<std::string> lenName = {
    "reality",
    "len=25",
    "len=50",
    "len=100",
    "len=200",
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
        return sum + util::GetLap(0, 1 / (e / a.size()));
    }
};

class algorithm2 : public PublishAlgorithm
{
public:
    using PublishAlgorithm::PublishAlgorithm;
    double add(double val)
    {
        a.push_back(val);
        sum += val + util::GetLap(0, 1 / e);
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
    algorithm3(double _e, int n) : PublishAlgorithm(_e), N(n), B(sqrt(n)), q(0), r(0), cnt(0), block(n / sqrt(n) + 2) { a = std::vector<double>(n); }
    virtual double add(double val)
    {
        a[cnt++] = val;
        block[q] += val;
        q += (r + 1 == B);
        r = (r + 1) % B;
        sum = 0;
        for (int i = 0; i < q; i++)
            sum += block[i] + util::GetLap(0, 1 / e);
        for (int i = q * B; i < cnt; i++)
            sum += a[i] + util::GetLap(0, 1 / e);
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
            sum += psum[i - 1] + util::GetLap(0, 1 / (e / log2(N)));
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
                block[pos].resize(t / B + 2);
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
                block[pos].resize(t / B + 2);

            int fl = (ql - L) / B;
            int fr = (qr - L) / B;
            double res = 0;
            if (fr - fl <= 1)
            {
                for (int i = ql; i <= qr; i++)
                    res += a[i] + util::GetLap(0, 1 / (e / log2(N / len)));
            }
            else
            {
                for (int i = fl + 1; i < fr; i++)
                    res += block[pos][i] + util::GetLap(0, 1 / (e / log2(N / len)));
                for (int i = fr * B + L; i <= qr; i++)
                    res += a[i] + util::GetLap(0, 1 / (e / log2(N / len)));
                for (int i = ql; i < (fl + 1) * B + L; i++)
                    res += a[i] + util::GetLap(0, 1 / (e / log2(N / len)));
            }

            return res;
        }
        else if (ql <= L && R <= qr)
        {
            return psum[pos] + util::GetLap(0, 1 / (e / log2(N / len)));
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
class algorithm6 : public PublishAlgorithm
{
private:
    std::vector<double> c;
    std::vector<double> NoiCn;
    MakeCoefficient coefficient;
    int N, M, cnt;

public:
    algorithm6(double _e, int m) : PublishAlgorithm(_e), cnt(0), M(m), c((1 << m) + 1), NoiCn((1 << m) + 1), N(1 << m), coefficient(m){};

    double add(double val)
    {
        a.push_back(val);
        cnt++;
        for (int i = cnt; i <= N; i += LOWBIT(i))
            c[i] += val;
        NoiCn[cnt] = coefficient.getCof(cnt) * c[cnt] + util::GetLap(0, 1 / e);
        return query();
    }
    double getReal()
    {
        double ans = 0;
        for (int i = cnt; i; i -= LOWBIT(i))
            ans += c[i];
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
void CensusIncomeTest(int M, double e, std::vector<std::vector<double>> &errorList, std::vector<std::vector<double>> &timeList)
{
    int T = (1 << M) - 1;
    auto CensusList = ReadCensusCsv();
    assert(int(CensusList.size()) >= T);
    std::vector<std::vector<double>>
        testData(NUM, std::vector<double>(T));

    int tt = 10;

    std::vector<double> error(NUM);
    std::vector<double> timeSpent(NUM);
    for (int t = 0; t < tt; t++)
    {
        std::vector<PublishAlgorithm *> algorithmMechanism(NUM);
        algorithmMechanism[0] = new algorithm0(e);
        algorithmMechanism[1] = new algorithm1(e);
        algorithmMechanism[2] = new algorithm2(e);
        algorithmMechanism[3] = new algorithm3(e, T);
        algorithmMechanism[4] = new algorithm4(e, T);
        algorithmMechanism[5] = new algorithm5(e, T, 100);
        algorithmMechanism[6] = new algorithm6(e, M);
        for (int i = 0; i < T; i++)
        {
            for (int j = 0; j < NUM; j++)
            {
                auto t = clock();
                testData[j][i] = algorithmMechanism[j]->add(CensusList[i].income);
                timeSpent[j] += (double)(clock() - t) / CLOCKS_PER_SEC;
            }
        }
        for (int j = 1; j < NUM; j++)
        {
            error[j] += util::CalcError(testData[0], testData[j]);
        }
    }

    for (int j = 1; j < NUM; j++)
    {
        errorList[j][M] = error[j] / tt;
        timeList[j][M] = timeSpent[j] / tt;
    }
}
void RetailTest(int M, double e, std::vector<std::vector<double>> &errorList, std::vector<std::vector<double>> &timeList)
{
    int T = (1 << M) - 1;
    auto RetailList = ReadRetailCsv();
    assert(int(RetailList.size()) >= T);
    std::vector<std::vector<double>>
        testData(NUM, std::vector<double>(T));

    int tt = 10;

    std::vector<double> error(NUM);
    std::vector<double> timeSpent(NUM);
    for (int t = 0; t < tt; t++)
    {
        std::vector<PublishAlgorithm *> algorithmMechanism(NUM);
        algorithmMechanism[0] = new algorithm0(e);
        algorithmMechanism[1] = new algorithm1(e);
        algorithmMechanism[2] = new algorithm2(e);
        algorithmMechanism[3] = new algorithm3(e, T);
        algorithmMechanism[4] = new algorithm4(e, T);
        algorithmMechanism[5] = new algorithm5(e, T, 100);
        algorithmMechanism[6] = new algorithm6(e, M);
        for (int i = 0; i < T; i++)
        {
            for (int j = 0; j < NUM; j++)
            {
                auto t = clock();
                testData[j][i] = algorithmMechanism[j]->add(RetailList[i].UnitPrice);
                timeSpent[j] += (double)(clock() - t) / CLOCKS_PER_SEC;
            }
        }
        for (int j = 1; j < NUM; j++)
        {
            error[j] += util::CalcError(testData[0], testData[j]);
        }
    }

    for (int j = 1; j < NUM; j++)
    {
        errorList[j][M] = error[j] / tt;
        timeList[j][M] = timeSpent[j] / tt;
    }
}
void lenTest(int M, double e, std::vector<std::vector<double>> &errorList)
{
    int T = (1 << M) - 1;
    int N = lenName.size();

    auto RetailList = ReadRetailCsv();
    assert(int(RetailList.size()) >= T);
    std::vector<std::vector<double>>
        testData(N, std::vector<double>(T));

    int tt = 10;

    std::vector<double> error(N);
    for (int t = 0; t < tt; t++)
    {
        std::vector<PublishAlgorithm *> algorithmMechanism(N);
        algorithmMechanism[0] = new algorithm0(e);
        algorithmMechanism[1] = new algorithm5(e, T, 25);
        algorithmMechanism[2] = new algorithm5(e, T, 50);
        algorithmMechanism[3] = new algorithm5(e, T, 100);
        algorithmMechanism[4] = new algorithm5(e, T, 200);
        for (int i = 0; i < T; i++)
        {
            for (int j = 0; j < N; j++)
            {
                testData[j][i] = algorithmMechanism[j]->add(RetailList[i].UnitPrice);
            }
        }
        for (int j = 1; j < N; j++)
        {
            error[j] += util::CalcError(testData[0], testData[j]);
        }
    }

    for (int j = 1; j < N; j++)
    {
        errorList[j][M] = error[j] / tt;
    }
}
void test1(double e)
{
    std::string outputName = std::to_string(int(e * 100)) + "Income.csv";
    int maxM = 14;
    std::vector<std::vector<double>> errorList(NUM, std::vector<double>(maxM));
    std::vector<std::vector<double>> timeList(NUM, std::vector<double>(maxM));
    for (int i = 8; i < maxM; i++)
    {
        CensusIncomeTest(i, e, errorList, timeList);
    }
    WirteCsv(outputName, algorithmName, errorList);
}
void test2(double e)
{
    std::string outputName = std::to_string(int(e * 100)) + "Retail.csv";
    int maxM = 14;
    std::vector<std::vector<double>> errorList(NUM, std::vector<double>(maxM));
    std::vector<std::vector<double>> timeList(NUM, std::vector<double>(maxM));

    for (int i = 8; i < maxM; i++)
    {
        RetailTest(i, e, errorList, timeList);
    }
    WirteCsv(outputName, algorithmName, errorList);
}
void test3(double e)
{
    std::string outputName = std::to_string(int(e * 100)) + "lenRetail.csv";
    int maxM = 14;
    std::vector<std::vector<double>> errorList(NUM, std::vector<double>(maxM));
    for (int i = 8; i < maxM; i++)
    {
        lenTest(i, e, errorList);
    }
    WirteCsv(outputName, lenName, errorList);
}
int main()
{
    // for (double e = 0.01; e <= 1; e *= 10)
    // {
    //     test1(e);
    //     test2(e);
    // }
    test3(0.01);
}