#include <iostream>
#include <ctime>
#include <random>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Core>
#include "util.h"

using namespace std;
using namespace Eigen;
const double eps = 1e-5;
const double inf = 1e10;

static default_random_engine e(time(0));
static normal_distribution<double> normal(0, 1);
double sign(double x)
{
    if (x < 0)
        return -1;
    if (x > 0)
        return 1;
    return 0;
}
double MyRank(MatrixXd x, double err_eps)
{
    FullPivLU<MatrixXd> lu_decomp(x);
    lu_decomp.setThreshold(err_eps);
    return lu_decomp.rank();
}

MatrixXd GetNormalMatrix(int n, int m)
{
    return MatrixXd::Zero(n, m).unaryExpr([&](double dummy)
                                          { return normal(e); });
}

MatrixXd DimDot(MatrixXd x, MatrixXd y, int p = 0)
{
    if (p == 0)
    {
        MatrixXd z(1, x.cols());
        for (int i = 0; i < x.cols(); i++)
        {
            z(0, i) = x.col(i).dot(y.col(i));
        }
        return z;
    }
    else
    {
        MatrixXd z(x.rows(), 1);
        for (int i = 0; i < x.rows(); i++)
        {
            z(i, 0) = x.row(i).dot(y.row(i));
        }
        return z;
    }
}

MatrixXd GenLaplace(int n, int m, double mu, double b)
{
    MatrixXd res(n, m);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            res(i, j) = util::GetLap(mu, b);
    return res;
}

double L1sensitivity(MatrixXd T)
{
    double ans = 0;
    for (int j = 0; j < T.cols(); j++)
    {
        double sum = 0;
        for (int i = 0; i < T.rows(); i++)
            sum += abs(T(i, j));
        ans = max(ans, sum);
    }
    return ans;
}

MatrixXd GenNoiseAns(MatrixXd x, MatrixXd T, double epsilon)
{
    double noiseScale = L1sensitivity(T) / (epsilon);
    MatrixXd noise = GenLaplace(T.rows(), 1, 0, noiseScale);
    MatrixXd xNoise = T * x + noise;
    return xNoise;
}

double ComputeRealErrorMF(MatrixXd W, MatrixXd domain_x, double epsilon, MatrixXd Q, MatrixXd T)
{
    MatrixXd trueAns = W * domain_x;
    MatrixXd noisyX1 = GenNoiseAns(domain_x, T, epsilon);
    MatrixXd myAns = Q * noisyX1;
    double error = pow((trueAns - myAns).norm(), 2);
    return error;
}

double Compute_QT_Eror(MatrixXd W, MatrixXd Q, MatrixXd T, MatrixXd domain_x, double epsilon, int times)
{
    MatrixXd es(1, times);
    for (int i = 0; i < times; i++)
    {
        double e = ComputeRealErrorMF(W, domain_x, epsilon, Q, T);
        es(0, i) = e;
    }
    return es.mean();
}

MatrixXd GenNoiseAns1(MatrixXd x, MatrixXd I, double epsilon, int times = 10)
{
    epsilon = epsilon / times;
    MatrixXd noise = MatrixXd::Zero(x.rows(), x.cols());
    for (int i = 0; i < times; i++)
    {
        double noiseScale = 1 / (epsilon);
        MatrixXd noise_tmp = GenLaplace(I.rows(), 1, 0, noiseScale);
        noise = noise + noise_tmp;
    }
    noise = noise / times;
    MatrixXd x_noise2 = x + noise;
    return x_noise2;
}

double ComputeIRealError(MatrixXd W, MatrixXd domain_x, double epsilon, MatrixXd I)
{
    MatrixXd trueAns = W * domain_x;
    MatrixXd noisyX1 = GenNoiseAns1(domain_x, I, epsilon);
    MatrixXd myAns = W * noisyX1;
    double error = pow((trueAns - myAns).norm(), 2);
    return error;
}
double Compute_I_Eror(MatrixXd W, MatrixXd domain_x, double epsilon, int times)
{
    MatrixXd I = MatrixXd::Identity(W.cols(), W.cols());
    MatrixXd es(1, times);
    for (int i = 0; i < times; i++)
    {
        double e = ComputeIRealError(W, domain_x, epsilon, I);
        es(0, i) = e;
    }
    return es.mean();
}

MatrixXd UpdateQ(MatrixXd W, MatrixXd Q, MatrixXd T, double C, double beta, MatrixXd pi, double WFnorm)
{
    MatrixXd A = beta * W * T.transpose() + pi * T.transpose();
    MatrixXd B = beta * T * T.transpose() + MatrixXd::Identity(T.rows(), T.rows());
    Q = A * (B.inverse());
    return Q;
}

void smallperm(MatrixXd &T)
{
    MatrixXd R = T.unaryExpr([](double val)
                             { return double(abs(val) < eps); });
    T = T + (R.cwiseProduct(GetNormalMatrix(R.rows(), R.cols()))) * 1e-5;
}

pair<double, MatrixXd> ComputeT(MatrixXd W, MatrixXd Q, MatrixXd T, double C, double beta, MatrixXd pi, double WFnorm)
{
    MatrixXd QQ = Q.transpose() * Q;
    MatrixXd TT = T * T.transpose();
    double sig = 1e-6;
    double F = 0.5 * beta * (-2 * (Q.transpose() * W * T.transpose()).trace() + (TT * QQ).trace()) - (T * pi.transpose() * Q).trace();
    F = F + 0.5 * sig * TT.trace();
    MatrixXd G = beta * (QQ * T - Q.transpose() * W) - Q.transpose() * pi;
    G = G + sig * T;
    return make_pair(F, G);
}
VectorXd L1Projection(VectorXd v, double b)
{
    // vector V = vector<double>(v.data(), v.data() + v.size());
    auto u = v;
    sort(u.begin(), u.end(), greater<double>());
    u = u.cwiseAbs();
    assert(u.cols() == 1);
    auto sv = u;
    for (int i = 1; i < sv.rows(); i++)
    {
        sv(i, 0) += sv(i - 1, 0);
    }
    int pos = 0;
    for (int i = 1; i < v.rows(); i++)
    {
        if (u(i, 0) > (sv(i, 0) - b) / (i + 1))
        {
            pos = i;
        }
    }
    double theta = max(0., (sv(pos, 0) - b) / (pos + 1));
    MatrixXd sign_v = v.unaryExpr([&](double x)
                                  { return sign(x); });
    MatrixXd vv = v.array() - theta;
    MatrixXd w = sign_v.cwiseProduct(vv.cwiseMax(1e-8));
    return w;
}
MatrixXd BallProjectionT(MatrixXd &T, double c)
{
    MatrixXd TT(T.rows(), T.cols());
    for (int i = 0; i < T.cols(); i++)
    {
        TT.col(i) = L1Projection(T.col(i), c);
    }
    return TT;
}

MatrixXd UpdateT(MatrixXd W, MatrixXd Q, MatrixXd T, double C, double beta, MatrixXd pi, double WFnorm)
{
    int n = T.rows();
    int m = T.cols();
    double threshold = T.size() * 1e-12;
    int maxIter = 5;
    bool bFlag = 0;
    MatrixXd xxp = MatrixXd::Zero(n, m);
    double alphap = 0, alpha = 1;
    double L = max(1., beta * (Q.transpose() * Q).cwiseAbs().maxCoeff() / 100);
    for (int t = 0; t < maxIter; t++)
    {
        double alpha1 = (alphap - 1) / alpha;
        MatrixXd s = T + alpha1 * xxp;
        auto [f_last, g] = ComputeT(W, Q, s, C, beta, pi, WFnorm);
        MatrixXd xp = T;
        int inner = 0;
        while (1)
        {
            inner++;
            MatrixXd v = s - g / L;
            T = BallProjectionT(v, C);
            double f_curr = ComputeT(W, Q, T, C, beta, pi, WFnorm).first;
            MatrixXd diff = T - s;
            double r_sum1 = 0.5 * DimDot(diff, diff, 1).sum();
            double l_sum1 = f_curr - f_last - DimDot(diff, g, 1).sum();
            if (r_sum1 <= threshold)
            {
                bFlag = 1;
                printf("the gradient step makes little improvement");
                break;
            }
            if (l_sum1 - r_sum1 * L <= 0)
            {
                printf("L smooth");
                break;
            }
            L = max(2 * L, l_sum1 / r_sum1);
            if (inner == 5)
            {
                auto QS = (Q.transpose() * Q).cwiseAbs().colwise().sum();
                L = beta * QS.maxCoeff();
            }
        }

        alphap = alpha;
        alpha = (1 + sqrt(4 * alpha * alpha + 1)) / 2;
        xxp = T - xp;
        if (bFlag)
        {
            printf("\n The program terminates as the gradient step changes the solution very small.");
            break;
        }
    }
    smallperm(T);
    return T;
}
pair<MatrixXd, MatrixXd> SolveApproximatly(MatrixXd W, MatrixXd Q, MatrixXd T, double C, double beta, MatrixXd pi, double WFnorm)
{
    for (int i = 0; i < 5; i++)
    {
        UpdateT(W, Q, T, C, beta, pi, WFnorm);
        UpdateQ(W, Q, T, C, beta, pi, WFnorm);
    }
}

double ComputeExpectedError(MatrixXd Q, double C)
{
    return (Q.transpose() * Q).trace() * C * C;
}

// This programme solve the following optimization problem:
// min 0.5 tr(Q'Q)
// s.t. W = QT
//      |T(:,i)| <= 1,
//      for all i = 1...n
// W: m x n
// Q: m x r
// T: r x n

pair<MatrixXd, MatrixXd> LowRankDP(MatrixXd W, double r = 0.0, double err_eps = 0.01, int max_iter = 100)
{
    FullPivLU<MatrixXd> lu_decomp(W);
    lu_decomp.setThreshold(err_eps);
    if (r == 0.0)
        r = 1.2 * lu_decomp.rank();
    int WiseInit = 1;
    int IncreaseRank = 0;
    int UpdateRule = 1;
    int factor = 5;
    int n = W.rows();
    int m = W.cols();
    double WFnorm = W.squaredNorm();
    MatrixXd Q, T;
    if (WiseInit)
    {
    }
    else
    {
        Q = GetNormalMatrix(n, r);
        T = GetNormalMatrix(r, m);
    }
    printf("n:%ld, m:%ld, r:%ld, r:%ld\n", Q.rows(), T.cols(), T.rows(), Q.cols());
    MatrixXd pi = MatrixXd::Zero(W.rows(), W.cols());
    double beta = 1;
    double C = 1;
    double M = 1e10;
    int increasement = 0; //!!!!!!
    if (IncreaseRank)
    {
        double rr = MyRank(W, 0.01);
        int increase = floor(0.2 * rr);
        if (increasement == 0)
            increasement = 1;
    }
    int flag = 1;
    int inner_iter = 0;
    int outer_iter = 0;
    double last_err = inf;
    double t1 = clock();
    double curr_err = 0.0;
    while (1)
    {
        inner_iter++;
        outer_iter++;
        double t11 = clock();
        auto [Q1, T1] = SolveApproximatly(W, Q, T, C, beta, pi, WFnorm);
        Q = Q1;
        T = T1;
        curr_err = (W - Q * T).squaredNorm();
        if (UpdateRule == 1)
        {
            if (outer_iter % 8 == 0)
                beta = min(100000.0, beta * factor);
            else
            {
                if (curr_err / last_err > 0.95)
                    beta = beta * 5;
                last_err = curr_err;
            }
        }
        double improve = WFnorm * WFnorm - Q.squaredNorm() * Q.squaredNorm();
        if (curr_err < err_eps)
        {
            if (IncreaseRank == 0)
                break;
            if (improve > 0)
                break;
        }
        pi = pi + beta * (W - Q * T);
        double t22 = clock();
        printf("iter: %d, |W-QT|_F: %.5f, beta:%.2f, improve: %.2e, cur rank:%ld, time(s): %f\n", outer_iter, curr_err, beta, improve, T.rows(), t22 - t11);
        if (inner_iter == max_iter || beta > M)
        {
            if (IncreaseRank == 0)
                break;
            inner_iter = 0;
            int cur_r = T.rows();
            int set_rank = cur_r + increasement;
            set_rank = min(set_rank, int(W.cols()));
            // Q(:,(cur_r+1):set_rank)=0.001;
            // T((cur_r+1):set_rank,:)=0.001;
            beta = 10;
            pi = MatrixXd::Zero(W.rows(), W.cols());
            last_err = inf;
            continue;
        }
    }
    double t2 = clock();
    double improve = WFnorm * WFnorm - Q.squaredNorm() * Q.squaredNorm();
    if (improve < 0 || curr_err > 0.1)
    {
        flag = -1;
        printf("LowRankDP Fail!\n");
    }

    // fprintf("gamma: %f, Sen: %f, W: %f, Q: %f, IError:%f, ExpectedError:%f, totaltime:%f\n", ... norm(W - Q * T, 'fro'), L1sensitivity(T), WFnorm ^ 2, norm(Q, 'fro') ^ 2, ComputeExpectedError(W, 1), ComputeExpectedError(Q, C), ts);
}

int main()
{
    // MatrixXd m = MatrixXd::Random(4, 4);

    // cout << m << '\n';
    // m = m.array() - 1.5;
    // cout << m << '\n';
    // auto v = m.col(1);
    // m = m.cwiseMax(0.0);
    // cout << m << endl;
    // // sort(v.begin(), v.end(), greater<double>());
    // // cout << v << endl;
    int N = 10;
    int M = 128;
    MatrixXd W = MatrixXd::Random(N, M).unaryExpr([&](double x) -> double
                                                  { return sign(x); });
    MatrixXd domain_x = MatrixXd::Random(M, 1).unaryExpr([](double x) -> double
                                                         {if(x<0)return 0.; return x; });
    domain_x = 10 * domain_x;
    auto [Q, T] = LowRankDP(W);

    int times = 20;
    double epsi = 0.1;

    double errI = Compute_I_Eror(W, domain_x, epsi, times);
    printf("Laplacian Mechanism (noise on data): %e\n", errI);

    double errO = Compute_QT_Eror(W, Q, T, domain_x, epsi, times);
    printf("Low-Rank Mechanism: %e\n", errO);
}
