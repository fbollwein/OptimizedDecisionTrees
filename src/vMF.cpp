/*Copyright (C) 2022  Ferdinand Bollwein

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.*/

#include "vMF.h"
#include <chrono>
#include "eigenmvn.h"
#include "boost/math/special_functions/bessel.hpp"
#include <thread>
using namespace Eigen;
using namespace std;
std::mutex mtx3;
vMF::vMF(Environment *env, int m)
{
    this->env = env;
    uniform = std::uniform_real_distribution<double>(0.0, 1.0);
    beta = boost::random::beta_distribution<>((double)(m - 1) / 2., (double)(m - 1) / 2.);
    normal = std::normal_distribution<double>(0.0, 1.0);
}
VectorXd vMF::rw(int n, double lambda, int m)
{
    VectorXd W(n);
    double d = m - 1;
    double b = d / (2. * lambda + sqrt(4. * lambda * lambda + d * d));
    double x = (1. - b) / (1. + b);
    double c = lambda * x + d * log(1. - x * x);
    double w, Z, U;
    for (int i = 0; i < n; ++i)
    {
        bool cont = true;
        while (cont)
        {
            cont = false;
            mtx3.lock();
            Z = beta(env->generator);
            mtx3.unlock();
            w = (1. - (1. + b) * Z) / (1. - (1. - b) * Z);
            mtx3.lock();
            U = uniform(env->generator);
            mtx3.unlock();
            while (U == 0. || U == 1.)
            {
                mtx3.lock();
                U = uniform(env->generator);
                mtx3.unlock();
            }
            if (lambda * w + d * log(1. - x * w) - c < log(1 - U))
            {
                cont = true;
            }
        }
        W(i) = w;
    }
    return W;
}
MatrixXd vMF::rvMF(int size, VectorXd &theta)
{
    double eps = 1e-10;
    int m = theta.size();
    double lambda = theta.lpNorm<2>();
    MatrixXd X(size, m);
    if (lambda < eps)
    {
        X = MatrixXd::Zero(size, m);
        bool initialized = false;
        while (X.isZero(eps))
        {
            for (int i = 0; i < size; i++)
            {
                for (int j = 0; j < m; j++)
                {
                    mtx3.lock();
                    X(i, j) = normal(env->generator);
                    mtx3.unlock();
                }
            }
            initialized = true;
        }
    }
    else
    {
        int d = m - 1;
        VectorXd W = rw(size, lambda, m);
        VectorXd mu = theta / lambda;
        MatrixXd V(size, d);
        bool initialized = false;
        while (V.isZero(eps) || !initialized)
        {
            for (int i = 0; i < size; i++)
            {
                for (int j = 0; j < d; j++)
                {
                    mtx3.lock();
                    V(i, j) = normal(env->generator);
                    mtx3.unlock();
                }
            }
            initialized = true;
        }
        V.rowwise().normalize();
        for (int i = 0; i < size; i++)
        {
            V.row(i) = sqrt(1.0 - (W(i) * W(i))) * V.row(i);
        }
        X << W, V;
        VectorXd basis(m);
        for (int i = 0; i < m; i++)
        {
            basis(i) = 0;
        }
        basis(0) = 1;
        VectorXd u = basis - mu;
        u.normalize();
        X = X - 2 * X * u * u.transpose();
    }
    X.rowwise().normalize();
    return X;
}
Eigen::MatrixXd rvMF_t(int id, vMF *vMF, int size, VectorXd theta)
{
    return vMF->rvMF(size, theta);
}
Eigen::MatrixXd vMF::rvMF_t2(int id, int size, VectorXd theta)
{
    return rvMF(size, theta);
}
Eigen::MatrixXd vMF::rvMF_parallel(int size, VectorXd theta, ctpl::thread_pool *p)
{
    int n_threads = p->size();
    int chunk_size = size / n_threads;
    vector<future<MatrixXd>> fs;
    int sum = 0;
    for (int k = 0; k < n_threads; k++)
    {
        if (k != n_threads - 1 && chunk_size != 0)
        {
            fs.push_back(p->push([&](int id)
                                 { return this->rvMF(chunk_size, theta); }));
            sum += chunk_size;
        }
        else
        {
            if (sum != size)
            {
                fs.push_back(p->push([&](int id)
                                     { return this->rvMF(size - sum, theta); }));
            }
        }
    }
    for (int k = 0; k < fs.size(); k++)
    {
        fs[k].wait();
    }
    int n = 0;
    MatrixXd X(size, theta.size());
    for (int k = 0; k < n_threads; k++)
    {
        MatrixXd tmp = fs[k].get();
        for (int i = 0; i < tmp.rows(); i++)
        {
            X.row(n) = tmp.row(i);
            n++;
        }
    }
    return X;
}
MatrixXd vMF::rvMF2(int size, VectorXd theta)
{
    int m = theta.size();
    MatrixXd X(size, m);
    for (int n = 0; n < size; n++)
    {
        Eigen::MatrixXd sampleM = rvMF(1, theta);
        for (int i = 0; i < m; i++)
        {
            X(n, i) = sampleM(0, i);
        }
    }
    return X;
}
double approxRatioOfBessel(double d, double kappa)
{
    double a1_n = d / kappa;
    double a2_n = ((d + 2) / kappa);
    double num = a1_n * (a2_n + 1. / a1_n);
    a1_n = (a2_n + 1. / a1_n);
    double denom = a2_n;
    int n = 3;
    double an = 1;
    double A = 0;
    while (n <= 100000)
    {
        an = (d + 2 * (n - 1)) / kappa;
        if (n % 2 == 0)
        {
            an = 1 * an;
        }
        a1_n = an + 1. / a1_n;
        a2_n = an + 1. / a2_n;
        if (abs(a1_n - a2_n) < 1e-14)
        {
            break;
        }
        num = num * a1_n;
        denom = denom * a2_n;
        n++;
    }
    A = num / denom;
    return 1. / A;
}
double vMF::approxKappa(double kappa, int d, double R, int iterations)
{
    for (int i = 0; i < iterations; i++)
    {
        double A = approxRatioOfBessel(d, kappa);
        kappa = kappa - (A - R) / (1 - A * A - ((d - 1) / kappa) * A);
    }
    return kappa;
}