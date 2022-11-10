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

#include "Eigen/Dense"
#include "boost/random.hpp"
#include <random>
#include "ctpl_stl.h"
#include "Environment.h"
using namespace Eigen;
class vMF
{
public:
    Environment *env;
    std::uniform_real_distribution<double> uniform;
    boost::random::beta_distribution<> beta;
    std::normal_distribution<double> normal;
    vMF(Environment *env, int m);
    Eigen::VectorXd rw(int n, double lambda, int m);
    Eigen::MatrixXd rvMF(int size, VectorXd &theta);
    static Eigen::MatrixXd rvMF_t(int id, vMF *vMF, int size, VectorXd theta);
    Eigen::MatrixXd rvMF_t2(int id, int size, VectorXd theta);
    Eigen::MatrixXd rvMF_parallel(int size, VectorXd theta, ctpl::thread_pool *p);
    Eigen::MatrixXd rvMF2(int size, VectorXd theta);
    static double approxKappa(double kappa, int d, double R, int iterations);
};