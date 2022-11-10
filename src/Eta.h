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

#define EIGEN_NO_DEBUG
#include <string>
#include "Eigen/Sparse"
using namespace std;
#ifndef ETA_H
#define ETA_H
class Eta
{
public:
        Eigen::VectorXd w;
        bool invert;
        int col;
        vector<int> indices;
        double tol = 1e-30;
        Eta(Eigen::VectorXd w, int col);
        Eta(int col);
        void solve(Eigen::VectorXd &a);
        void solve_transpose(Eigen::VectorXd &c);
};
#endif