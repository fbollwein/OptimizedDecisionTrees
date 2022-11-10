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
#include "Eta.h"
#include <iostream>
using namespace std;
Eta::Eta(Eigen::VectorXd w, int col)
{
    this->w = w;
    this->col = col;
    invert = false;
    for (int i = 0; i < w.size(); i++)
    {
        if (w(i) != 0)
        {
            indices.push_back(i);
        }
    }
}
Eta::Eta(int col)
{
    this->col = col;
    invert = true;
}
void Eta::solve(Eigen::VectorXd &a)
{
    if (invert)
    {
        a(col) = -a(col);
    }
    else
    {
        a(col) = a(col) / w(col);
        for (int k = 0; k < indices.size(); k++)
        {
            if (indices[k] == col)
            {
                continue;
            }
            a(indices[k]) = a(indices[k]) - w(indices[k]) * a(col);
        }
    }
}
void Eta::solve_transpose(Eigen::VectorXd &c)
{
    if (invert)
    {
        c(col) = -c(col);
    }
    else
    {
        double sum = 0;
        int j;
        for (int k = 0; k < indices.size(); k++)
        {
            j = indices[k];
            if (j != col)
            {
                sum += c(j) * w(j);
            }
        }
        c(col) = (c(col) - sum) / w(col);
    }
}