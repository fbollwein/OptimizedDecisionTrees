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

#include "NNode_reg.h"
#include <limits>
#include <algorithm>
#include "Environment.h"
NNode_reg::NNode_reg(){
};
NNode_reg::NNode_reg(vector<Vector_reg> A, vector<vector<int>> fixed, vector<int> inner, vector<Vector_reg> fixed_y, vector<double> maxL, vector<double> minL,
                     int N, int L, int depth, int criterion)
{
    this->A = A;
    this->fixed = fixed;
    this->inner = inner;
    this->fixed_y = fixed_y;
    this->N = N;
    this->L = L;
    this->depth = depth;
    this->criterion = criterion;
    this->maxL = maxL;
    this->minL = minL;
    this->maxi = maxi;
    this->mini = mini;
    this->meani = meani;
};
vector<NNode_reg> NNode_reg::branch()
{
    vector<NNode_reg> nodes;
    if (inner.size() == 0)
    {
        return nodes;
    }
    int i = inner[inner.size() - 1];
    double mean = 0;
    int fixed_l = -1;
    if (criterion == Environment::mse)
    {
        if (A[i].size() > 0)
        {
            mean = A[i].sum_t / A[i].size();
            for (int l = 0; l < L; l++)
            {
                if (mean < maxL[l] && mean > minL[l])
                {
                    fixed_l = l;
                }
            }
        }
    }
    for (int l = 0; l < L; l++)
    {
        if (fixed_l != -1)
        {
            l = fixed_l;
        }
        else
        {
        }
        bool flag = false;
        if (fixed[l].size() == 0)
        {
            flag = true;
        }
        vector<int> inner2;
        for (int i2 = 0; i2 < inner.size(); i2++)
        {
            if (inner[i2] != i)
            {
                inner2.push_back(inner[i2]);
            }
        }
        vector<vector<int>> fixed2;
        vector<Vector_reg> fixed_y2;
        for (int l2 = 0; l2 < L; l2++)
        {
            fixed_y2.push_back(Vector_reg(fixed_y[l2]));
            vector<int> f;
            fixed2.push_back(f);
            for (int i2 = 0; i2 < fixed[l2].size(); i2++)
            {
                fixed2[l2].push_back(fixed[l2][i2]);
            }
        }
        fixed2[l].push_back(i);
        fixed_y2[l].merge(A[i]);
        vector<double> maxL2 = maxL;
        vector<double> minL2 = minL;
        if (criterion == Environment::mse && A[i].size() > 0)
        {
            if (mean > maxL2[l])
            {
                maxL2[l] = mean;
            }
            if (mean < minL2[l])
            {
                minL2[l] = mean;
            }
        }
        NNode_reg node(A, fixed2, inner2, fixed_y2, maxL2, minL2, N, L, depth + 1, criterion);
        nodes.push_back(node);
        if (fixed_l != -1)
        {
            break;
        }
        if (flag)
        {
            break;
        }
    }
    if (nodes.size() > 0)
    {
    }
    return nodes;
}
void NNode_reg::lowerBound()
{
    LB = 0;
    for (int l = 0; l < fixed_y.size(); l++)
    {
        LB += fixed_y[l].error(criterion);
    }
    for (int i2 = 0; i2 < inner.size(); i2++)
    {
        LB += A[inner[i2]].error(criterion);
    }
}
void NNode_reg::upperBound()
{
    UB = std::numeric_limits<double>::max();
    if (inner.size() == 0)
    {
        UB = 0;
        for (int l = 0; l < fixed_y.size(); l++)
        {
            UB += fixed_y[l].error(criterion);
        }
    }
    fixed_best = fixed;
}
int NNode_reg::partition_size()
{
    int size = 0;
    for (int l = 0; l < L; l++)
    {
        if (fixed[l].size() > 0)
        {
            size += 1;
        }
    }
    return size;
}