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

#include <stdlib.h>
#include <vector>
#include "Vector_reg.h"
#ifndef NNODEREG
#define NNODEREG
class NNode_reg
{
public:
    vector<Vector_reg> A;
    int depth;
    int N;
    int L;
    int criterion;
    double LB;
    double UB;
    NNode_reg();
    NNode_reg(vector<Vector_reg> A, vector<vector<int>> fixed, vector<int> inner, vector<Vector_reg> fixed_y, vector<double> maxL, vector<double> minL, int N, int L, int depth, int criterion);
    vector<vector<int>> fixed;
    vector<Vector_reg> fixed_y;
    vector<int> inner;
    vector<vector<int>> fixed_best;
    vector<double> maxL;
    vector<double> minL;
    vector<double> maxi;
    vector<double> mini;
    vector<double> meani;
    void upperBound();
    void lowerBound();
    bool isLinearSeparable();
    vector<NNode_reg> branch();
    int partition_size();
};
#endif