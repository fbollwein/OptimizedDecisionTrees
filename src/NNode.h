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
#include "General.h"
#include "Environment.h"
#ifndef NNODE
#define NNODE
struct partitionpair
{
  vector<vector<int>> part;
  double impurity;
};
class NNode
{
public:
  int depth;
  double LB;
  double UB;
  int N;
  int K;
  int L;
  vector<vector<int>> fixed;
  int **A;
  vector<int> inner;
  vector<double> hyp;
  vector<vector<double>> distr;
  vector<double> distr_i;
  int criterion;
  Environment *env;
  vector<double> maxL;
  vector<double> minL;
  bool is_binary;
  int sortclass;
  NNode();
  NNode(int **A, vector<vector<int>> fixed, vector<int> &inner,
        vector<vector<double>> &distr, vector<double> &distr_i, int N, int L, int K,
        vector<double> hyp, int criterion, int depth,
        vector<double> minL, vector<double> maxL, bool is_binary, int sortclass, Environment *env);
  void upperBound();
  void lowerBound();
  #ifdef use_gurobi
  bool isLinearSeparable();
  #endif
  vector<NNode> branch();
  void lowerbound(vector<double> impurities);
  void upperbound(vector<double> impurities);
  int partition_size();
  pair<double, vector<vector<int>>> heuristic();
};
#endif