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

#include "Rule.h"
#include "Environment.h"
#include "NNode.h"
#include "NNode_reg.h"
#include "ctpl_stl.h"
using namespace std;
void NBranchAndBound_class(double **Xnominal, int *Y, int ***data, int total_points, int n_classes, int N, int which_feature, vector<vector<bool>> is_nominal_val_in_dataset, Rule *rule, Environment *env, ctpl::thread_pool *p);
void NBranchAndBound_reg(double **Xnominal, double *Y_reg, double ***data_reg, int total_points, int N, int which_feature, vector<vector<bool>> is_nominal_val_in_dataset, Rule *rule, Environment *env);