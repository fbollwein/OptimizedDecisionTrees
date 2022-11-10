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

#include <vector>
#include "Result.h"
#include "Environment.h"
#include "ctpl_stl.h"
using namespace std;
void BrightSide(ctpl::thread_pool *p, int ***data, double ***data_reg, double **X2, vector<int> order, double **X, int *Y, double *Y_reg, int total_points, int d, int n_classes, double *distr_S_G, Result *res, Environment *env, int norm, double start_time);
void SimpleCrossEntropy(ctpl::thread_pool *p, int ***data, double ***data_reg, double **X2, vector<int> order, double **X, int *Y, double *Y_reg, int total_points, int d, int d_orig, int n_classes, double *distr_S_G, Result *res, Environment *env, int norm, double start_time);
