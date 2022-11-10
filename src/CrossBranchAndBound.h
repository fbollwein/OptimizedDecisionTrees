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
#include "PointValuePair.h"
#include "Result.h"
#include "Environment.h"
using namespace std;
void CrossBranchAndBound(int id, double X1[], double X2[], int ind1, int ind2, int n_classes, int total_points,
                         int **dataX2, int **dataY2, double **dataX1_reg, double **dataX2_reg, int cX[], int cY[], int countX, int countY, Result *result, Environment *env);
