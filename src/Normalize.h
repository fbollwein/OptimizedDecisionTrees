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

#include "Result.h"
using namespace std;
void normalize(int ***data, double **X2, double **X, int total_points, int d, double *offset, double *denominator, int type);
void normalize_reg(double ***data, double **X2, double **X, int total_points, int d, double *offset, double *denominator, int type);
void denormalize(Result *res, int d, double *offset, double *denominator);