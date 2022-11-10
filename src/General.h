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
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
using namespace std;
typedef struct
{
	double dist;
	int sect;
	int y;
} doubleIntPair;
double gini(double *N1, int n_classes);
double gini(vector<double> &N1, int n_classes);
double gini(int *N1, int *N2, int n_classes);
double gini1(double *N1, int n_classes);
double lowerboundgini(double *N1, int n_classes);
double entropy(double *N1, int n_classes);
double entropy(vector<double> &N1, int n_classes);
double entropy(int *N1, int *N2, int n_classes);
double entropy1(double *N1, int n_classes);
double twoing(double *N1, int n_classes);
double twoing(vector<double> &N1, int n_classes);
double twoing(int *N1, int *N2, int n_classes);
double twoing1(double *N1, int n_classes);
double imp(double *N1, int n_classes, int type);
double imp(vector<double> &N1, int n_classes, int type);
double imp(int *N1, int *N2, int n_classes, int type);
double imp1(double *N1, int n_classes, int type);
string itos(int i);
string dtos(double d);
