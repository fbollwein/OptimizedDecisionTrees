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
using namespace std;
class Costs
{
public:
    vector<double> costsX_u;
    vector<double> unitsX_u;
    vector<double> deltaX_u;
    vector<double> costsN_u;
    vector<double> unitsN_u;
    vector<double> deltaN_u;
    vector<double> costsX_d;
    vector<double> unitsX_d;
    vector<double> deltaX_d;
    vector<double> costsN_d;
    vector<double> unitsN_d;
    vector<double> deltaN_d;
    vector<double> costs;
    vector<double> units;
    vector<double> delta;
    double avgX_u;
    double avgN_u;
    double avgX_d;
    double avgN_d;
    double *X;
    int *Y;
    int total_points;
    int n_classes;
    Costs(){};
    Costs(double *X, int *Y, int total_points, int n_classes);
    void update(int ind, double delta, double f, bool x, bool up);
};