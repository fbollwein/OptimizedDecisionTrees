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

#include "Costs.h"
#include <iostream>
Costs::Costs(double *X, int *Y, int total_points, int n_classes)
{
    this->X = X;
    this->Y = Y;
    this->total_points = total_points;
    this->n_classes = n_classes;
    costsX_u = vector<double>(total_points, 1);
    unitsX_u = vector<double>(total_points, 0);
    deltaX_u = vector<double>(total_points);
    costsX_d = vector<double>(total_points, 1);
    unitsX_d = vector<double>(total_points, 0);
    deltaX_d = vector<double>(total_points, 0);
    costsN_u = vector<double>(n_classes, 1);
    unitsN_u = vector<double>(n_classes, 0);
    deltaN_u = vector<double>(n_classes, 0);
    costsN_d = vector<double>(n_classes, 1);
    unitsN_d = vector<double>(n_classes, 0);
    deltaN_d = vector<double>(n_classes, 0);
    avgX_u = 1;
    avgX_d = 1;
    avgN_u = 1;
    avgN_d = 1;
}
void Costs::update(int ind, double delta, double f, bool x, bool up)
{
    if (x)
    {
        if (up)
        {
            if (unitsX_u[ind] == 0)
            {
                deltaX_u[ind] = 0;
            }
            deltaX_u[ind] += (delta / f);
            unitsX_u[ind] += 1;
            costsX_u[ind] = deltaX_u[ind] / unitsX_u[ind];
        }
        else
        {
            if (unitsX_d[ind] == 0)
            {
                deltaX_d[ind] = 0;
            }
            deltaX_d[ind] += (delta / f);
            unitsX_d[ind] += 1;
            costsX_d[ind] = deltaX_d[ind] / unitsX_d[ind];
        }
    }
    else
    {
        if (up)
        {
            if (unitsN_u[ind] == 0)
            {
                deltaN_u[ind] = 0;
            }
            deltaN_u[ind] += (delta / f);
            unitsN_u[ind] += 1;
            costsN_u[ind] = deltaN_u[ind] / unitsN_u[ind];
        }
        else
        {
            if (unitsN_d[ind] == 0)
            {
                deltaN_d[ind] = 0;
            }
            deltaN_d[ind] += (delta / f);
            unitsN_d[ind] += 1;
            costsN_d[ind] = deltaN_d[ind] / unitsN_d[ind];
        }
    }
    if (up)
    {
        if (x)
        {
            avgX_u = 0;
            int n = 0;
            for (int i = 0; i < total_points; i++)
            {
                if (unitsX_u[i] > 0)
                {
                    avgX_u += costsX_u[i];
                    n++;
                }
            }
            avgX_u = avgX_u / n;
        }
        else
        {
            avgN_u = 0;
            int n = 0;
            for (int i = 0; i < n_classes; i++)
            {
                if (unitsN_u[i] > 0)
                {
                    avgN_u += costsN_u[i];
                    n++;
                }
            }
            avgN_u = avgN_u / n;
        }
    }
    else
    {
        if (x)
        {
            avgX_d = 0;
            int n = 0;
            for (int i = 0; i < total_points; i++)
            {
                if (unitsX_d[i] > 0)
                {
                    avgX_d += costsX_d[i];
                    n++;
                }
            }
            avgX_d = avgX_d / n;
        }
        else
        {
            avgN_u = 0;
            int n = 0;
            for (int i = 0; i < n_classes; i++)
            {
                if (unitsN_u[i] > 0)
                {
                    avgN_d += costsN_d[i];
                    n++;
                }
            }
            avgN_d = avgN_d / n;
        }
    }
}
