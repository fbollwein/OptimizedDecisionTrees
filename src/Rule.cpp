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
#include <math.h>
Rule::Rule(){};
Rule::Rule(vector<double> a, double b)
{
    is_nominal = false;
    is_cross = false;
    this->a = a;
    this->b = b;
    no_children = 2;
};
Rule::Rule(int ind1, int ind2, double x1, double x2)
{
    is_nominal = false;
    is_cross = true;
    this->ind1 = ind1;
    this->ind2 = ind2;
    this->x1 = x1;
    this->x2 = x2;
    no_children = 4;
};
Rule::Rule(int nom_feature, vector<int> which_partition)
{
    is_nominal = true;
    is_cross = false;
    nominal_ind = nom_feature;
    no_nominal_partitions = 0;
    for (int i = 0; i < which_partition.size(); i++)
    {
        if (which_partition[i] > no_nominal_partitions)
        {
            no_nominal_partitions = which_partition[i];
        }
    }
    no_nominal_partitions += 1;
    no_children = no_nominal_partitions;
    this->which_partition = which_partition;
};
void Rule::normalize()
{
    if (is_nominal || is_cross)
    {
        return;
    }
    else
    {
        double norm = b * b;
        int used_features = 0;
        for (int j = 0; j < a.size(); j++)
        {
            norm += a[j] * a[j];
            if (abs(a[j]) != 0)
            {
                used_features++;
            }
        }
        if (used_features > 2)
        {
            norm = sqrt(norm);
            for (int j = 0; j < a.size(); j++)
            {
                a[j] = a[j] / norm;
            }
            b = b / norm;
        }
    }
}
int Rule::evaluate(int i, double **X, double **Xnominal)
{
    if (is_nominal)
    {
        return evaluateNominal(i, Xnominal);
    }
    if (is_cross)
    {
        return evaluateCross(i, X);
    }
    return evaluateNumeric(i, X);
};
int Rule::evaluateCross(int i, double **X)
{
    if (X[ind1][i] <= x1)
    {
        if (X[ind2][i] <= x2)
        {
            return 0;
        }
        else
        {
            return 2;
        }
    }
    else
    {
        if (X[ind2][i] <= x2)
        {
            return 1;
        }
        else
        {
            return 3;
        }
    }
};
int Rule::evaluateNumeric(int i, double **X)
{
    if (!leq)
    {
        double sum = -b;
        for (int j = 0; j < a.size(); j++)
        {
            sum += a[j] * X[j][i];
        }
        if (sum >= 0)
        {
            return 0;
        }
        else
        {
            return 1;
        }
    }
    else
    {
        double sum = -b;
        for (int j = 0; j < a.size(); j++)
        {
            sum += a[j] * X[j][i];
        }
        if (sum > 0)
        {
            return 0;
        }
        else
        {
            return 1;
        }
    }
};
int Rule::evaluateNominal(int i, double **Xnominal)
{
    return which_partition[Xnominal[i][nominal_ind]];
}