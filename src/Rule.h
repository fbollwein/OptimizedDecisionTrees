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

#include <string>
#include <vector>
#include <iostream>
#include <limits>
using namespace std;
#ifndef RULE
#define RULE
class Rule
{
public:
        bool split_found = false;
        double imp = std::numeric_limits<double>::max();
        bool is_nominal = false;
        bool is_cross = false;
        int no_children = 0;
        int no_nominal_partitions;
        vector<double> a;
        double b;
        int nominal_ind;
        vector<vector<int>> partition;
        vector<int> which_partition;
        vector<int> attribute_order;
        vector<bool> single_considered;
        int ind1 = 0;
        int ind2 = 0;
        double x1 = 0;
        double x2 = 0;
        bool divides = false;
        bool leq = false;
        double opt_gap = 0;
        double time = -1;
        Rule();
        Rule(vector<double> a, double b);
        Rule(int nom_feature, vector<int> which_partition);
        Rule(int ind1, int ind2, double x1, double x2);
        int evaluate(int i, double **X, double **Xnominal);
        int evaluateNumeric(int i, double **X);
        int evaluateCross(int i, double **X);
        int evaluateNominal(int i, double **Xnominal);
        void normalize();
};
#endif