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
#include <iostream>
using namespace std;
class Vector_reg
{
public:
    Vector_reg();
    Vector_reg(vector<double> a);
    Vector_reg(const Vector_reg &a);
    int size();
    double sum();
    double sum2();
    void merge(Vector_reg b);
    void push_back(double b);
    double pop_back();
    double squaredError();
    double absoluteError();
    double squaredError(Vector_reg b);
    void copy(Vector_reg *reg);
    void clear();
    void reserve(int n);
    void sort();
    void insert_inplace(double b);
    Vector_reg add(Vector_reg &b);
    void add(Vector_reg *b);
    double error(int criterion);
    bool remove(double b);
    void use_mae();
    bool is_sorted = false;
    double sum_t = 0;
    double sum2_t = 0;
    vector<double> vec;
    int length = 0;
    double sq_error = -1;
    double a_error = -1;
    bool median = false;
    bool keep_vec = true;
};