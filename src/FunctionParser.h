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

#ifndef FPARSE
#define FPARSE
using namespace std;
#include <string>
#include "exprtk.hpp"
class FunctionParser
{
public:
    double d_expr;
    double n_expr;
    double *d_expr_addr;
    double *n_expr_addr;
    string expression_string;
    exprtk::expression<double> expression;
    exprtk::symbol_table<double> symbol_table;
    bool is_compiled = false;
    void parse(string expression_string);
    double evaluate(double n, double d);
    void destroy();
    static int evaluateFunction(string expression_string, int n, int d);
};
#endif
