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

#include "FunctionParser.h"
#include <math.h>
#include <iostream>
using namespace std;
void FunctionParser::parse(string expression_string)
{
    this->expression_string = expression_string;
    d_expr = 5;
    n_expr = 5;
    d_expr_addr = &d_expr;
    n_expr_addr = &n_expr;
    symbol_table.add_variable("n", n_expr);
    symbol_table.add_variable("d", d_expr);
    expression.register_symbol_table(symbol_table);
    exprtk::parser<double> parser;
    if (!parser.compile(expression_string, expression))
    {
        cout<<"Expression compilation error...\n";
        is_compiled = false;
    }
    else
    {
        is_compiled = true;
    }
}
double FunctionParser::evaluate(double d, double n)
{
    *d_expr_addr = d;
    *n_expr_addr = n;
    return round(expression.value());
}
int FunctionParser::evaluateFunction(string expression_string, int n, int d)
{
    double d_expr = d;
    double n_expr = n;
    exprtk::symbol_table<double> symbol_table;
    symbol_table.add_variable("n", n_expr);
    symbol_table.add_variable("d", d_expr);
    exprtk::expression<double> expression;
    expression.register_symbol_table(symbol_table);
    exprtk::parser<double> parser;
    if (!parser.compile(expression_string, expression))
    {
        cout<<"Expression compilation error...\n";
        return -1;
    }
    else
    {
        return round(expression.value());
    }
}
void FunctionParser::destroy(){
    is_compiled = false;
}