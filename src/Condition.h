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
#include <map>
using namespace std;
class Condition
{
public:
  vector<double> a;
  double b;
  int ind1;
  int ind2;
  double x1;
  double x2;
  string condition;
  static const int leq = -1;
  static const int l = -2;
  static const int eq = 0;
  static const int geq = 1;
  static const int g = 2;
  static const int t = 3;
  int op;
  int op2;
  Condition();
  Condition(vector<double> a, double b, int op);
  Condition(int ind1, int ind2, double x1, double x2, int op, int op2);
  Condition(int nom_feature, vector<bool> in_partition, vector<vector<bool>> is_nominal_val_in_dataset);
  Condition copy();
  bool evaluate(vector<double> &x, vector<double> &x_nominal);
  string toString();
  void invert();
  vector<string> names;
  vector<string> nominalNames;
  map<int, map<int, string>> cat_to_name;
  bool isNominal = false;
  bool isCross = false;
  int nom_feature;
  vector<vector<bool>> is_nominal_val_in_dataset;
  vector<bool> in_partition;
};
