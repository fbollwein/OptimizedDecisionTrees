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

#ifndef RES
#define RES
#include <vector>
using namespace std;
class Result
{
public:
	int id;
	int nextid;
	bool split_found;
	double bestImp;
	double bestA;
	double bestB;
	double bestInd1;
	double bestInd2;
	vector<double> a;
	bool divides = false;
	bool isCross = false;
	double bestX1;
	double bestX2;
	Result();
	Result(int total_features);
	Result(bool split_found, double bestImp, double bestA, double bestB,
		   int bestInd1, int bestInd2, int total_features);
	Result(bool split_found, double bestImp, vector<double> a, double b);
	int getId();
	int generateId();
};
#endif
