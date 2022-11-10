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

#include "Result.h"
#include <limits>
Result::Result()
{
	this->a = vector<double>(0, 0);
}
Result::Result(int total_features)
{
	this->id = 0;
	this->nextid = 0;
	this->split_found = false;
	this->a = vector<double>(0, 0);
	this->bestImp = std::numeric_limits<double>::max();
	for (int j = 0; j < total_features; j++)
	{
		a.push_back(1);
	}
	this->bestB = 1;
}
Result::Result(bool split_found, double bestImp, double bestA, double bestB,
			   int bestInd1, int bestInd2, int total_features)
{
	this->id = 0;
	this->nextid = 0;
	this->split_found = split_found;
	this->bestImp = bestImp;
	this->bestA = bestA;
	this->bestB = bestB;
	this->bestInd1 = bestInd1;
	this->bestInd2 = bestInd2;
	this->a = vector<double>(0, 0);
	for (int j = 0; j < total_features; j++)
	{
		if (j == bestInd1)
		{
			if (j == bestInd2)
			{
				a.push_back(1);
			}
			else
			{
				a.push_back(-bestA);
			}
		}
		else if (j == bestInd2)
		{
			a.push_back(1);
		}
		else
		{
			a.push_back(0);
		}
	}
}
Result::Result(bool split_found, double bestImp, vector<double> a, double b)
{
	this->nextid = 0;
	this->id = 0;
	this->split_found = split_found;
	this->bestImp = bestImp;
	for (int j = 0; j < a.size(); j++)
	{
		this->a.push_back(a[j]);
	}
	this->bestB = b;
}
int Result::getId()
{
	return id;
}
int Result::generateId()
{
	nextid += 1;
	return nextid;
}
