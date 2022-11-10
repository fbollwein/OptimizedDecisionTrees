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

#include "General.h"
#include <stdlib.h>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include "Environment.h"
using namespace std;
double gini(double *N1, int n_classes)
{
	double gini = 0;
	double No = 0;
	double No2 = 0;
	double imp1 = 0;
	double imp2 = 0;
	for (int m = 0; m < n_classes; m++)
	{
		No += N1[m];
		No2 += N1[n_classes + m];
	}
	if (No != 0)
	{
		for (int m = 0; m < n_classes; m++)
		{
			if (N1[m] > 0)
			{
				imp1 += (N1[m] / No) * N1[m];
			}
		}
		gini += No - (imp1);
	}
	if (No2 != 0)
	{
		for (int m = 0; m < n_classes; m++)
		{
			if (N1[n_classes + m] > 0)
			{
				imp2 += (N1[n_classes + m] / No2) * N1[n_classes + m];
			}
		}
		gini += No2 - (imp2);
	}
	return gini;
}
double gini(vector<double> &N1, int n_classes)
{
	double gini = 0;
	double No = 0;
	double No2 = 0;
	double imp1 = 0;
	double imp2 = 0;
	if (N1.size() == n_classes)
	{
		for (int m = 0; m < n_classes; m++)
		{
			No += N1[m];
		}
	}
	else
	{
		for (int m = 0; m < n_classes; m++)
		{
			No += N1[m];
			No2 += N1[n_classes + m];
		}
	}
	if (No != 0)
	{
		for (int m = 0; m < n_classes; m++)
		{
			if (N1[m] > 0)
			{
				imp1 += (N1[m] / No) * N1[m];
			}
		}
		gini += No - (imp1);
	}
	if (No2 != 0)
	{
		for (int m = 0; m < n_classes; m++)
		{
			if (N1[n_classes + m] > 0)
			{
				imp2 += (N1[n_classes + m] / No2) * N1[n_classes + m];
			}
		}
		gini += No2 - (imp2);
	}
	return gini;
}
double gini(int *N1, int *N2, int n_classes)
{
	double gini = 0;
	double No = 0;
	double No2 = 0;
	double imp1 = 0;
	double imp2 = 0;
	for (int m = 0; m < n_classes; m++)
	{
		No += N1[m];
		No2 += N2[m];
	}
	if (No != 0)
	{
		for (int m = 0; m < n_classes; m++)
		{
			if (N1[m] > 0)
			{
				imp1 += (N1[m] / No) * N1[m];
			}
		}
		gini += No - (imp1);
	}
	if (No2 != 0)
	{
		for (int m = 0; m < n_classes; m++)
		{
			if (N2[m] > 0)
			{
				imp2 += (N2[m] / No2) * N2[m];
			}
		}
		gini += No2 - (imp2);
	}
	return gini;
}
double gini1(double *N1, int n_classes)
{
	double gini = 0;
	double No = 0;
	double imp1 = 0;
	double qu = 0;
	for (int m = 0; m < n_classes; m++)
	{
		No += N1[m];
	}
	if (No != 0)
	{
		for (int m = 0; m < n_classes; m++)
		{
			if (N1[m] > 0)
			{
				imp1 += (N1[m] / No) * N1[m];
			}
		}
		gini += No - (imp1);
	}
	return gini;
}
double lowerboundgini(double *N1, int n_classes)
{
	int max = 0;
	int max_ind = 0;
	for (int k = 0; k < n_classes; k++)
	{
		if (N1[2 * n_classes + k] > max)
		{
			max_ind = k;
			max = N1[2 * n_classes + k];
		}
	}
	N1[2 * n_classes + max_ind] = 0;
	double tmp = gini1(N1, n_classes) + gini1(&N1[n_classes], n_classes);
	N1[2 * n_classes + max_ind] = max;
	return tmp;
}
double entropy(double *N1, int n_classes)
{
	double entropy = 0;
	double No = 0;
	double No2 = 0;
	for (int m = 0; m < n_classes; m++)
	{
		No += N1[m];
		No2 += N1[n_classes + m];
	}
	if (No != 0)
	{
		for (int m = 0; m < n_classes; m++)
		{
			if (N1[m] > 0)
			{
				entropy -= N1[m] * log2((N1[m] / No));
			}
		}
	}
	if (No2 != 0)
	{
		for (int m = 0; m < n_classes; m++)
		{
			if (N1[n_classes + m] > 0)
			{
				entropy -= N1[n_classes + m] * log2(N1[n_classes + m] / No2);
			}
		}
	}
	return entropy;
}
double entropy(vector<double> &N1, int n_classes)
{
	double entropy = 0;
	double No = 0;
	double No2 = 0;
	if (N1.size() == n_classes)
	{
		for (int m = 0; m < n_classes; m++)
		{
			No += N1[m];
		}
	}
	else
	{
		for (int m = 0; m < n_classes; m++)
		{
			No += N1[m];
			No2 += N1[n_classes + m];
		}
	}
	if (No != 0)
	{
		for (int m = 0; m < n_classes; m++)
		{
			if (N1[m] > 0)
			{
				entropy -= N1[m] * log2((N1[m] / No));
			}
		}
	}
	if (No2 != 0)
	{
		for (int m = 0; m < n_classes; m++)
		{
			if (N1[n_classes + m] > 0)
			{
				entropy -= N1[n_classes + m] * log2(N1[n_classes + m] / No2);
			}
		}
	}
	return entropy;
}
double entropy(int *N1, int *N2, int n_classes)
{
	double entropy = 0;
	double No = 0;
	double No2 = 0;
	for (int m = 0; m < n_classes; m++)
	{
		No += N1[m];
		No2 += N2[m];
	}
	if (No != 0)
	{
		for (int m = 0; m < n_classes; m++)
		{
			if (N1[m] > 0)
			{
				entropy -= N1[m] * log2((N1[m] / No));
			}
		}
	}
	if (No2 != 0)
	{
		for (int m = 0; m < n_classes; m++)
		{
			if (N2[m] > 0)
			{
				entropy -= N2[m] * log2(N2[m] / No2);
			}
		}
	}
	return entropy;
}
double entropy1(double *N1, int n_classes)
{
	double entropy = 0;
	double No = 0;
	for (int m = 0; m < n_classes; m++)
	{
		No += N1[m];
	}
	if (No != 0)
	{
		for (int m = 0; m < n_classes; m++)
		{
			if (N1[m] > 0)
			{
				entropy -= (N1[m]) * log2((N1[m] / No));
			}
		}
	}
	return entropy;
}
double imp(double *N1, int n_classes, int type)
{
	if (type == Environment::entropy_imp)
	{
		return entropy(N1, n_classes);
	}
	else if (type == Environment::twoing)
	{
		return twoing(N1, n_classes);
	}
	else
	{
		return gini(N1, n_classes);
	}
}
double imp(vector<double> &N1, int n_classes, int type)
{
	if (type == Environment::entropy_imp)
	{
		return entropy(N1, n_classes);
	}
	else if (type == Environment::twoing)
	{
		return twoing(N1, n_classes);
	}
	else
	{
		return gini(N1, n_classes);
	}
}
double imp(int *N1, int *N2, int n_classes, int type)
{
	if (type == Environment::entropy_imp)
	{
		return entropy(N1, N2, n_classes);
	}
	else if (type == Environment::twoing)
	{
		return twoing(N1, N2, n_classes);
	}
	else
	{
		return gini(N1, N2, n_classes);
	}
}
double imp1(double *N1, int n_classes, int type)
{
	if (type == Environment::entropy_imp)
	{
		return entropy1(N1, n_classes);
	}
	else if (type == Environment::twoing)
	{
		return twoing1(N1, n_classes);
	}
	else
	{
		return gini1(N1, n_classes);
	}
}
double twoing(double *N1, int n_classes)
{
	double imp = 0;
	double No = 0;
	double No2 = 0;
	for (int m = 0; m < n_classes; m++)
	{
		No += N1[m];
		No2 += N1[n_classes + m];
	}
	if (No == 0 || No2 == 0)
	{
		return n_classes * n_classes;
	}
	double N = No + No2;
	for (int m = 0; m < n_classes; m++)
	{
		imp += abs(N1[m] / No - N1[n_classes + m] / No2);
	}
	imp = imp * imp;
	imp = (((No / N) * (No2 / N)) / 4.) * imp;
	imp = (double)(n_classes * n_classes) - imp;
	return imp;
}
double twoing(vector<double> &N1, int n_classes)
{
	double imp = 0;
	double No = 0;
	double No2 = 0;
	if (N1.size() == n_classes)
	{
		for (int m = 0; m < n_classes; m++)
		{
			No += N1[m];
		}
	}
	else
	{
		for (int m = 0; m < n_classes; m++)
		{
			No += N1[m];
			No2 += N1[n_classes + m];
		}
	}
	if (No == 0 || No2 == 0)
	{
		return n_classes * n_classes;
	}
	double N = No + No2;
	double tmp;
	for (int m = 0; m < n_classes; m++)
	{
		imp += abs((N1[m] / No) - (N1[n_classes + m] / No2));
	}
	imp = imp * imp;
	imp = (((No / N) * (No2 / N)) / 4.) * imp;
	imp = ((double)(n_classes * n_classes)) - imp;
	return imp;
}
double twoing(int *N1, int *N2, int n_classes)
{
	double imp = 0;
	double No = 0;
	double No2 = 0;
	for (int m = 0; m < n_classes; m++)
	{
		No += N1[m];
		No2 += N2[m];
	}
	if (No == 0 || No2 == 0)
	{
		return n_classes * n_classes;
	}
	double N = No + No2;
	for (int m = 0; m < n_classes; m++)
	{
		imp += abs(N1[m] / No - N2[m] / No2);
	}
	imp = imp * imp;
	imp = (((No / N) * (No2 / N)) / 4.) * imp;
	imp = (double)(n_classes * n_classes) - imp;
	return imp;
}
double twoing1(double *N1, int n_classes)
{
	return n_classes * n_classes;
}
string itos(int i)
{
	stringstream s;
	s << i;
	return s.str();
}
string dtos(double d)
{
	stringstream s;
	s << d;
	return s.str();
}
