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

#include <stdlib.h>
#include <vector>
#include <iostream>
#include "CrossNode.h"
#include "General.h"
#include "Environment.h"
using namespace std;
CrossNode::CrossNode() {}
CrossNode::CrossNode(double X1[], double X2[], double *N, Vector_reg **N_reg, int *P,
					 int n_classes, int total_points, int **dataX1, int **dataX2, double **dataX1_reg, double **dataX2_reg,
					 int *cX1, int *cX2, bool regression, int criterion)
{
	this->N = N;
	this->N_reg = N_reg;
	this->n_classes = n_classes;
	this->total_points = total_points;
	this->cX1 = cX1;
	this->cX2 = cX2;
	this->dataX1 = dataX1;
	this->dataX2 = dataX2;
	this->dataX1_reg = dataX1_reg;
	this->dataX2_reg = dataX2_reg;
	this->criterion = criterion;
	this->regression = regression;
	this->X1 = X1;
	this->X2 = X2;
	for (int i = 0; i < 4; i++)
	{
		this->P[i] = P[i];
	}
	if (P[0] == 0)
	{
		if (P[2] == 0)
		{
			pointer[0] = cX1[P[0]];
			pointer[1] = cX1[P[1]];
			pointer[2] = cX2[P[2]];
			pointer[3] = cX2[P[3]];
		}
		else
		{
			pointer[0] = cX1[P[0]];
			pointer[1] = cX1[P[1]];
			pointer[2] = cX2[P[2]] + 1;
			pointer[3] = cX2[P[3]];
		}
	}
	else
	{
		if (P[2] == 0)
		{
			pointer[0] = cX1[P[0]] + 1;
			pointer[1] = cX1[P[1]];
			pointer[2] = cX2[P[2]];
			pointer[3] = cX2[P[3]];
		}
		else
		{
			pointer[0] = cX1[P[0]] + 1;
			pointer[1] = cX1[P[1]];
			pointer[2] = cX2[P[2]] + 1;
			pointer[3] = cX2[P[3]];
		}
	}
	if (!regression)
	{
		xmin = X1[dataX1[pointer[0]][0]];
		xmax = X1[dataX1[pointer[1]][0]];
		ymin = X2[dataX2[pointer[2]][0]];
		ymax = X2[dataX2[pointer[3]][0]];
	}
	else
	{
		xmin = X1[(int)dataX1_reg[pointer[0]][0]];
		xmax = X1[(int)dataX1_reg[pointer[1]][0]];
		ymin = X2[(int)dataX2_reg[pointer[2]][0]];
		ymax = X2[(int)dataX2_reg[pointer[3]][0]];
	}
}
void CrossNode::getSectors_opt(int *Q, Vector_reg **Q_reg)
{
	this->Q = Q;
	if (regression)
	{
		this->Q_reg = Q_reg;
		for (int m = 0; m < 12 * n_classes; m++)
		{
			Q_reg[m]->clear();
		}
	}
	else
	{
		for (int m = 0; m < 12 * n_classes; m++)
		{
			Q[m] = 0;
		}
	}
	int indx;
	int indy;
	if (P[1] - P[0] <= 1)
	{
		indx = cX1[P[1]];
	}
	else
	{
		indx = cX1[(P[0] + P[1]) / 2];
	}
	if (P[3] - P[2] <= 1)
	{
		indy = cX2[P[3]];
	}
	else
	{
		indy = cX2[(P[2] + P[3]) / 2];
	}
	if (!regression)
	{
		x = X1[dataX1[indx][0]];
		y = X2[dataX2[indy][0]];
		x = (X1[dataX1[indx][0]] + X1[dataX1[indx + 1][0]]) / 2;
		y = (X2[dataX2[indy][0]] + X2[dataX2[indy + 1][0]]) / 2;
	}
	else
	{
		x = X1[(int)dataX1_reg[indx][0]];
		y = X2[(int)dataX2_reg[indy][0]];
		x = (X1[(int)dataX1_reg[indx][0]] + X1[(int)dataX1_reg[indx + 1][0]]) / 2;
		y = (X2[(int)dataX2_reg[indy][0]] + X2[(int)dataX2_reg[indy + 1][0]]) / 2;
	}
	for (int i = pointer[0]; i < indx + 1; i++)
	{
		double d;
		if (!regression)
		{
			d = X2[dataX1[i][0]];
		}
		else
		{
			d = X2[(int)dataX1_reg[i][0]];
		}
		if (d <= y)
		{
			if (d < ymin)
			{
				if (!regression)
				{
					Q[0 * n_classes + dataX1[i][1]] += 1;
				}
				else
				{
					Q_reg[0]->push_back(dataX1_reg[i][1]);
				}
			}
			else
			{
				if (!regression)
				{
					Q[3 * n_classes + dataX1[i][1]] += 1;
				}
				else
				{
					Q_reg[3]->push_back(dataX1_reg[i][1]);
				}
			}
		}
		else
		{
			if (d <= ymax)
			{
				if (!regression)
				{
					Q[7 * n_classes + dataX1[i][1]] += 1;
				}
				else
				{
					Q_reg[7]->push_back(dataX1_reg[i][1]);
				}
			}
			else
			{
				if (!regression)
				{
					Q[10 * n_classes + dataX1[i][1]] += 1;
				}
				else
				{
					Q_reg[10]->push_back(dataX1_reg[i][1]);
				}
			}
		}
	}
	for (int i = indx + 1; i < pointer[1] + 1; i++)
	{
		double d;
		if (!regression)
		{
			d = X2[dataX1[i][0]];
		}
		else
		{
			d = X2[(int)dataX1_reg[i][0]];
		}
		if (d <= y)
		{
			if (d < ymin)
			{
				if (!regression)
				{
					Q[1 * n_classes + dataX1[i][1]] += 1;
				}
				else
				{
					Q_reg[1]->push_back(dataX1_reg[i][1]);
				}
			}
			else
			{
				if (!regression)
				{
					Q[4 * n_classes + dataX1[i][1]] += 1;
				}
				else
				{
					Q_reg[4]->push_back(dataX1_reg[i][1]);
				}
			}
		}
		else
		{
			if (d <= ymax)
			{
				if (!regression)
				{
					Q[8 * n_classes + dataX1[i][1]] += 1;
				}
				else
				{
					Q_reg[8]->push_back(dataX1_reg[i][1]);
				}
			}
			else
			{
				if (!regression)
				{
					Q[11 * n_classes + dataX1[i][1]] += 1;
				}
				else
				{
					Q_reg[11]->push_back(dataX1_reg[i][1]);
				}
			}
		}
	}
	for (int i = pointer[2]; i < indy + 1; i++)
	{
		double d;
		if (!regression)
		{
			d = X1[dataX2[i][0]];
		}
		else
		{
			d = X1[(int)dataX2_reg[i][0]];
		}
		if (d < xmin)
		{
			if (!regression)
			{
				Q[2 * n_classes + dataX2[i][1]] += 1;
			}
			else
			{
				Q_reg[2]->push_back(dataX2_reg[i][1]);
			}
		}
		else if (d > xmax)
		{
			if (!regression)
			{
				Q[5 * n_classes + dataX2[i][1]] += 1;
			}
			else
			{
				Q_reg[5]->push_back(dataX2_reg[i][1]);
			}
		}
	}
	for (int i = indy + 1; i < pointer[3] + 1; i++)
	{
		double d;
		if (!regression)
		{
			d = X1[dataX2[i][0]];
		}
		else
		{
			d = X1[(int)dataX2_reg[i][0]];
		}
		if (d < xmin)
		{
			if (!regression)
			{
				Q[6 * n_classes + dataX2[i][1]] += 1;
			}
			else
			{
				Q_reg[6]->push_back(dataX2_reg[i][1]);
			}
		}
		else if (d > xmax)
		{
			if (!regression)
			{
				Q[9 * n_classes + dataX2[i][1]] += 1;
			}
			else
			{
				Q_reg[9]->push_back(dataX2_reg[i][1]);
			}
		}
	}
}
void CrossNode::upperBound()
{
	if (!regression)
	{
		UB = 0;
		vector<double> Ns(n_classes, 0);
		for (int k = 0; k < 4; k++)
		{
			for (int m = 0; m < n_classes; m++)
			{
				if (k == 0)
				{
					Ns[m] = N[0 * n_classes + m] + Q[0 * n_classes + m] + Q[2 * n_classes + m] + Q[3 * n_classes + m];
				}
				else if (k == 1)
				{
					Ns[m] = N[1 * n_classes + m] + Q[1 * n_classes + m] + Q[4 * n_classes + m] + Q[5 * n_classes + m];
				}
				else if (k == 2)
				{
					Ns[m] = N[2 * n_classes + m] + Q[6 * n_classes + m] + Q[7 * n_classes + m] + Q[10 * n_classes + m];
				}
				else if (k == 3)
				{
					Ns[m] = N[3 * n_classes + m] + Q[8 * n_classes + m] + Q[9 * n_classes + m] + Q[11 * n_classes + m];
				}
			}
			UB += imp(Ns, n_classes, criterion);
		}
	}
	else
	{
		UB = 0;
		Vector_reg reg_tmp;
		if (criterion == Environment::mse)
		{
			reg_tmp.keep_vec = false;
		}
		else
		{
			reg_tmp.use_mae();
		}
		reg_tmp.reserve(total_points);
		for (int k = 0; k < 4; k++)
		{
			if (k == 0)
			{
				reg_tmp.add(N_reg[0]);
				reg_tmp.add(Q_reg[0]);
				reg_tmp.add(Q_reg[2]);
				reg_tmp.add(Q_reg[3]);
			}
			else if (k == 1)
			{
				reg_tmp.clear();
				reg_tmp.add(N_reg[1]);
				reg_tmp.add(Q_reg[1]);
				reg_tmp.add(Q_reg[4]);
				reg_tmp.add(Q_reg[5]);
			}
			else if (k == 2)
			{
				reg_tmp.clear();
				reg_tmp.add(N_reg[2]);
				reg_tmp.add(Q_reg[6]);
				reg_tmp.add(Q_reg[7]);
				reg_tmp.add(Q_reg[10]);
			}
			else if (k == 3)
			{
				reg_tmp.clear();
				reg_tmp.add(N_reg[3]);
				reg_tmp.add(Q_reg[8]);
				reg_tmp.add(Q_reg[9]);
				reg_tmp.add(Q_reg[11]);
			}
			UB += reg_tmp.error(criterion);
		}
	}
}
void CrossNode::simpleBound()
{
	if (!regression)
	{
		LB = 0;
		for (int k = 0; k < 4; k++)
		{
			LB += imp1(N + k * n_classes, n_classes, criterion);
		}
	}
	else
	{
		LB = 0;
		for (int k = 0; k < 4; k++)
		{
			LB += N_reg[k]->error(criterion);
		}
	}
}
void CrossNode::getPoint(int countX1, int countX2)
{
	x = xmin;
	y = ymin;
	if (!regression)
	{
		if (pointer[1] < total_points - 1)
		{
			x = (X1[dataX1[pointer[1]][0]] + X1[dataX1[pointer[1] + 1][0]]) / 2;
		}
		if (pointer[3] < total_points - 1)
		{
			y = (X2[dataX2[pointer[3]][0]] + X2[dataX2[pointer[3] + 1][0]]) / 2;
		}
	}
	else
	{
		if (pointer[1] < total_points - 1)
		{
			x = (X1[(int)dataX1_reg[pointer[1]][0]] + X1[(int)dataX1_reg[pointer[1] + 1][0]]) / 2;
		}
		if (pointer[3] < total_points - 1)
		{
			y = (X2[(int)dataX2_reg[pointer[3]][0]] + X2[(int)dataX2_reg[pointer[3] + 1][0]]) / 2;
		}
	}
	if (!regression)
	{
		for (int i = pointer[2]; i < pointer[3] + 1; i++)
		{
			double d = X1[dataX2[i][0]];
			if (d < xmin)
			{
				N[0 * n_classes + dataX2[i][1]]++;
			}
			else if (d > xmax)
			{
				N[1 * n_classes + dataX2[i][1]]++;
			}
		}
		for (int i = pointer[0]; i < pointer[1] + 1; i++)
		{
			double d = X2[dataX1[i][0]];
			if (d > ymin)
			{
				N[2 * n_classes + dataX1[i][1]]++;
			}
			else
			{
				N[0 * n_classes + dataX1[i][1]]++;
			}
		}
	}
	else
	{
		for (int i = pointer[2]; i < pointer[3] + 1; i++)
		{
			double d = X1[(int)dataX2_reg[i][0]];
			if (d < xmin)
			{
				N_reg[0]->push_back(dataX2_reg[i][1]);
			}
			else if (d > xmax)
			{
				N_reg[1]->push_back(dataX2_reg[i][1]);
			}
		}
		for (int i = pointer[0]; i < pointer[1] + 1; i++)
		{
			double d = X2[(int)dataX1_reg[i][0]];
			if (d > ymin)
			{
				N_reg[2]->push_back(dataX1_reg[i][1]);
			}
			else
			{
				N_reg[0]->push_back(dataX1_reg[i][1]);
			}
		}
	}
	UB = getImpN();
}
void CrossNode::getBestXonLine(int countY)
{
	if (!regression)
	{
		for (int i = pointer[2]; i < pointer[3] + 1; i++)
		{
			double d = X1[dataX2[i][0]];
			if (d < xmin)
			{
				N[0 * n_classes + dataX2[i][1]]++;
			}
			else if (d > xmax)
			{
				N[1 * n_classes + dataX2[i][1]]++;
			}
		}
	}
	else
	{
		for (int i = pointer[2]; i < pointer[3] + 1; i++)
		{
			double d = X1[(int)dataX2_reg[i][0]];
			if (d < xmin)
			{
				N_reg[0]->push_back(dataX2_reg[i][1]);
			}
			else if (d > xmax)
			{
				N_reg[1]->push_back(dataX2_reg[i][1]);
			}
		}
	}
	if (!regression)
	{
		for (int i = pointer[0]; i < pointer[1] + 1; i++)
		{
			double d = X2[dataX1[i][0]];
			if (d > ymin)
			{
				N[3 * n_classes + dataX1[i][1]]++;
			}
			else
			{
				N[1 * n_classes + dataX1[i][1]]++;
			}
		}
	}
	else
	{
		for (int i = pointer[0]; i < pointer[1] + 1; i++)
		{
			double d = X2[(int)dataX1_reg[i][0]];
			if (d > ymin)
			{
				N_reg[3]->push_back(dataX1_reg[i][1]);
			}
			else
			{
				N_reg[1]->push_back(dataX1_reg[i][1]);
			}
		}
	}
	double bestSolution = std::numeric_limits<double>::max();
	int indy = cX2[P[3]];
	double y_line;
	if (!regression)
	{
		y_line = X2[dataX2[indy][0]];
		if (indy + 1 < total_points)
		{
			y_line = (X2[dataX2[indy][0]] + X2[dataX2[indy + 1][0]]) / 2;
		}
		else
		{
			y_line = X2[dataX2[indy][0]];
		}
	}
	else
	{
		y_line = X2[(int)dataX2_reg[indy][0]];
		if (indy + 1 < total_points)
		{
			y_line = (X2[(int)dataX2_reg[indy][0]] + X2[(int)dataX2_reg[indy + 1][0]]) / 2;
		}
		else
		{
			y_line = X2[(int)dataX2_reg[indy][0]];
		}
	}
	double bestX;
	double bestY = y_line;
	for (int i = pointer[0]; i < pointer[1] + 1;)
	{
		if (!regression)
		{
			x = X1[dataX1[i][0]];
			while (i < pointer[1] + 1 && X1[dataX1[i][0]] == x)
			{
				y = X2[dataX1[i][0]];
				int m = dataX1[i][1];
				if (y > y_line)
				{
					N[2 * n_classes + m] += 1;
					N[3 * n_classes + m] -= 1;
				}
				else
				{
					N[0 * n_classes + m] += 1;
					N[1 * n_classes + m] -= 1;
				}
				i++;
			}
		}
		else
		{
			x = X1[(int)dataX1_reg[i][0]];
			while (i < pointer[1] + 1 && X1[(int)dataX1_reg[i][0]] == x)
			{
				y = X2[(int)dataX1_reg[i][0]];
				double m = dataX1_reg[i][1];
				if (y > y_line)
				{
					N_reg[2]->push_back(m);
					N_reg[3]->remove(m);
				}
				else
				{
					N_reg[0]->push_back(m);
					N_reg[1]->remove(m);
				}
				i++;
			}
		}
		double gini = getImpN();
		if (gini < bestSolution)
		{
			if (!regression)
			{
				if (i < total_points)
				{
					bestX = (X1[dataX1[i - 1][0]] + X1[dataX1[i][0]]) / 2;
				}
				else
				{
					bestX = X1[dataX1[i - 1][0]];
					if (indy + 1 == total_points)
					{
						goto jmp;
					}
				}
			}
			else
			{
				if (i < total_points)
				{
					bestX = (X1[(int)dataX1_reg[i - 1][0]] + X1[(int)dataX1_reg[i][0]]) / 2;
				}
				else
				{
					bestX = X1[(int)dataX1_reg[i - 1][0]];
					if (indy + 1 == total_points)
					{
						goto jmp;
					}
				}
			}
			bestSolution = gini;
		jmp:;
		}
	}
	x = bestX;
	y = bestY;
	UB = bestSolution;
}
void CrossNode::getBestYonLine(int countX)
{
	if (!regression)
	{
		for (int i = pointer[0]; i < pointer[1] + 1; i++)
		{
			double d = X2[dataX1[i][0]];
			if (d < ymin)
			{
				N[0 * n_classes + dataX1[i][1]]++;
			}
			else if (d > ymax)
			{
				N[2 * n_classes + dataX1[i][1]]++;
			}
		}
	}
	else
	{
		for (int i = pointer[0]; i < pointer[1] + 1; i++)
		{
			double d = X2[(int)dataX1_reg[i][0]];
			if (d < ymin)
			{
				N_reg[0]->push_back(dataX1_reg[i][1]);
			}
			else if (d > ymax)
			{
				N_reg[2]->push_back(dataX1_reg[i][1]);
			}
		}
	}
	if (!regression)
	{
		for (int i = pointer[2]; i < pointer[3] + 1; i++)
		{
			double d = X1[dataX2[i][0]];
			if (d > xmin)
			{
				N[3 * n_classes + dataX2[i][1]]++;
			}
			else
			{
				N[2 * n_classes + dataX2[i][1]]++;
			}
		}
	}
	else
	{
		for (int i = pointer[2]; i < pointer[3] + 1; i++)
		{
			double d = X1[(int)dataX2_reg[i][0]];
			if (d > xmin)
			{
				N_reg[3]->push_back(dataX2_reg[i][1]);
			}
			else
			{
				N_reg[2]->push_back(dataX2_reg[i][1]);
			}
		}
	}
	double bestSolution = std::numeric_limits<double>::max();
	int indx = cX1[P[1]];
	double x_line;
	if (!regression)
	{
		x_line = X1[dataX1[indx][0]];
		if (indx + 1 < total_points)
		{
			x_line = (X1[dataX1[indx][0]] + X1[dataX1[indx + 1][0]]) / 2;
		}
		else
		{
			x_line = X1[dataX1[indx][0]];
		}
	}
	else
	{
		x_line = X1[(int)dataX1_reg[indx][0]];
		if (indx + 1 < total_points)
		{
			x_line = (X1[(int)dataX1_reg[indx][0]] + X1[(int)dataX1_reg[indx + 1][0]]) / 2;
		}
		else
		{
			x_line = X1[(int)dataX1_reg[indx][0]];
		}
	}
	double bestY;
	double bestX = x_line;
	;
	for (int i = pointer[2]; i < pointer[3] + 1;)
	{
		if (!regression)
		{
			y = X2[dataX2[i][0]];
			while (i < pointer[3] + 1 && X2[dataX2[i][0]] == y)
			{
				x = X1[dataX2[i][0]];
				int m = dataX2[i][1];
				if (x <= x_line)
				{
					N[0 * n_classes + m] += 1;
					N[2 * n_classes + m] -= 1;
				}
				else
				{
					N[1 * n_classes + m] += 1;
					N[3 * n_classes + m] -= 1;
				}
				i++;
			}
		}
		else
		{
			y = X2[(int)dataX2_reg[i][0]];
			while (i < pointer[3] + 1 && X2[(int)dataX2_reg[i][0]] == y)
			{
				x = X1[(int)dataX2_reg[i][0]];
				double m = dataX2_reg[i][1];
				if (x <= x_line)
				{
					N_reg[0]->push_back(m);
					N_reg[2]->remove(m);
				}
				else
				{
					N_reg[1]->push_back(m);
					N_reg[3]->remove(m);
				}
				i++;
			}
		}
		double gini = getImpN();
		if (gini < bestSolution)
		{
			if (!regression)
			{
				if (i != total_points)
				{
					bestY = (X2[dataX2[i - 1][0]] + X2[dataX2[i][0]]) / 2;
				}
				else
				{
					bestY = X2[dataX2[i - 1][0]];
					if (indx + 1 == total_points)
					{
						goto jmp;
					}
				}
			}
			else
			{
				if (i != total_points)
				{
					bestY = (X2[(int)dataX2_reg[i - 1][0]] + X2[(int)dataX2_reg[i][0]]) / 2;
				}
				else
				{
					bestY = X2[(int)dataX2_reg[i - 1][0]];
					if (indx + 1 == total_points)
					{
						goto jmp;
					}
				}
			}
			bestSolution = gini;
		jmp:;
		}
	}
	x = bestX;
	y = bestY;
	UB = bestSolution;
}
double CrossNode::getImpN()
{
	double imp = 0;
	if (!regression)
	{
		for (int k = 0; k < 4; k++)
		{
			imp += imp1(N + k * n_classes, n_classes, criterion);
		}
	}
	else
	{
		for (int k = 0; k < 4; k++)
		{
			imp += N_reg[k]->error(criterion);
		}
	}
	return imp;
}
bool CrossNode::checkN()
{
	for (int k = 0; k < 4; k++)
	{
		for (int m = 0; m < n_classes; m++)
		{
			if (N[k * n_classes + m] < 0)
			{
				return false;
			}
		}
	}
	return true;
}
