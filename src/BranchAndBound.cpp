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

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <queue>
#include <stack>
#include <deque>
#include <mutex>
#include <math.h>
#include <cmath>
#include "BBNode.h"
#include "PointValuePair.h"
#include "Result.h"
#include "Environment.h"
using namespace std;
std::mutex mtx;

bool verify(int total_points, bool regression, double *X, double X2[], double X1[], int **dataX1, int **dataX2, double **dataX1_reg, double **dataX2_reg, double a, double b, double d_x, double d_y, double D_x, double D_y)
{
	bool above = false;
	bool below = false;
	double tmp;
	bool splits1 = false;
	bool splits2 = false;
	for (int i = 0; i < total_points; i++)
	{
		tmp = X[total_points + i] - (a * X[i] + b);
		if (tmp >= 0)
		{
			above = true;
		}
		else
		{
			below = true;
		}
		if (above && below)
		{
			splits1 = true;
			break;
		}
	}
	double a2 = a * (D_y / D_x);
	double b2 = D_y * b + d_y - a * d_x * (D_y / D_x);
	above = false;
	below = false;
	for (int i = 0; i < total_points; i++)
	{
		if (!regression)
		{
			tmp = X2[dataX1[i][0]] - (a2 * X1[dataX1[i][0]] + b2);
		}
		else
		{
			tmp = X2[(int)dataX1_reg[i][0]] - (a2 * X1[(int)dataX1_reg[i][0]] + b2);
		}
		if (tmp >= 0)
		{
			above = true;
		}
		else
		{
			below = true;
		}
		if (above && below)
		{
			splits2 = true;
			break;
		}
	}
	if (splits1 && splits2)
	{
		return true;
	}
	return false;
}
void getData(double X1[], double X2[], int ind1, int ind2, int n_classes, int total_points, int **dataX1, int **dataX2, double *X, int *Y, int *P)
{
	double min_x = std::numeric_limits<double>::max();
	double max_x = std::numeric_limits<double>::min();
	double mean_x = 0;
	double min_y = std::numeric_limits<double>::max();
	double max_y = std::numeric_limits<double>::min();
	double mean_y = 0;
	double median_x = X1[dataX1[total_points / 2][0]];
	double median_y = X2[dataX2[total_points / 2][0]];
	double uq_x;
	double lq_x;
	double uq_y;
	double lq_y;
	if (total_points % 2 == 0)
	{
		if ((total_points / 2) % 2 == 0)
		{
			int ind = (total_points) / 4;
			lq_x = 0.5 * (X1[dataX1[ind - 1][0]] + X1[dataX1[ind][0]]);
			lq_y = 0.5 * (X2[dataX2[ind - 1][0]] + X2[dataX2[ind][0]]);
			uq_x = 0.5 * (X1[dataX1[total_points / 2 + ind - 1][0]] + X1[dataX1[total_points / 2 + ind][0]]);
			uq_y = 0.5 * (X2[dataX2[total_points / 2 + ind - 1][0]] + X2[dataX2[total_points / 2 + ind][0]]);
		}
		else
		{
			lq_x = X1[dataX1[(total_points / 2 - 1) / 2][0]];
			lq_y = X2[dataX2[(total_points / 2 - 1) / 2][0]];
			uq_x = X1[dataX1[(total_points / 2 - 1) / 2 + total_points / 2][0]];
			uq_y = X2[dataX2[(total_points / 2 - 1) / 2 + total_points / 2][0]];
		}
	}
	else
	{
		if (((total_points + 1) / 2) % 2 == 0)
		{
			int ind = ((total_points + 1) / 2) / 2;
			lq_x = 0.5 * (X1[dataX1[ind - 1][0]] + X1[dataX1[ind][0]]);
			lq_y = 0.5 * (X2[dataX2[ind - 1][0]] + X2[dataX2[ind][0]]);
			ind = ((total_points + 1) / 2) / 2 + (total_points + 1) / 2 - 1;
			uq_x = 0.5 * (X1[dataX1[ind - 1][0]] + X1[dataX1[ind][0]]);
			uq_y = 0.5 * (X2[dataX2[ind - 1][0]] + X2[dataX2[ind][0]]);
		}
		else
		{
			int ind = ((total_points + 1) / 2 - 1) / 2;
			lq_x = X1[dataX1[ind][0]];
			lq_y = X2[dataX2[ind][0]];
			ind = ((total_points + 1) / 2 - 1) / 2 + (total_points + 1) / 2 - 1;
			;
			uq_x = X1[dataX1[ind][0]];
			uq_y = X2[dataX2[ind][0]];
		}
	}
	if (total_points % 2 == 0)
	{
		median_x = 0.5 * (X1[dataX1[(total_points) / 2 - 1][0]] + X1[dataX1[(total_points) / 2][0]]);
		median_y = 0.5 * (X2[dataX2[(total_points) / 2 - 1][0]] + X2[dataX2[(total_points) / 2][0]]);
	}
	else
	{
		median_x = X1[dataX1[(total_points - 1) / 2][0]];
		median_y = X2[dataX2[(total_points - 1) / 2][0]];
	}
	for (int i = 0; i < total_points; i++)
	{
		mean_x += X[i];
		if (X[i] < min_x)
		{
			min_x = X[i];
		}
		if (X[i] > max_x)
		{
			max_x = X[i];
		}
		mean_y += X[1 * total_points + i];
		if (X[1 * total_points + i] < min_y)
		{
			min_y = X[1 * total_points + i];
		}
		if (X[1 * total_points + i] > max_y)
		{
			max_y = X[1 * total_points + i];
		}
	}
	double absmax_x = max_x;
	if (abs(min_x) > max_x)
	{
		absmax_x = abs(min_x);
	}
	double absmax_y = max_y;
	if (abs(min_y) > max_y)
	{
		absmax_y = abs(min_y);
	}
	mean_x /= total_points;
	mean_y /= total_points;
	double std_x = 0;
	double std_y = 0;
	for (int i = 0; i < total_points; i++)
	{
		std_x += (X[i] - mean_x) * (X[i] - mean_x);
		std_y += (X[1 * total_points + i] - mean_y) * (X[1 * total_points + i] - mean_y);
	}
	std_x = sqrt(std_x / (total_points - 1));
	std_y = sqrt(std_y / (total_points - 1));
	double D_x = 2 * std_x;
	double D_y = 2 * std_y;
	D_x = uq_x - lq_x;
	D_y = uq_y - lq_y;
	if (D_x == 0)
	{
		D_x = 2 * std_x;
	}
	if (D_y == 0)
	{
		D_y = 2 * std_y;
	}
	double d_x = median_x;
	double d_y = median_y;
	if (D_x == 0)
	{
		D_x = 1;
	}
	if (D_y == 0)
	{
		D_y = 1;
	}
	for (int i = 0; i < total_points; i++)
	{
		X[i] = (X[i] - d_x) / D_x;
		X[1 * total_points + i] = (X[1 * total_points + i] - d_y) / D_y;
	}
}
void branchAndBound(int id, double X1[], double X2[], int ind1, int ind2, int n_classes, int total_points, int **dataX1, int **dataX2, double **dataX1_reg, double **dataX2_reg, Result *res, Environment *env)
{
	bool regression = env->regression;
	int counter[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	double *X = new double[2 * total_points];
	double *Y_reg;
	int *Y;
	if (!regression)
	{
		Y = new int[total_points];
	}
	else
	{
		Y_reg = new double[total_points];
	}
	BBNode *nodes = new BBNode[4];
	double median;
	if (!regression)
	{
		median = X1[dataX1[total_points / 2][0]];
	}
	else
	{
		median = X1[(int)dataX1_reg[total_points / 2][0]];
	}
	double amin = -1;
	double amax = 1;
	double eps = std::numeric_limits<double>::epsilon();
	for (int i = 0; i < total_points; i++)
	{
		if (!regression)
		{
			X[i] = X1[dataX1[i][0]];
			X[1 * total_points + i] = X2[dataX1[i][0]];
			Y[i] = dataX1[i][1];
		}
		else
		{
			X[i] = X1[(int)dataX1_reg[i][0]];
			X[1 * total_points + i] = X2[(int)dataX1_reg[i][0]];
			Y_reg[i] = dataX1_reg[i][1];
		}
	}
	double D_x = 1;
	double D_y = 1;
	double d_x = 0;
	double d_y = 0;
	double min_x = std::numeric_limits<double>::max();
	double max_x = std::numeric_limits<double>::min();
	double min_y = std::numeric_limits<double>::max();
	double max_y = std::numeric_limits<double>::min();
	if (env->normalize)
	{
		double mean_x = 0;
		double mean_y = 0;
		double median_x;
		double median_y;
		if (!regression)
		{
			median_x = X1[dataX1[total_points / 2][0]];
			median_y = X2[dataX2[total_points / 2][0]];
		}
		else
		{
			median_x = X1[(int)dataX1_reg[total_points / 2][0]];
			median_y = X2[(int)dataX2_reg[total_points / 2][0]];
		}
		double uq_x;
		double lq_x;
		double uq_y;
		double lq_y;
		if (!regression)
		{
			if (total_points % 2 == 0)
			{
				if ((total_points / 2) % 2 == 0)
				{
					int ind = (total_points) / 4;
					lq_x = 0.5 * (X1[dataX1[ind - 1][0]] + X1[dataX1[ind][0]]);
					lq_y = 0.5 * (X2[dataX2[ind - 1][0]] + X2[dataX2[ind][0]]);
					uq_x = 0.5 * (X1[dataX1[total_points / 2 + ind - 1][0]] + X1[dataX1[total_points / 2 + ind][0]]);
					uq_y = 0.5 * (X2[dataX2[total_points / 2 + ind - 1][0]] + X2[dataX2[total_points / 2 + ind][0]]);
				}
				else
				{
					lq_x = X1[dataX1[(total_points / 2 - 1) / 2][0]];
					lq_y = X2[dataX2[(total_points / 2 - 1) / 2][0]];
					uq_x = X1[dataX1[(total_points / 2 - 1) / 2 + total_points / 2][0]];
					uq_y = X2[dataX2[(total_points / 2 - 1) / 2 + total_points / 2][0]];
				}
			}
			else
			{
				if (((total_points + 1) / 2) % 2 == 0)
				{
					int ind = ((total_points + 1) / 2) / 2;
					lq_x = 0.5 * (X1[dataX1[ind - 1][0]] + X1[dataX1[ind][0]]);
					lq_y = 0.5 * (X2[dataX2[ind - 1][0]] + X2[dataX2[ind][0]]);
					ind = ((total_points + 1) / 2) / 2 + (total_points + 1) / 2 - 1;
					uq_x = 0.5 * (X1[dataX1[ind - 1][0]] + X1[dataX1[ind][0]]);
					uq_y = 0.5 * (X2[dataX2[ind - 1][0]] + X2[dataX2[ind][0]]);
				}
				else
				{
					int ind = ((total_points + 1) / 2 - 1) / 2;
					lq_x = X1[dataX1[ind][0]];
					lq_y = X2[dataX2[ind][0]];
					ind = ((total_points + 1) / 2 - 1) / 2 + (total_points + 1) / 2 - 1;
					;
					uq_x = X1[dataX1[ind][0]];
					uq_y = X2[dataX2[ind][0]];
				}
			}
			if (total_points % 2 == 0)
			{
				median_x = 0.5 * (X1[dataX1[(total_points) / 2 - 1][0]] + X1[dataX1[(total_points) / 2][0]]);
				median_y = 0.5 * (X2[dataX2[(total_points) / 2 - 1][0]] + X2[dataX2[(total_points) / 2][0]]);
			}
			else
			{
				median_x = X1[dataX1[(total_points - 1) / 2][0]];
				median_y = X2[dataX2[(total_points - 1) / 2][0]];
			}
		}
		else
		{
			if (total_points % 2 == 0)
			{
				if ((total_points / 2) % 2 == 0)
				{
					int ind = (total_points) / 4;
					lq_x = 0.5 * (X1[(int)dataX1_reg[ind - 1][0]] + X1[(int)dataX1_reg[ind][0]]);
					lq_y = 0.5 * (X2[(int)dataX2_reg[ind - 1][0]] + X2[(int)dataX2_reg[ind][0]]);
					uq_x = 0.5 * (X1[(int)dataX1_reg[total_points / 2 + ind - 1][0]] + X1[(int)dataX1_reg[total_points / 2 + ind][0]]);
					uq_y = 0.5 * (X2[(int)dataX2_reg[total_points / 2 + ind - 1][0]] + X2[(int)dataX2_reg[total_points / 2 + ind][0]]);
				}
				else
				{
					lq_x = X1[(int)dataX1_reg[(total_points / 2 - 1) / 2][0]];
					lq_y = X2[(int)dataX2_reg[(total_points / 2 - 1) / 2][0]];
					uq_x = X1[(int)dataX1_reg[(total_points / 2 - 1) / 2 + total_points / 2][0]];
					uq_y = X2[(int)dataX2_reg[(total_points / 2 - 1) / 2 + total_points / 2][0]];
				}
			}
			else
			{
				if (((total_points + 1) / 2) % 2 == 0)
				{
					int ind = ((total_points + 1) / 2) / 2;
					lq_x = 0.5 * (X1[(int)dataX1_reg[ind - 1][0]] + X1[(int)dataX1_reg[ind][0]]);
					lq_y = 0.5 * (X2[(int)dataX2_reg[ind - 1][0]] + X2[(int)dataX2_reg[ind][0]]);
					ind = ((total_points + 1) / 2) / 2 + (total_points + 1) / 2 - 1;
					uq_x = 0.5 * (X1[(int)dataX1_reg[ind - 1][0]] + X1[(int)dataX1_reg[ind][0]]);
					uq_y = 0.5 * (X2[(int)dataX2_reg[ind - 1][0]] + X2[(int)dataX2_reg[ind][0]]);
				}
				else
				{
					int ind = ((total_points + 1) / 2 - 1) / 2;
					lq_x = X1[(int)dataX1_reg[ind][0]];
					lq_y = X2[(int)dataX2_reg[ind][0]];
					ind = ((total_points + 1) / 2 - 1) / 2 + (total_points + 1) / 2 - 1;
					;
					uq_x = X1[(int)dataX1_reg[ind][0]];
					uq_y = X2[(int)dataX2_reg[ind][0]];
				}
			}
			if (total_points % 2 == 0)
			{
				median_x = 0.5 * (X1[(int)dataX1_reg[(total_points) / 2 - 1][0]] + X1[(int)dataX1_reg[(total_points) / 2][0]]);
				median_y = 0.5 * (X2[(int)dataX2_reg[(total_points) / 2 - 1][0]] + X2[(int)dataX2_reg[(total_points) / 2][0]]);
			}
			else
			{
				median_x = X1[(int)dataX1_reg[(total_points - 1) / 2][0]];
				median_y = X2[(int)dataX2_reg[(total_points - 1) / 2][0]];
			}
		}
		for (int i = 0; i < total_points; i++)
		{
			mean_x += X[i];
			if (X[i] < min_x)
			{
				min_x = X[i];
			}
			if (X[i] > max_x)
			{
				max_x = X[i];
			}
			mean_y += X[1 * total_points + i];
			if (X[1 * total_points + i] < min_y)
			{
				min_y = X[1 * total_points + i];
			}
			if (X[1 * total_points + i] > max_y)
			{
				max_y = X[1 * total_points + i];
			}
		}
		double absmax_x = max_x;
		if (abs(min_x) > max_x)
		{
			absmax_x = abs(min_x);
		}
		double absmax_y = max_y;
		if (abs(min_y) > max_y)
		{
			absmax_y = abs(min_y);
		}
		mean_x /= total_points;
		mean_y /= total_points;
		double std_x = 0;
		double std_y = 0;
		for (int i = 0; i < total_points; i++)
		{
			std_x += (X[i] - mean_x) * (X[i] - mean_x);
			std_y += (X[1 * total_points + i] - mean_y) * (X[1 * total_points + i] - mean_y);
		}
		std_x = sqrt(std_x / (total_points - 1));
		std_y = sqrt(std_y / (total_points - 1));
		D_x = 2 * std_x;
		D_y = 2 * std_y;
		D_x = uq_x - lq_x;
		D_y = uq_y - lq_y;
		if (D_x == 0)
		{
			D_x = 2 * std_x;
		}
		if (D_y == 0)
		{
			D_y = 2 * std_y;
		}
		d_x = median_x;
		d_y = median_y;
		if (median_x == max_x || median_x == min_x)
		{
			d_x = mean_x;
		}
		if (median_y == max_y || median_y == min_y)
		{
			d_y = mean_y;
		}
		if (D_x == 0)
		{
			D_x = 1;
		}
		if (D_y == 0)
		{
			D_y = 1;
		}
		for (int i = 0; i < total_points; i++)
		{
			X[i] = (X[i] - d_x) / D_x;
			X[1 * total_points + i] = (X[1 * total_points + i] - d_y) / D_y;
		}
	}
	else
	{
		D_x = 1;
		D_y = 1;
		d_x = 0;
		d_y = 0;
		for (int i = 0; i < total_points; i++)
		{
			X[i] = (X[i] - d_x) / D_x;
			X[1 * total_points + i] = (X[1 * total_points + i] - d_y) / D_y;
		}
	}
	double bmin = std::numeric_limits<double>::max();
	double bmax = std::numeric_limits<double>::min();
	for (int i = 0; i < total_points; i++)
	{
		double tmp = X[1 * total_points + i] - amin * X[i];
		if (tmp > bmax)
		{
			bmax = tmp;
		}
		if (tmp < bmin)
		{
			bmin = tmp;
		}
		tmp = X[1 * total_points + i] - amax * X[i];
		if (tmp > bmax)
		{
			bmax = tmp;
		}
		if (tmp < bmin)
		{
			bmin = tmp;
		}
	}
	bmin -= eps;
	bmax += eps;
	double *O;
	double *o;
	double *o1;
	double *o2;
	if (!regression)
	{
		O = new double[3 * n_classes];
		for (int m = 0; m < n_classes; m++)
		{
			O[m] = 0;
			O[n_classes + m] = 0;
			O[2 * n_classes + m] = 0;
		}
		for (int i = 0; i < total_points; i++)
		{
			O[2 * n_classes + Y[i]]++;
		}
		o = new double[2 * n_classes];
		o1 = new double[3 * n_classes];
		o2 = new double[3 * n_classes];
	}
	else
	{
		O = new double[9];
		o = new double[6];
		for (int m = 0; m < 9; m++)
		{
			O[m] = 0;
			if (m < 6)
			{
				o[m] = 0;
			}
		}
	}
	doubleIntPair *dists = new doubleIntPair[2 * total_points];
	deque<BBNode> que;
	double *newX = new double[2 * total_points];
	int *newY;
	double *newY_reg;
	Vector_reg *new_reg1;
	Vector_reg *new_reg2;
	Vector_reg *new_reg3;
	if (!regression)
	{
		newY = new int[total_points];
	}
	else
	{
		newY_reg = new double[total_points];
		if (env->criterion == Environment::mae)
		{
			new_reg1 = new Vector_reg();
			new_reg2 = new Vector_reg();
			new_reg3 = new Vector_reg();
			new_reg1->use_mae();
			new_reg2->use_mae();
			new_reg3->use_mae();
		}
	}
	for (int i = 0; i < total_points; i++)
	{
		newX[2 * i] = X[i];
		newX[2 * i + 1] = X[total_points + i];
		if (!regression)
		{
			newY[i] = Y[i];
		}
		else
		{
			newY_reg[i] = Y_reg[i];
			if (env->criterion == Environment::mae)
			{
				new_reg3->push_back(Y_reg[i]);
			}
		}
	}
	BBNode root;
	if (!regression)
	{
		int type = 3;
		root = BBNode(newX, newY, newY_reg, amin, amax, bmin, bmax, O, o, new_reg1, new_reg2, new_reg3, total_points, false, 0, n_classes, total_points, dists, type, regression, env->criterion);
	}
	else
	{
		int type = 3;
		root = BBNode(newX, newY, newY_reg, amin, amax, bmin, bmax, O, o, new_reg1, new_reg2, new_reg3, total_points, false, 0, n_classes, total_points, dists, type, regression, env->criterion);
	}
	que.push_back(root);
	int maxNodes = 0;
	double t_amin = 0;
	double t_amax = 0;
	double t_bmin = 0;
	double t_bmax = 0;
	int t_n_P = 0;
	int times = 10000;
	int g = 0;
	bool switch_Rule;
	int max_free_size = 2048;
	double *freeXs[2048];
	int *freeYs[2048];
	double *freeYs_reg[2048];
	double *freeOs[2048];
	Vector_reg *free_reg1[2048];
	Vector_reg *free_reg2[2048];
	Vector_reg *free_reg3[2048];
	int free_size = 0;
	double bestA = 0;
	double bestB = (bmin + bmax) / 2;
	double bestImp = std::numeric_limits<double>::max();
	int herkunft = -1;
	bool update = false;
	int count = 0;
	BBNode bestNode = root;
	bool verification = false;
	while (que.size() > 0)
	{
		if (que.size() > maxNodes)
		{
			maxNodes = que.size();
		}
		BBNode v = que.back();
		que.pop_back();
		if (v.LB < res->bestImp)
		{
			count++;
			v.process();
			if (v.status == 1)
			{
				if (v.type == 3)
				{
					if (v.amax - v.amin > v.eps || v.bmax - v.bmin > v.eps)
					{
						v.computeBranchingRule();
						if (v.n_P <= 1 * total_points)
						{
							v.upperBound();
						}
					}
					else
					{
						v.branch_a = (v.amin + v.amax) / 2.0;
						v.branch_b = (v.bmin + v.bmax) / 2.0;
						v.upperBound();
					}
				}
				else
				{
					if (v.type == 0)
					{
						if (!regression)
						{
							v.computeBranchingRule_B_Bound(false);
						}
						else
						{
							cout << "NOT IMPLEMENTED";
							exit(0);
						}
					}
					else
					{
						if (!regression)
						{
							v.computeBranchingRule_A_Bound(false);
						}
						else
						{
							cout << "NOT IMPLEMENTED";
							exit(0);
						}
					}
				}
				bool flag = false;
				mtx.lock();
				if (v.UB < res->bestImp)
				{
					bool verified = false;
					if (verification && !verify(total_points, regression, X, X2, X1, dataX1, dataX2, dataX1_reg, dataX2_reg, v.a, v.b, d_x, d_y, D_x, D_y))
					{
						v.upperBound_mid();
						if (v.UB < res->bestImp)
						{
							if (verify(total_points, regression, X, X2, X1, dataX1, dataX2, dataX1_reg, dataX2_reg, v.a, v.b, d_x, d_y, D_x, D_y))
							{
								verified = true;
							}
						}
					}
					else
					{
						verified = true;
					}
					if (verified)
					{
						res->split_found = true;
						res->bestImp = v.UB;
						res->bestA = v.a * (D_y / D_x);
						res->bestB = D_y * v.b + d_y - v.a * d_x * (D_y / D_x);
						res->bestInd1 = ind1;
						res->bestInd2 = ind2;
						flag = true;
					}
				}
				mtx.unlock();
				if (flag)
				{
					bestA = v.a;
					bestB = v.b;
					bestImp = v.UB;
					herkunft = 1;
					update = true;
				}
				if (v.LB < v.UB && (v.amax - v.amin > v.eps || v.bmax - v.bmin > v.eps))
				{
					if (v.type == 3)
					{
						for (int l = 0; l < 4; l++)
						{
							t_n_P = 0;
							if (l == 0)
							{
								t_amin = v.branch_a;
								t_amax = v.amax;
								t_bmin = v.branch_b;
								t_bmax = v.bmax;
							}
							else if (l == 1)
							{
								t_amin = v.amin;
								t_amax = v.branch_a;
								t_bmin = v.branch_b;
								t_bmax = v.bmax;
							}
							else if (l == 2)
							{
								t_amin = v.branch_a;
								t_amax = v.amax;
								t_bmin = v.bmin;
								t_bmax = v.branch_b;
							}
							else
							{
								t_amin = v.amin;
								t_amax = v.branch_a;
								t_bmin = v.bmin;
								t_bmax = v.branch_b;
							}
							double *newO;
							newX;
							newY;
							newY_reg;
							new_reg1;
							new_reg2;
							new_reg3;
							if (free_size > 0)
							{
								newO = freeOs[free_size - 1];
								newX = freeXs[free_size - 1];
								if (!regression)
								{
									newY = freeYs[free_size - 1];
								}
								else
								{
									newY_reg = freeYs_reg[free_size - 1];
									if (env->criterion == Environment::mae)
									{
										new_reg1 = free_reg1[free_size - 1];
										new_reg2 = free_reg2[free_size - 1];
										new_reg3 = free_reg3[free_size - 1];
										new_reg1->clear();
										new_reg2->clear();
										new_reg3->clear();
									}
								}
								free_size--;
							}
							else
							{
								if (!regression)
								{
									newO = new double[3 * n_classes];
									newX = new double[2 * total_points];
									newY = new int[total_points];
								}
								else
								{
									newO = new double[9];
									newX = new double[2 * total_points];
									newY_reg = new double[total_points];
									if (env->criterion == Environment::mae)
									{
										new_reg1 = new Vector_reg();
										new_reg2 = new Vector_reg();
										new_reg3 = new Vector_reg();
										new_reg1->reserve(v.reg1->length);
										new_reg2->reserve(v.reg2->length);
										new_reg3->reserve(v.reg3->length);
										new_reg1->use_mae();
										new_reg2->use_mae();
										new_reg3->use_mae();
									}
								}
							}
							bool same_as_parent = true;
							if (!regression)
							{
								for (int m = 0; m < n_classes; m++)
								{
									newO[m] = v.O[m];
									newO[n_classes + m] = v.O[n_classes + m];
									newO[2 * n_classes + m] = 0;
								}
							}
							else
							{
								for (int m = 0; m < 3; m++)
								{
									newO[m] = v.O[m];
									newO[3 + m] = v.O[3 + m];
									newO[6 + m] = 0;
								}
								if (env->criterion == Environment::mae)
								{
									v.reg1->copy(new_reg1);
									v.reg2->copy(new_reg2);
								}
							}
							for (int i = 0; i < v.n_P; i++)
							{
								if (v.X[2 * i] >= 0)
								{
									if (v.X[2 * i + 1] >= t_amax * v.X[2 * i] + t_bmax)
									{
										if (!regression)
										{
											newO[v.Y[i]] += 1;
										}
										else
										{
											newO[0] += v.Y_reg[i];
											newO[1] += v.Y_reg[i] * v.Y_reg[i];
											newO[2] += 1;
											if (env->criterion == Environment::mae)
											{
												new_reg1->push_back(v.Y_reg[i]);
											}
										}
										same_as_parent = false;
									}
									else if (v.X[2 * i + 1] < t_amin * v.X[2 * i] + t_bmin)
									{
										if (!regression)
										{
											newO[n_classes + v.Y[i]] += 1;
										}
										else
										{
											newO[3] += v.Y_reg[i];
											newO[4] += v.Y_reg[i] * v.Y_reg[i];
											newO[5] += 1;
											if (env->criterion == Environment::mae)
											{
												new_reg2->push_back(v.Y_reg[i]);
											}
										}
										same_as_parent = false;
									}
									else
									{
										newX[2 * t_n_P] = v.X[2 * i];
										newX[2 * t_n_P + 1] = v.X[2 * i + 1];
										if (!regression)
										{
											newY[t_n_P] = v.Y[i];
										}
										else
										{
											newY_reg[t_n_P] = v.Y_reg[i];
										}
										if (!regression)
										{
											newO[2 * n_classes + v.Y[i]] += 1;
										}
										else
										{
											newO[6] += v.Y_reg[i];
											newO[7] += v.Y_reg[i] * v.Y_reg[i];
											newO[8] += 1;
											if (env->criterion == Environment::mae)
											{
												new_reg3->push_back(v.Y_reg[i]);
											}
										}
										t_n_P += 1;
									}
								}
								else
								{
									if (v.X[1 * 2 * i + 1] >= t_amin * v.X[2 * i] + t_bmax)
									{
										if (!regression)
										{
											newO[v.Y[i]] += 1;
										}
										else
										{
											newO[0] += v.Y_reg[i];
											newO[1] += v.Y_reg[i] * v.Y_reg[i];
											newO[2] += 1;
											if (env->criterion == Environment::mae)
											{
												new_reg1->push_back(v.Y_reg[i]);
											}
										}
										same_as_parent = false;
									}
									else if (v.X[1 * 2 * i + 1] < t_amax * v.X[2 * i] + t_bmin)
									{
										if (!regression)
										{
											newO[n_classes + v.Y[i]] += 1;
										}
										else
										{
											newO[3] += v.Y_reg[i];
											newO[4] += v.Y_reg[i] * v.Y_reg[i];
											newO[5] += 1;
											if (env->criterion == Environment::mae)
											{
												new_reg2->push_back(v.Y_reg[i]);
											}
										}
										same_as_parent = false;
									}
									else
									{
										newX[2 * t_n_P] = v.X[2 * i];
										newX[2 * t_n_P + 1] = v.X[2 * i + 1];
										if (!regression)
										{
											newY[t_n_P] = v.Y[i];
										}
										else
										{
											newY_reg[t_n_P] = v.Y_reg[i];
										}
										if (!regression)
										{
											newO[2 * n_classes + v.Y[i]] += 1;
										}
										else
										{
											newO[6] += v.Y_reg[i];
											newO[7] += v.Y_reg[i] * v.Y_reg[i];
											newO[8] += 1;
											if (env->criterion == Environment::mae)
											{
												new_reg3->push_back(v.Y_reg[i]);
											}
										}
										t_n_P += 1;
									}
								}
							}
							int ind = 0;
							BBNode node = BBNode(newX, newY, newY_reg, t_amin, t_amax, t_bmin, t_bmax, newO, o, new_reg1, new_reg2, new_reg3, t_n_P, same_as_parent, v.LB, n_classes, total_points, dists, 3, regression, env->criterion);
							while (ind < l && node.LB > nodes[ind].LB)
							{
								ind++;
							}
							for (int ind2 = l - 1; ind2 >= ind; ind2--)
							{
								nodes[ind2 + 1] = nodes[ind2];
							}
							nodes[ind] = node;
						}
						for (int l = 3; l >= 0; l--)
						{
							que.push_back(nodes[l]);
						}
					}
					else if (v.type == 0)
					{
						for (int l = 0; l < 2; l++)
						{
							t_n_P = 0;
							if (l == 0)
							{
								t_amin = v.amin;
								t_amax = v.amax;
								t_bmin = v.branch_b;
								t_bmax = v.bmax;
							}
							else if (l == 1)
							{
								t_amin = v.amin;
								t_amax = v.amax;
								t_bmin = v.bmin;
								t_bmax = v.branch_b;
							}
							double *newO;
							newX;
							newY;
							new_reg1;
							new_reg2;
							new_reg3;
							if (free_size > 0)
							{
								newO = freeOs[free_size - 1];
								newX = freeXs[free_size - 1];
								if (!regression)
								{
									newY = freeYs[free_size - 1];
								}
								else
								{
									newY_reg = freeYs_reg[free_size - 1];
									if (env->criterion == Environment::mae)
									{
										new_reg1 = free_reg1[free_size - 1];
										new_reg2 = free_reg2[free_size - 1];
										new_reg3 = free_reg3[free_size - 1];
										new_reg1->clear();
										new_reg2->clear();
										new_reg3->clear();
									}
								}
								free_size--;
							}
							else
							{
								if (!regression)
								{
									newO = new double[3 * n_classes];
									newX = new double[2 * total_points];
									newY = new int[total_points];
								}
								else
								{
									newO = new double[9];
									newX = new double[2 * total_points];
									newY_reg = new double[total_points];
									if (env->criterion == Environment::mae)
									{
										new_reg1 = new Vector_reg();
										new_reg2 = new Vector_reg();
										new_reg3 = new Vector_reg();
										new_reg1->reserve(v.reg1->length);
										new_reg2->reserve(v.reg2->length);
										new_reg3->reserve(v.reg3->length);
										new_reg1->use_mae();
										new_reg2->use_mae();
										new_reg3->use_mae();
									}
								}
							}
							bool same_as_parent = true;
							if (!regression)
							{
								for (int m = 0; m < n_classes; m++)
								{
									newO[m] = v.O[m];
									newO[n_classes + m] = v.O[n_classes + m];
									newO[2 * n_classes + m] = 0;
								}
							}
							else
							{
								for (int m = 0; m < 3; m++)
								{
									newO[m] = v.O[m];
									newO[3 + m] = v.O[3 + m];
									newO[6 + m] = 0;
								}
								if (env->criterion == Environment::mae)
								{
									v.reg1->copy(new_reg1);
									v.reg2->copy(new_reg2);
								}
							}
							for (int i = 0; i < v.n_P; i++)
							{
								if (v.X[2 * i] >= 0)
								{
									if (v.X[2 * i + 1] >= t_amax * v.X[2 * i] + t_bmax)
									{
										if (!regression)
										{
											newO[v.Y[i]] += 1;
										}
										else
										{
											newO[0] += v.Y_reg[i];
											newO[1] += v.Y_reg[i] * v.Y_reg[i];
											newO[2] += 1;
											if (env->criterion == Environment::mae)
											{
												new_reg1->push_back(v.Y_reg[i]);
											}
										}
										same_as_parent = false;
									}
									else if (v.X[2 * i + 1] < t_amin * v.X[2 * i] + t_bmin)
									{
										if (!regression)
										{
											newO[n_classes + v.Y[i]] += 1;
										}
										else
										{
											newO[3] += v.Y_reg[i];
											newO[4] += v.Y_reg[i] * v.Y_reg[i];
											newO[5] += 1;
											if (env->criterion == Environment::mae)
											{
												new_reg2->push_back(v.Y_reg[i]);
											}
										}
										same_as_parent = false;
									}
									else
									{
										newX[2 * t_n_P] = v.X[2 * i];
										newX[2 * t_n_P + 1] = v.X[2 * i + 1];
										if (!regression)
										{
											newY[t_n_P] = v.Y[i];
										}
										else
										{
											newY_reg[t_n_P] = v.Y_reg[i];
										}
										if (!regression)
										{
											newO[2 * n_classes + v.Y[i]] += 1;
										}
										else
										{
											newO[6] += v.Y_reg[i];
											newO[7] += v.Y_reg[i] * v.Y_reg[i];
											newO[8] += 1;
											if (env->criterion == Environment::mae)
											{
												new_reg3->push_back(v.Y_reg[i]);
											}
										}
										t_n_P += 1;
									}
								}
								else
								{
									if (v.X[1 * 2 * i + 1] >= t_amin * v.X[2 * i] + t_bmax)
									{
										if (!regression)
										{
											newO[v.Y[i]] += 1;
										}
										else
										{
											newO[0] += v.Y_reg[i];
											newO[1] += v.Y_reg[i] * v.Y_reg[i];
											newO[2] += 1;
											if (env->criterion == Environment::mae)
											{
												new_reg1->push_back(v.Y_reg[i]);
											}
										}
										same_as_parent = false;
									}
									else if (v.X[1 * 2 * i + 1] < t_amax * v.X[2 * i] + t_bmin)
									{
										if (!regression)
										{
											newO[n_classes + v.Y[i]] += 1;
										}
										else
										{
											newO[3] += v.Y_reg[i];
											newO[4] += v.Y_reg[i] * v.Y_reg[i];
											newO[5] += 1;
											if (env->criterion == Environment::mae)
											{
												new_reg2->push_back(v.Y_reg[i]);
											}
										}
										same_as_parent = false;
									}
									else
									{
										newX[2 * t_n_P] = v.X[2 * i];
										newX[2 * t_n_P + 1] = v.X[2 * i + 1];
										if (!regression)
										{
											newY[t_n_P] = v.Y[i];
										}
										else
										{
											newY_reg[t_n_P] = v.Y_reg[i];
										}
										if (!regression)
										{
											newO[2 * n_classes + v.Y[i]] += 1;
										}
										else
										{
											newO[6] += v.Y_reg[i];
											newO[7] += v.Y_reg[i] * v.Y_reg[i];
											newO[8] += 1;
											if (env->criterion == Environment::mae)
											{
												new_reg3->push_back(v.Y_reg[i]);
											}
										}
										t_n_P += 1;
									}
								}
							}
							int ind = 0;
							BBNode node = BBNode(newX, newY, newY_reg, t_amin, t_amax, t_bmin, t_bmax, newO, o, new_reg1, new_reg2, new_reg3, t_n_P, same_as_parent, v.LB, n_classes, total_points, dists, (v.type + 1) % 2, regression, env->criterion);
							while (ind < l && node.LB > nodes[ind].LB)
							{
								ind++;
							}
							for (int ind2 = l - 1; ind2 >= ind; ind2--)
							{
								nodes[ind2 + 1] = nodes[ind2];
							}
							nodes[ind] = node;
						}
						for (int l = 1; l >= 0; l--)
						{
							que.push_back(nodes[l]);
						}
					}
					else
					{
						for (int l = 0; l < 2; l++)
						{
							t_n_P = 0;
							if (l == 0)
							{
								t_amin = v.branch_a;
								t_amax = v.amax;
								t_bmin = v.bmin;
								t_bmax = v.bmax;
							}
							else
							{
								t_amin = v.amin;
								t_amax = v.branch_a;
								t_bmin = v.bmin;
								t_bmax = v.bmax;
							}
							double *newO;
							newX;
							newY;
							new_reg1;
							new_reg2;
							new_reg3;
							if (free_size > 0)
							{
								newO = freeOs[free_size - 1];
								newX = freeXs[free_size - 1];
								if (!regression)
								{
									newY = freeYs[free_size - 1];
								}
								else
								{
									newY_reg = freeYs_reg[free_size - 1];
									if (env->criterion == Environment::mae)
									{
										new_reg1 = free_reg1[free_size - 1];
										new_reg2 = free_reg2[free_size - 1];
										new_reg3 = free_reg3[free_size - 1];
										new_reg1->clear();
										new_reg2->clear();
										new_reg3->clear();
									}
								}
								free_size--;
							}
							else
							{
								if (!regression)
								{
									newO = new double[3 * n_classes];
									newX = new double[2 * total_points];
									newY = new int[total_points];
								}
								else
								{
									newO = new double[9];
									newX = new double[2 * total_points];
									newY_reg = new double[total_points];
									if (env->criterion == Environment::mae)
									{
										new_reg1 = new Vector_reg();
										new_reg2 = new Vector_reg();
										new_reg3 = new Vector_reg();
										new_reg1->reserve(v.reg1->length);
										new_reg2->reserve(v.reg2->length);
										new_reg3->reserve(v.reg3->length);
										new_reg1->use_mae();
										new_reg2->use_mae();
										new_reg3->use_mae();
									}
								}
							}
							bool same_as_parent = true;
							if (!regression)
							{
								for (int m = 0; m < n_classes; m++)
								{
									newO[m] = v.O[m];
									newO[n_classes + m] = v.O[n_classes + m];
									newO[2 * n_classes + m] = 0;
								}
							}
							else
							{
								for (int m = 0; m < 3; m++)
								{
									newO[m] = v.O[m];
									newO[3 + m] = v.O[3 + m];
									newO[6 + m] = 0;
								}
								if (env->criterion == Environment::mae)
								{
									v.reg1->copy(new_reg1);
									v.reg2->copy(new_reg2);
								}
							}
							for (int i = 0; i < v.n_P; i++)
							{
								if (v.X[2 * i] >= 0)
								{
									if (v.X[2 * i + 1] >= t_amax * v.X[2 * i] + t_bmax)
									{
										if (!regression)
										{
											newO[v.Y[i]] += 1;
										}
										else
										{
											newO[0] += v.Y_reg[i];
											newO[1] += v.Y_reg[i] * v.Y_reg[i];
											newO[2] += 1;
											if (env->criterion == Environment::mae)
											{
												new_reg1->push_back(v.Y_reg[i]);
											}
										}
										same_as_parent = false;
									}
									else if (v.X[2 * i + 1] < t_amin * v.X[2 * i] + t_bmin)
									{
										if (!regression)
										{
											newO[n_classes + v.Y[i]] += 1;
										}
										else
										{
											newO[3] += v.Y_reg[i];
											newO[4] += v.Y_reg[i] * v.Y_reg[i];
											newO[5] += 1;
											if (env->criterion == Environment::mae)
											{
												new_reg2->push_back(v.Y_reg[i]);
											}
										}
										same_as_parent = false;
									}
									else
									{
										newX[2 * t_n_P] = v.X[2 * i];
										newX[2 * t_n_P + 1] = v.X[2 * i + 1];
										newY[t_n_P] = v.Y[i];
										if (!regression)
										{
											newY[t_n_P] = v.Y[i];
										}
										else
										{
											newY_reg[t_n_P] = v.Y_reg[i];
										}
										if (!regression)
										{
											newO[2 * n_classes + v.Y[i]] += 1;
										}
										else
										{
											newO[6] += v.Y_reg[i];
											newO[7] += v.Y_reg[i] * v.Y_reg[i];
											newO[8] += 1;
											if (env->criterion == Environment::mae)
											{
												new_reg3->push_back(v.Y_reg[i]);
											}
										}
										t_n_P += 1;
									}
								}
								else
								{
									if (v.X[1 * 2 * i + 1] >= t_amin * v.X[2 * i] + t_bmax)
									{
										if (!regression)
										{
											newO[v.Y[i]] += 1;
										}
										else
										{
											newO[0] += v.Y_reg[i];
											newO[1] += v.Y_reg[i] * v.Y_reg[i];
											newO[2] += 1;
											if (env->criterion == Environment::mae)
											{
												new_reg1->push_back(v.Y_reg[i]);
											}
										}
										same_as_parent = false;
									}
									else if (v.X[1 * 2 * i + 1] < t_amax * v.X[2 * i] + t_bmin)
									{
										if (!regression)
										{
											newO[n_classes + v.Y[i]] += 1;
										}
										else
										{
											newO[3] += v.Y_reg[i];
											newO[4] += v.Y_reg[i] * v.Y_reg[i];
											newO[5] += 1;
											if (env->criterion == Environment::mae)
											{
												new_reg2->push_back(v.Y_reg[i]);
											}
										}
										same_as_parent = false;
									}
									else
									{
										newX[2 * t_n_P] = v.X[2 * i];
										newX[2 * t_n_P + 1] = v.X[2 * i + 1];
										if (!regression)
										{
											newY[t_n_P] = v.Y[i];
										}
										else
										{
											newY_reg[t_n_P] = v.Y_reg[i];
										}
										if (!regression)
										{
											newO[2 * n_classes + v.Y[i]] += 1;
										}
										else
										{
											newO[6] += v.Y_reg[i];
											newO[7] += v.Y_reg[i] * v.Y_reg[i];
											newO[8] += 1;
											if (env->criterion == Environment::mae)
											{
												new_reg3->push_back(v.Y_reg[i]);
											}
										}
										t_n_P += 1;
									}
								}
							}
							int ind = 0;
							BBNode node = BBNode(newX, newY, newY_reg, t_amin, t_amax, t_bmin, t_bmax, newO, o, new_reg1, new_reg2, new_reg3, t_n_P, same_as_parent, v.LB, n_classes, total_points, dists, (v.type + 1) % 2, regression, env->criterion);
							while (ind < l && node.LB > nodes[ind].LB)
							{
								ind++;
							}
							for (int ind2 = l - 1; ind2 >= ind; ind2--)
							{
								nodes[ind2 + 1] = nodes[ind2];
							}
							nodes[ind] = node;
						}
						for (int l = 1; l >= 0; l--)
						{
							que.push_back(nodes[l]);
						}
					}
				}
			}
			else if (v.status == 0)
			{
				v.test_all();
				bool flag = false;
				mtx.lock();
				if (v.UB < res->bestImp)
				{
					bool verified = false;
					if (verification && !verify(total_points, regression, X, X2, X1, dataX1, dataX2, dataX1_reg, dataX2_reg, v.a, v.b, d_x, d_y, D_x, D_y))
					{
						v.upperBound_mid();
						if (v.UB < res->bestImp)
						{
							if (verify(total_points, regression, X, X2, X1, dataX1, dataX2, dataX1_reg, dataX2_reg, v.a, v.b, d_x, d_y, D_x, D_y))
							{
								verified = true;
							}
							else
							{
							}
						}
					}
					else
					{
						verified = true;
					}
					if (verified)
					{
						res->split_found = true;
						res->bestImp = v.UB;
						res->bestA = v.a * (D_y / D_x);
						res->bestB = D_y * v.b + d_y - v.a * d_x * (D_y / D_x);
						res->bestInd1 = ind1;
						res->bestInd2 = ind2;
						flag = true;
					}
				}
				mtx.unlock();
				if (flag)
				{
					bestA = v.a;
					bestB = v.b;
					bestImp = v.UB;
					herkunft = 2;
					update = true;
				}
			}
			else if (v.status == -1)
			{
				v.test_4();
				bool flag = false;
				mtx.lock();
				if (v.UB < res->bestImp)
				{
					bool verified = false;
					if (verification && !verify(total_points, regression, X, X2, X1, dataX1, dataX2, dataX1_reg, dataX2_reg, v.a, v.b, d_x, d_y, D_x, D_y))
					{
						v.upperBound_mid();
						if (v.UB < res->bestImp)
						{
							if (verify(total_points, regression, X, X2, X1, dataX1, dataX2, dataX1_reg, dataX2_reg, v.a, v.b, d_x, d_y, D_x, D_y))
							{
								verified = true;
							}
							else
							{
							}
						}
					}
					else
					{
						verified = true;
					}
					if (verified)
					{
						res->split_found = true;
						res->bestImp = v.UB;
						res->bestA = v.a * (D_y / D_x);
						res->bestB = D_y * v.b + d_y - v.a * d_x * (D_y / D_x);
						res->bestInd1 = ind1;
						res->bestInd2 = ind2;
						flag = true;
					}
				}
				mtx.unlock();
				if (flag)
				{
					bestA = v.a;
					bestB = v.b;
					bestImp = v.UB;
					herkunft = 3;
					update = true;
				}
			}
			else if (v.status == -2)
			{
			}
			else if (v.status == -3)
			{
				bool flag = false;
				mtx.lock();
				if (v.UB < res->bestImp)
				{
					bool verified = false;
					if (verification && !verify(total_points, regression, X, X2, X1, dataX1, dataX2, dataX1_reg, dataX2_reg, v.a, v.b, d_x, d_y, D_x, D_y))
					{
						v.upperBound_mid();
						if (v.UB < res->bestImp)
						{
							if (verify(total_points, regression, X, X2, X1, dataX1, dataX2, dataX1_reg, dataX2_reg, v.a, v.b, d_x, d_y, D_x, D_y))
							{
								verified = true;
							}
							else
							{
							}
						}
					}
					else
					{
						verified = true;
					}
					if (verified)
					{
						res->split_found = true;
						res->bestImp = v.UB;
						res->bestA = v.a * (D_y / D_x);
						res->bestB = D_y * v.b + d_y - v.a * d_x * (D_y / D_x);
						res->bestInd1 = ind1;
						res->bestInd2 = ind2;
						flag = true;
					}
				}
				mtx.unlock();
				if (flag)
				{
					bestA = v.a;
					bestB = v.b;
					bestImp = v.UB;
					herkunft = 4;
					update = true;
				}
			}
		}
		if (free_size >= max_free_size)
		{
			delete[] v.O;
			delete[] v.X;
			if (!regression)
			{
				delete[] v.Y;
			}
			else
			{
				delete[] v.Y_reg;
				if (env->criterion == Environment::mae)
				{
					delete v.reg1;
					delete v.reg2;
					delete v.reg3;
				}
			}
		}
		else
		{
			freeXs[free_size] = v.X;
			if (!regression)
			{
				freeYs[free_size] = v.Y;
			}
			else
			{
				freeYs_reg[free_size] = v.Y_reg;
				if (env->criterion == Environment::mae)
				{
					free_reg1[free_size] = v.reg1;
					free_reg2[free_size] = v.reg2;
					free_reg3[free_size] = v.reg3;
				}
			}
			freeOs[free_size] = v.O;
			free_size++;
		}
	}
	for (int i = 0; i < free_size; i++)
	{
		delete[] freeOs[i];
		delete[] freeXs[i];
		if (!regression)
		{
			delete[] freeYs[i];
		}
		else
		{
			delete[] freeYs_reg[i];
			if (env->criterion == Environment::mae)
			{
				delete free_reg1[i];
				delete free_reg2[i];
				delete free_reg3[i];
			}
		}
	}
	bool flag = false;
	mtx.lock();
	if (res->bestInd1 == ind1 && res->bestInd2 == ind2)
	{
		flag = true;
	}
	mtx.unlock();
	bool above = false;
	bool below = false;
	if (flag)
	{
		double tmp;
		for (int i = 0; i < total_points; i++)
		{
			tmp = X[total_points + i] - (bestA * X[i] + bestB);
			if (tmp >= 0)
			{
				above = true;
			}
			else
			{
				below = true;
			}
			if (above && below)
			{
				break;
			}
		}
	}
	mtx.lock();
	if (update && res->bestInd1 == ind1 && res->bestInd2 == ind2)
	{
		res->split_found = true;
		res->bestImp = bestImp;
		res->bestA = bestA * (D_y / D_x);
		res->bestB = D_y * bestB + d_y - bestA * d_x * (D_y / D_x);
		res->bestInd1 = ind1;
		res->bestInd2 = ind2;
		if (above && below)
		{
			res->divides = true;
		}
	}
	mtx.unlock();
	delete[] X;
	if (!regression)
	{
		delete[] Y;
		delete[] o;
		delete[] o1;
		delete[] o2;
	}
	else
	{
		delete[] Y_reg;
		delete[] o;
	}
	delete[] nodes;
	delete[] dists;
}
