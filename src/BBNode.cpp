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
#include <mutex>
#include "BBNode.h"
#include "Result.h"
#include "Environment.h"
#include <cmath>
#include <bits/stdc++.h>
using namespace std;
int compare(const void *a, const void *b)
{
	if (*(double *)a > *(double *)b)
		return 1;
	if (*(double *)a > *(double *)b)
		return -1;
	return 0;
}
BBNode::BBNode() {}
BBNode::BBNode(double *X, int *Y, double *Y_reg, double amin, double amax, double bmin, double bmax,
			   double *O, double *o, Vector_reg *reg1, Vector_reg *reg2, Vector_reg *reg3, int n_P, bool same_as_parent, double parentLB, int n_classes, int total_points, doubleIntPair *dists, int type, bool regression, int criterion)
{
	this->O = O;
	this->o = o;
	this->dists = dists;
	this->n_P = n_P;
	this->amin = amin;
	this->amax = amax;
	this->bmin = bmin;
	this->bmax = bmax;
	this->n_classes = n_classes;
	this->total_points = total_points;
	this->X = X;
	this->Y = Y;
	this->regression = regression;
	this->Y_reg = Y_reg;
	this->reg1 = reg1;
	this->reg2 = reg2;
	this->reg3 = reg3;
	this->status = 1;
	this->same_as_parent = same_as_parent;
	this->LB = parentLB;
	this->UB = std::numeric_limits<double>::max();
	this->type = type;
	this->criterion = criterion;
	eps = std::numeric_limits<double>::epsilon();
	eps = 1e-7;
	if (same_as_parent == false && !regression)
	{
		lowerBound2();
	}
	else if (same_as_parent == false)
	{
		lowerBound_reg();
	}
}
bool BBNode::singleClassTest()
{
	double t_b = 0;
	double t_a = 0;
	double t = 0;
	double a_u = amin;
	double b_u = bmax;
	double a_l = amin;
	double b_l = bmin;
	bool above = true;
	bool below = true;
	if (X[0] >= 0)
	{
		a_u = amax;
		a_l = amin;
		b_u = bmax;
		b_l = bmin;
	}
	else if (X[2 * (n_P - 1)] <= 0)
	{
		a_u = amin;
		a_l = amax;
		b_u = bmax;
		b_l = bmin;
	}
	else
	{
		double ub_u = amax;
		double lb_u = amin;
		double ub_l = amax;
		double lb_l = amin;
		above = false;
		below = false;
		for (int i = 0; i < n_P; i++)
		{
			if (X[2 * i] > 0)
			{
				t = (X[2 * i + 1] - bmax) / X[2 * i];
				if (t > lb_u)
				{
					lb_u = t;
				}
				t = (X[2 * i + 1] - bmin) / X[2 * i];
				if (t < ub_l)
				{
					ub_l = t;
				}
			}
			if (X[2 * i] < 0)
			{
				t = (X[2 * i + 1] - bmax) / X[2 * i];
				if (t < ub_u)
				{
					ub_u = t;
				}
				t = (X[2 * i + 1] - bmin) / X[2 * i];
				if (t > lb_l)
				{
					lb_l = t;
				}
			}
			if (lb_u >= ub_u)
			{
				return false;
			}
			if (lb_l > ub_l)
			{
				return false;
			}
		}
		if (lb_u < ub_u - eps)
		{
			a_u = (lb_u + ub_u) / 2;
			b_u = bmax;
			above = true;
		}
		if (lb_l <= ub_l)
		{
			a_l = (lb_l + ub_l) / 2;
			b_l = bmin;
			below = true;
		}
	}
	if (above && below)
	{
		a = a_u;
		b = b_u;
		for (int m = 0; m < n_classes; m++)
		{
			o[m] = O[m];
			o[n_classes + m] = O[n_classes + m];
			o[n_classes + m] += O[2 * n_classes + m];
		}
		double t_imp = imp(o, n_classes, criterion);
		if (t_imp < UB)
		{
			a = a_u;
			b = b_u;
			UB = t_imp;
		}
		for (int m = 0; m < n_classes; m++)
		{
			o[m] = O[m];
			o[n_classes + m] = O[n_classes + m];
			o[m] += O[2 * n_classes + m];
		}
		t_imp = imp(o, n_classes, criterion);
		if (t_imp < UB)
		{
			a = a_l;
			b = b_l;
			UB = t_imp;
		}
		return true;
	}
	return false;
}
bool BBNode::computeBranchingRule()
{
	branch_a = (amin + amax) / 2;
	branch_b = (bmin + bmax) / 2;
	branch_a = 0;
	branch_b = 0;
	double b1 = 0;
	double b2 = 0;
	for (int i = 0; i < n_P; i++)
	{
		b1 = X[2 * i + 1] - (amin * X[2 * i]);
		b2 = X[2 * i + 1] - (amax * X[2 * i]);
		if (b1 > bmax)
		{
			b1 = bmax;
		}
		if (b2 > bmax)
		{
			b2 = bmax;
		}
		if (b1 < bmin)
		{
			b1 = bmin;
		}
		if (b2 < bmin)
		{
			b2 = bmin;
		}
		double mid_b = (b1 + b2) / 2.0;
		double a_2 = 0;
		double b_2 = mid_b;
		if (X[2 * i] == 0)
		{
			a_2 = (amin + amax) / 2.0;
			b_2 = X[2 * i + 1];
		}
		else
		{
			a_2 = (X[2 * i + 1] - mid_b) / (X[2 * i]);
		}
		if (a_2 < amin + eps)
		{
			a_2 = (amin + amax) / 2.0;
		}
		if (a_2 > amax - eps)
		{
			a_2 = (amin + amax) / 2.0;
		}
		branch_a += a_2;
		branch_b += b_2;
	}
	branch_a /= n_P;
	branch_b /= n_P;
	return false;
}
int comparePairs(const void *a, const void *b)
{
	doubleIntPair *pairA = (doubleIntPair *)a;
	doubleIntPair *pairB = (doubleIntPair *)b;
	if (pairA->dist < pairB->dist)
	{
		return -1;
	}
	if (pairA->dist > pairB->dist)
	{
		return 1;
	}
	else
	{
		return 0;
	}
	return (pairB->dist - pairA->dist);
}
bool comparePairs2(doubleIntPair a, doubleIntPair b)
{
	return (b.dist > a.dist);
}
bool BBNode::computeBranchingRule_B(double *o0, double *o1)
{
	for (int m = 0; m < n_classes; m++)
	{
		o0[m] = 0;
		o0[n_classes + m] = 0;
		o0[2 * n_classes + m] = 0;
		o1[m] = 0;
		o1[n_classes + m] = 0;
		o1[2 * n_classes + m] = 0;
	}
	branch_a = (amin + amax) / 2.0;
	branch_b = 0;
	int k = 0;
	double tmp0;
	int n_1 = n_P;
	int n_0 = 0;
	double tmp1;
	for (int i = 0; i < n_P; i++)
	{
		o1[Y[i]] += 1;
		if (X[2 * i] > 0)
		{
			tmp1 = X[2 * i + 1] - (amin * X[2 * i] + bmin);
			tmp0 = X[2 * i + 1] - (amax * X[2 * i] + bmin);
		}
		else if (X[2 * i] < 0)
		{
			tmp1 = X[2 * i + 1] - (amax * X[2 * i] + bmin);
			tmp0 = X[2 * i + 1] - (amin * X[2 * i] + bmin);
		}
		else
		{
			tmp1 = X[2 * i + 1] - bmin;
			tmp0 = X[2 * i + 1] - bmin;
		}
		if (tmp0 <= eps)
		{
			n_0 += 1;
			if (!regression)
			{
				o0[Y[i]] += 1;
			}
			else
			{
				o0[0] += Y_reg[i];
				o0[1] += Y_reg[i] * Y_reg[i];
				o0[2] += 1;
			}
		}
		if (tmp1 <= eps)
		{
			n_1 -= 1;
			if (!regression)
			{
				o1[Y[i]] -= 1;
			}
			else
			{
				o1[0] -= Y_reg[i];
				o1[1] -= Y_reg[i] * Y_reg[i];
				o1[2] -= 1;
			}
		}
		if (tmp1 > eps && bmin + tmp1 < bmax)
		{
			dists[k].dist = tmp1;
			dists[k].sect = 1;
			dists[k].y = Y[i];
			k++;
		}
		if (tmp0 > eps && bmin + tmp0 < bmax)
		{
			dists[k].dist = tmp0;
			dists[k].sect = 0;
			dists[k].y = Y[i];
			k++;
		}
	}
	sort(dists, dists + k, comparePairs2);
	double t_max;
	int t_max2;
	double max_dist = (bmin + bmax) / 2.0 - bmin;
	double max = std::numeric_limits<double>::max();
	int sum = 2 * n_P;
	int between_0 = 0;
	int between_1 = 0;
	double dist_before = eps;
	double t_dist;
	int maxit = 0;
	int max_n_0 = n_0;
	int max_n_1 = n_1;
	int start = 0;
	for (int i = start; i < k; i++)
	{
		t_dist = (dists[i].dist + dist_before) / 2.0;
		t_max2 = n_1;
		if (n_0 > n_1)
		{
			t_max2 = n_0;
		}
		t_max = t_max2;
		if ((t_max <= max && i > 0) || (true && t_max <= max && i == 0 && t_max2 < n_P))
		{
			if (t_max == max)
			{
				if (n_0 + n_1 < sum)
				{
					max = t_max;
					max_dist = t_dist;
					sum = n_0 + n_1;
					maxit = i;
					max_n_0 = n_0;
					max_n_1 = n_1;
				}
			}
			else
			{
				max = t_max;
				max_dist = t_dist;
				sum = n_0 + n_1;
				maxit = i;
				max_n_0 = n_0;
				max_n_1 = n_1;
			}
		}
		between_0 = 0;
		between_1 = 0;
		if (dists[i].sect == 0)
		{
			n_0++;
			o0[dists[i].y] += 1;
		}
		if (dists[i].sect == 1)
		{
			n_1--;
			o1[dists[i].y] -= 1;
		}
		while (i + 1 < k && dists[i].dist == dists[i + 1].dist)
		{
			i++;
			if (dists[i].sect == 0)
			{
				n_0++;
				o0[dists[i].y] += 1;
			}
			if (dists[i].sect == 1)
			{
				n_1--;
				o1[dists[i].y] -= 1;
			}
		}
		dist_before = dists[i].dist;
	}
	if (max_n_0 == n_P || max_n_1 == n_P)
	{
		max_dist = (bmin + bmax) / 2.0 - bmin;
	}
	branch_b = bmin + max_dist;
	branch_a = 0;
	for (int i = 0; i < n_P; i++)
	{
		if (X[2 * i] != 0.0)
		{
			tmp0 = (X[2 * i + 1] - branch_b) / X[2 * i];
		}
		else
		{
			tmp0 = (amin + amax) / 2.0;
		}
		if (tmp0 < amin)
		{
			tmp0 = amin;
		}
		if (tmp0 > amax)
		{
			tmp0 = amax;
		}
		branch_a += tmp0;
	}
	branch_a /= n_P;
	return max == n_P;
}
bool BBNode::computeBranchingRule_A(double *o1, double *o2)
{
	branch_b = (bmin + bmax) / 2.0;
	branch_a = 0;
	branch_b = 0;
	int k = 0;
	double tmp;
	int n_1 = n_P;
	int n_0 = 0;
	double tmp2;
	int between_0 = 0;
	int between_1 = 0;
	for (int i = 0; i < n_P; i++)
	{
		if (X[2 * i] > 0)
		{
			tmp = (X[2 * i + 1] - bmin) / (X[2 * i]) - amin;
			tmp2 = (X[2 * i + 1] - bmax) / (X[2 * i]) - amin;
		}
		else if (X[2 * i] < 0)
		{
			tmp = (X[2 * i + 1] - bmax) / (X[2 * i]) - amin;
			tmp2 = (X[2 * i + 1] - bmin) / (X[2 * i]) - amin;
		}
		if (tmp2 <= eps)
		{
			n_0 += 1;
		}
		if (tmp <= eps)
		{
			n_1 -= 1;
		}
		if (tmp > eps && amax - (amin + tmp) > eps)
		{
			dists[k].dist = tmp;
			dists[k].sect = 1;
			k++;
		}
		if (tmp2 > 0 && amax - (amin + tmp2) > eps)
		{
			dists[k].dist = tmp2;
			dists[k].sect = 0;
			k++;
		}
	}
	sort(dists, dists + k, comparePairs2);
	double t_max;
	int t_max2;
	double max_dist = (amin + amax) / 2.0 - amin;
	double max = std::numeric_limits<double>::max();
	int sum = 2 * n_P;
	int max_n_0 = n_0;
	int max_n_1 = n_1;
	double dist_before = eps;
	double t_dist;
	int start = 0;
	int maxit = 0;
	for (int i = start; i < k; i++)
	{
		t_dist = (dists[i].dist + dist_before) / 2.0;
		t_max2 = n_1;
		if (n_0 > n_1)
		{
			t_max2 = n_0;
		}
		t_max = t_max2;
		if ((t_max <= max && i > 0) || (true && t_max <= max && i == 0 && t_max2 < n_P))
		{
			if (t_max == max)
			{
				if (n_0 + n_1 < sum)
				{
					max = t_max;
					max_dist = t_dist;
					sum = n_0 + n_1;
					max_n_0 = n_0;
					max_n_1 = n_1;
				}
			}
			else
			{
				max = t_max;
				max_dist = t_dist;
				sum = n_0 + n_1;
				max_n_0 = n_0;
				max_n_1 = n_1;
			}
		}
		if (dists[i].sect == 0)
		{
			n_0++;
		}
		if (dists[i].sect == 1)
		{
			between_1++;
			n_1--;
		}
		while (i + 1 < k && dists[i].dist == dists[i + 1].dist)
		{
			i++;
			if (dists[i].sect == 0)
			{
				n_0++;
			}
			if (dists[i].sect == 1)
			{
				n_1--;
			}
		}
		dist_before = dists[i].dist;
	}
	if (max == n_P)
	{
		max_dist = (amin + amax) / 2.0 - amin;
	}
	branch_a = amin + max_dist;
	branch_b = 0;
	for (int i = 0; i < n_P; i++)
	{
		tmp = X[2 * i + 1] - branch_a * X[2 * i];
		if (tmp < bmin)
		{
			tmp = bmin;
		}
		if (tmp > bmax)
		{
			tmp = bmax;
		}
		branch_b += tmp;
	}
	branch_b /= n_P;
	return max == n_P;
}
void BBNode::computeBranchingRule_B_Bound(bool switched)
{
	if (switched && bmax - bmin <= eps)
	{
		computeBranchingRule();
		upperBound();
		return;
	}
	if (false && bmax - bmin <= eps)
	{
		type = 1;
		computeBranchingRule_A_Bound(true);
		return;
	}
	for (int m = 0; m < n_classes; m++)
	{
		o[m] = O[m];
		o[m] += O[2 * n_classes + m];
		o[n_classes + m] = O[n_classes + m];
	}
	computeBranchingRule();
	a = branch_a;
	b = bmin;
	int k = 0;
	bool flag_below = false;
	bool flag_above = false;
	double tmp;
	for (int i = 0; i < n_P; i++)
	{
		tmp = X[2 * i + 1] - (a * X[2 * i] + bmin);
		if (tmp > eps && bmin + tmp < bmax)
		{
			dists[k].dist = tmp;
			dists[k].sect = 1;
			dists[k].y = Y[i];
			k++;
		}
		else if (tmp <= eps)
		{
			o[Y[i]]--;
			o[n_classes + Y[i]]++;
			flag_below = true;
		}
		else
		{
			flag_above = true;
		}
	}
	if ((k == 0 || dists[0].dist == dists[k - 1].dist) && switched)
	{
		computeBranchingRule();
		upperBound();
		return;
	}
	if (false && (k == 0 || dists[0].dist == dists[k - 1].dist))
	{
		type = 1;
		computeBranchingRule_A_Bound(true);
		return;
	}
	sort(dists, dists + k, comparePairs2);
	bool flag = false;
	int i = 0;
	double delta = 0;
	int n1 = 0;
	int n2 = 0;
	for (int j = 0; j < n_classes; j++)
	{
		n1 += o[j];
		n2 += o[n_classes + j];
	}
	while (i < k)
	{
		o[dists[i].y]--;
		o[n_classes + dists[i].y]++;
		n1--;
		n2++;
		if (i == k - 1)
		{
			break;
		}
		if (dists[i].dist < dists[i + 1].dist - eps)
		{
			delta = (dists[i].dist + dists[i + 1].dist) / 2;
			tmp = imp(o, n_classes, criterion);
			if (tmp < UB && n1 > 0 && n2 > 0)
			{
				UB = tmp;
				b = bmin + delta;
				flag = true;
			}
		}
		i++;
	}
	if (!flag)
	{
		upperBound();
	}
}
void BBNode::computeBranchingRule_A_Bound(bool switched)
{
	if (switched && amax - amin <= eps)
	{
		computeBranchingRule();
		upperBound();
		return;
	}
	if (false && amax - amin <= eps)
	{
		type = 0;
		computeBranchingRule_B_Bound(true);
		return;
	}
	for (int m = 0; m < n_classes; m++)
	{
		o[m] = O[m];
		o[m] += O[2 * n_classes + m];
		o[n_classes + m] = O[n_classes + m];
	}
	computeBranchingRule();
	b = branch_b;
	a = amin;
	int k = 0;
	bool flag_below = false;
	bool flag_above = false;
	double tmp;
	for (int i = 0; i < n_P; i++)
	{
		if (X[2 * i] == 0)
		{
			if (X[2 * i + 1] < branch_b)
			{
				o[Y[i]]--;
				o[n_classes + Y[i]]++;
			}
			continue;
		}
		tmp = (X[2 * i + 1] - b) / X[2 * i] - amin;
		if (tmp > eps && amin + tmp < amax)
		{
			dists[k].dist = tmp;
			dists[k].sect = 1;
			if (X[2 * i] < 0)
			{
				dists[k].sect = 0;
				o[Y[i]]--;
				o[n_classes + Y[i]]++;
			}
			dists[k].y = Y[i];
			k++;
		}
		else if (tmp <= eps)
		{
			if (X[2 * i] > 0)
			{
				o[Y[i]]--;
				o[n_classes + Y[i]]++;
				flag_below = true;
			}
		}
		else
		{
			if (X[2 * i] < 0)
			{
				o[Y[i]]--;
				o[n_classes + Y[i]]++;
				flag_below = true;
			}
		}
	}
	if ((k == 0 || dists[0].dist == dists[k - 1].dist) && switched)
	{
		computeBranchingRule();
		upperBound();
		return;
	}
	if (false && (k == 0 || dists[0].dist == dists[k - 1].dist))
	{
		type = 0;
		computeBranchingRule_B_Bound(true);
		return;
	}
	sort(dists, dists + k, comparePairs2);
	bool flag = false;
	int i = 0;
	double delta = 0;
	int n1 = 0;
	int n2 = 0;
	for (int j = 0; j < n_classes; j++)
	{
		n1 += o[j];
		n2 += o[n_classes + j];
	}
	while (i < k)
	{
		if (dists[i].sect == 1)
		{
			o[dists[i].y]--;
			o[n_classes + dists[i].y]++;
			n1--;
			n2++;
		}
		else
		{
			o[dists[i].y]++;
			o[n_classes + dists[i].y]--;
			n1++;
			n2--;
		}
		if (i == k - 1)
		{
			break;
		}
		if (dists[i].dist < dists[i + 1].dist - eps)
		{
			delta = (dists[i].dist + dists[i + 1].dist) / 2;
			tmp = imp(o, n_classes, criterion);
			if (tmp < UB && n1 > 0 && n2 > 0)
			{
				UB = tmp;
				a = amin + delta;
				flag = true;
			}
		}
		i++;
	}
	if (!flag)
	{
		upperBound();
	}
}
void BBNode::upperBound()
{
	if (!regression)
	{
		for (int m = 0; m < n_classes; m++)
		{
			o[m] = O[m];
			o[n_classes + m] = O[n_classes + m];
		}
	}
	else
	{
		if (criterion == Environment::mse)
		{
			for (int m = 0; m < 6; m++)
			{
				o[m] = O[m];
			}
		}
		else
		{
			reg1_tmp = new Vector_reg();
			reg2_tmp = new Vector_reg();
			reg1_tmp->reserve(total_points);
			reg2_tmp->reserve(total_points);
			reg1_tmp->use_mae();
			reg2_tmp->use_mae();
			reg1_tmp->add(reg1);
			reg2_tmp->add(reg2);
		}
	}
	double ub_a = branch_a;
	double ub_b = branch_b;
	for (int i = 0; i < n_P; i++)
	{
		if (abs(X[2 * i + 1] - ub_a * X[2 * i] - ub_b) < 1e-10)
		{
			ub_a = amin + (double)(rand()) / ((double)(RAND_MAX / (amax - amin)));
			ub_b = bmin + (double)(rand()) / ((double)(RAND_MAX / (bmax - bmin)));
		}
	}
	for (int i = 0; i < n_P; i++)
	{
		if (X[2 * i + 1] >= ub_a * X[2 * i] + ub_b)
		{
			if (!regression)
			{
				o[Y[i]] += 1;
			}
			else
			{
				if (criterion == Environment::mse)
				{
					o[0] += Y_reg[i];
					o[1] += Y_reg[i] * Y_reg[i];
					o[2] += 1;
				}
				else
				{
					reg1_tmp->insert_inplace(Y_reg[i]);
				}
			}
		}
		else
		{
			if (!regression)
			{
				o[n_classes + Y[i]] += 1;
			}
			else
			{
				if (criterion == Environment::mse)
				{
					o[3] += Y_reg[i];
					o[4] += Y_reg[i] * Y_reg[i];
					o[5] += 1;
				}
				else
				{
					reg2_tmp->insert_inplace(Y_reg[i]);
				}
			}
		}
	}
	double tmp;
	bool splits = false;
	if (!regression)
	{
		int n1 = 0;
		int n2 = 0;
		for (int k = 0; k < n_classes; k++)
		{
			n1 += o[k];
			n2 += o[n_classes + k];
		}
		if (n1 > 0 && n2 > 0)
		{
			splits = true;
			tmp = imp(o, n_classes, criterion);
		}
	}
	else
	{
		if (criterion == Environment::mse)
		{
			if (o[2] > 0 && o[5] > 0)
			{
				tmp = (o[1] - 2 * o[0] * (o[0] / o[2]) + o[2] * (o[0] / o[2]) * (o[0] / o[2]));
				tmp += (o[4] - 2 * o[3] * (o[3] / o[5]) + o[5] * (o[3] / o[5]) * (o[3] / o[5]));
				splits = true;
			}
			else
			{
				tmp = std::numeric_limits<double>::max();
				splits = false;
			}
		}
		else
		{
			if (reg1_tmp->size() > 0 && reg2_tmp->size() > 0)
			{
				tmp = reg1_tmp->error(criterion);
				tmp += reg2_tmp->error(criterion);
				splits = true;
			}
			else
			{
				tmp = std::numeric_limits<double>::max();
				splits = false;
			}
		}
	}
	if (tmp < UB && splits)
	{
		UB = tmp;
		a = ub_a;
		b = ub_b;
	}
	if (regression && criterion != Environment::mse)
	{
		delete reg1_tmp;
		delete reg2_tmp;
	}
}
void BBNode::upperBound_mid()
{
	UB = std::numeric_limits<double>::max();
	double a_tmp = (amin + amax) / 2;
	double b_tmp = (bmin + bmax) / 2;
	a_tmp = amin + (double)(rand()) / ((double)(RAND_MAX / (amax - amin)));
	b_tmp = bmin + (double)(rand()) / ((double)(RAND_MAX / (bmax - bmin)));
	if (!regression)
	{
		for (int m = 0; m < n_classes; m++)
		{
			o[m] = O[m];
			o[n_classes + m] = O[n_classes + m];
		}
	}
	else
	{
		if (criterion == Environment::mse)
		{
			for (int m = 0; m < 6; m++)
			{
				o[m] = O[m];
			}
		}
		else
		{
			reg1_tmp = new Vector_reg();
			reg2_tmp = new Vector_reg();
			reg1_tmp->reserve(total_points);
			reg2_tmp->reserve(total_points);
			reg1_tmp->use_mae();
			reg2_tmp->use_mae();
			reg1_tmp->add(reg1);
			reg2_tmp->add(reg2);
		}
	}
	for (int i = 0; i < n_P; i++)
	{
		if (X[2 * i + 1] >= a_tmp * X[2 * i] + b_tmp)
		{
			if (!regression)
			{
				o[Y[i]] += 1;
			}
			else
			{
				if (criterion == Environment::mse)
				{
					o[0] += Y_reg[i];
					o[1] += Y_reg[i] * Y_reg[i];
					o[2] += 1;
				}
				else
				{
					reg1_tmp->insert_inplace(Y_reg[i]);
				}
			}
		}
		else
		{
			if (!regression)
			{
				o[n_classes + Y[i]] += 1;
			}
			else
			{
				if (criterion == Environment::mse)
				{
					o[3] += Y_reg[i];
					o[4] += Y_reg[i] * Y_reg[i];
					o[5] += 1;
				}
				else
				{
					reg2_tmp->insert_inplace(Y_reg[i]);
				}
			}
		}
	}
	double tmp;
	bool splits = false;
	if (!regression)
	{
		int n1 = 0;
		int n2 = 0;
		for (int k = 0; k < n_classes; k++)
		{
			n1 += o[k];
			n2 += o[n_classes + k];
		}
		if (n1 > 0 && n2 > 0)
		{
			splits = true;
			tmp = imp(o, n_classes, criterion);
		}
	}
	else
	{
		if (criterion == Environment::mse)
		{
			if (o[2] > 0 && o[5] > 0)
			{
				tmp = (o[1] - 2 * o[0] * (o[0] / o[2]) + o[2] * (o[0] / o[2]) * (o[0] / o[2]));
				tmp += (o[4] - 2 * o[3] * (o[3] / o[5]) + o[5] * (o[3] / o[5]) * (o[3] / o[5]));
				splits = true;
			}
			else
			{
				tmp = std::numeric_limits<double>::max();
				splits = false;
			}
		}
		else
		{
			if (reg1_tmp->size() > 0 && reg2_tmp->size() > 0)
			{
				tmp = reg1_tmp->error(criterion);
				tmp += reg2_tmp->error(criterion);
				splits = true;
			}
			else
			{
				tmp = std::numeric_limits<double>::max();
				splits = false;
			}
		}
	}
	if (tmp < UB && splits)
	{
		UB = tmp;
		a = a_tmp;
		b = b_tmp;
	}
}
void BBNode::upperBound(double a_t, double b_t)
{
	for (int m = 0; m < n_classes; m++)
	{
		o[m] = O[m];
		o[n_classes + m] = O[n_classes + m];
	}
	for (int i = 0; i < n_P; i++)
	{
		if (X[2 * i + 1] >= a_t * X[2 * i] + b_t)
		{
			o[Y[i]] += 1;
		}
		else
		{
			o[n_classes + Y[i]] += 1;
		}
	}
	int n1 = 0;
	int n2 = 0;
	for (int k = 0; k < n_classes; k++)
	{
		n1 += o[k];
		n2 += o[n_classes + k];
	}
	if (n1 == 0 || n2 == 0)
	{
		return;
	}
	double tmp = imp(o, n_classes, criterion);
	if (tmp < UB)
	{
		UB = tmp;
		this->a = a_t;
		this->b = b_t;
	}
}
void BBNode::find_Bound(double *o1)
{
	for (int m = 0; m < n_classes; m++)
	{
		o1[m] = bmax;
		o1[n_classes + m] = bmin;
	}
	double tmp;
	for (int i = 0; i < n_P; i++)
	{
		tmp = X[2 * i + 1] - branch_a * X[2 * i];
		if (tmp < o1[Y[i]] && tmp >= bmin)
		{
			o1[Y[i]] = tmp;
		}
		if (tmp > o1[n_classes + Y[i]] && tmp <= bmax)
		{
			o1[n_classes + Y[i]] = tmp;
		}
	}
	for (int m = 0; m < n_classes; m++)
	{
		o1[m] = o1[m] + (bmax - o1[m]) / 1024;
	}
	for (int m = 0; m < n_classes; m++)
	{
		if (O[2 * n_classes + m] != 0)
		{
			upperBound(branch_a, o1[m]);
			upperBound(branch_a, o1[n_classes + m]);
		}
	}
}
void BBNode::simpleBound()
{
	LB = imp(O, n_classes, criterion);
}
void BBNode::lowerBound2()
{
	LB = imp(O, n_classes, criterion);
	if (criterion == Environment::gini_imp)
	{
		int maxInd = 0;
		int max = 0;
		for (int m = 0; m < n_classes; m++)
		{
			o[m] = O[2 * n_classes + m];
			o[n_classes + m] = 0;
		}
		for (int m = 0; m < n_classes; m++)
		{
			if (O[2 * n_classes + m] > max)
			{
				max = O[2 * n_classes + m];
				maxInd = m;
			}
		}
		o[maxInd] = 0;
		o[n_classes + maxInd] = max;
		LB += imp(o, n_classes, criterion);
	}
}
void BBNode::lowerBound3()
{
	LB = imp(O, n_classes, criterion);
	if (criterion == Environment::gini_imp)
	{
		int maxInd = 0;
		int max = 0;
		int maxInd2 = 0;
		int max2 = 0;
		for (int m = 0; m < n_classes; m++)
		{
			if (O[2 * n_classes + m] > max)
			{
				max = O[2 * n_classes + m];
				maxInd = m;
			}
		}
		for (int m = 0; m < n_classes; m++)
		{
			if (O[2 * n_classes + m] > max2 && m != maxInd)
			{
				max2 = O[2 * n_classes + m];
				maxInd2 = m;
			}
		}
		for (int m = 0; m < n_classes; m++)
		{
			o[m] = O[2 * n_classes + m];
			o[n_classes + m] = 0;
		}
		o[maxInd] = 0;
		o[maxInd2] = 0;
		o[n_classes + maxInd2] = max2;
		LB = imp(o, n_classes, criterion);
		double tmp;
		for (int i = 0; i < 2; i++)
		{
			for (int m = 0; m < n_classes; m++)
			{
				o[m] = O[m];
				o[n_classes + m] = O[n_classes + m];
			}
			if (i == 0)
			{
				o[maxInd] += max;
			}
			else
			{
				o[n_classes + maxInd] += max;
			}
			if (i == 0)
			{
				tmp = imp(o, n_classes, criterion);
			}
			else
			{
				double tmp2 = imp(o, n_classes, criterion);
				if (tmp2 < tmp)
				{
					tmp = tmp2;
				}
			}
		}
		LB += tmp;
	}
}
void BBNode::lowerBound_reg()
{
	LB = 0;
	if (criterion == Environment::mse)
	{
		if (O[2] > 0)
		{
			LB = (O[1] - 2 * O[0] * (O[0] / O[2]) + O[2] * (O[0] / O[2]) * (O[0] / O[2]));
		}
		if (O[5] > 0)
		{
			LB += (O[4] - 2 * O[3] * (O[3] / O[5]) + O[5] * (O[3] / O[5]) * (O[3] / O[5]));
		}
	}
	else
	{
		LB = reg1->error(Environment::mae) + reg2->error(Environment::mae);
	}
}
void BBNode::lowerBound()
{
	int N1 = 0;
	int N2 = 0;
	int maxInd1 = 0;
	int maxInd2 = 0;
	int maxN1 = 0;
	int maxN2 = 1;
	for (int m = 0; m < n_classes; m++)
	{
		o[m] = O[m];
		o[n_classes + m] = O[n_classes + m];
	}
	for (int m = 0; m < n_classes; m++)
	{
		N1 += O[m];
		N2 += O[n_classes + m];
		if (O[m] > maxN1)
		{
			maxN1 = O[m];
			maxInd1 = m;
		}
		if (O[n_classes + m] > N2)
		{
			maxN2 = O[n_classes + m];
			maxInd2 = m;
		}
	}
	bool all_in_one = false;
	int n = n_P;
	double t_imp = std::numeric_limits<double>::max();
	double tmp = t_imp;
	double best = UB;
	for (int i = 0; i <= n; i++)
	{
		o[maxInd1] = O[maxInd1] + i;
		o[n_classes + maxInd2] = O[n_classes + maxInd2] + (n - i);
		tmp = imp(o, n_classes, criterion);
		if (tmp < t_imp)
		{
			t_imp = tmp;
		}
	}
	LB = t_imp;
}
void BBNode::reduce()
{
	double t_bmin = std::numeric_limits<double>::max();
	double t_bmax = std::numeric_limits<double>::min();
	for (int i = 0; i < n_P; i++)
	{
		double tmp = X[2 * i + 1] - amin * X[2 * i];
		if (tmp > t_bmax)
		{
			t_bmax = tmp;
		}
		if (tmp < t_bmin)
		{
			t_bmin = tmp;
		}
		tmp = X[2 * i + 1] - amax * X[2 * i];
		if (tmp > t_bmax)
		{
			t_bmax = tmp;
		}
		if (tmp < t_bmin)
		{
			t_bmin = tmp;
		}
	}
	t_bmin = t_bmin;
	double e = (bmax - t_bmax) / 1024;
	if (eps < e / 10)
	{
		e = eps / 10;
	}
	t_bmax = t_bmax + e;
	if (t_bmax < bmax)
	{
		bmax = t_bmax;
	}
	if (t_bmin > bmin)
	{
		bmin = t_bmin;
	}
}
void BBNode::process()
{
	reduce();
	bool two_different = false;
	bool aligned = true;
	if (n_P == 0)
	{
		status = -2;
		return;
	}
	if (n_P == 1)
	{
		status = -1;
		return;
	}
	else
	{
		double p1_x = X[0];
		double p1_y = X[1];
		double p2_x = X[0];
		double p2_y = X[1];
		bool found_line = false;
		double tmp_a = 1;
		double tmp_b = 0;
		double tmp_a2 = 1;
		double tmp_b2 = 0;
		for (int i = 1; i < n_P; i++)
		{
			p2_x = X[2 * i];
			p2_y = X[2 * i + 1];
			if (abs(p2_x - p1_x) > eps || abs(p2_y - p1_y) > eps)
			{
				two_different = true;
				break;
			}
		}
		if (!two_different)
		{
			status = -1;
			return;
		}
		else
		{
			if (!regression)
			{
				int n_innerClasses = 0;
				for (int m = 0; m < n_classes; m++)
				{
					if (O[2 * n_classes + m] > 0)
					{
						n_innerClasses += 1;
					}
				}
				if (n_innerClasses == 1)
				{
					bool flag = singleClassTest();
					if (flag)
					{
						status = -3;
						return;
					}
				}
			}
			for (int i = 1; i < n_P; i++)
			{
				p2_x = X[2 * i];
				p2_y = X[2 * i + 1];
				if (abs(p2_x - p1_x) > eps || abs(p2_y - p1_y) > eps)
				{
					if (abs(p2_x - p1_x) <= eps)
					{
						aligned = false;
						break;
					}
					if (!found_line)
					{
						found_line = true;
						tmp_a = (p2_y - p1_y) / (p2_x - p1_x);
						tmp_b = p1_y - tmp_a * p1_x;
					}
					else
					{
						tmp_a2 = (p2_y - p1_y) / (p2_x - p1_x);
						tmp_b2 = p1_y - tmp_a2 * p1_x;
						if (abs(tmp_a2 - tmp_a) > eps || abs(tmp_b2 - tmp_b) > eps)
						{
							aligned = false;
							break;
						}
					}
				}
			}
		}
		if (!aligned)
		{
			status = 1;
			return;
		}
		else
		{
			status = 0;
		}
	}
}
void BBNode::test(double t_a, double t_b)
{
	if (!regression)
	{
		for (int m = 0; m < n_classes; m++)
		{
			o[m] = O[m];
			o[n_classes + m] = O[n_classes + m];
		}
	}
	else
	{
		for (int m = 0; m < 6; m++)
		{
			o[m] = O[m];
		}
	}
	for (int i = 0; i < n_P; i++)
	{
		if (X[2 * i + 1] >= t_a * X[2 * i] + t_b)
		{
			if (!regression)
			{
				o[Y[i]] += 1;
			}
			else
			{
				o[0] += Y_reg[i];
				o[1] += Y_reg[i] * Y_reg[i];
				o[2] += 1;
			}
		}
		else
		{
			if (!regression)
			{
				o[n_classes + Y[i]] += 1;
			}
			else
			{
				o[3] += Y_reg[i];
				o[4] += Y_reg[i] * Y_reg[i];
				o[5] += 1;
			}
		}
	}
	double t_imp;
	bool splits = false;
	if (!regression)
	{
		int n1 = 0;
		int n2 = 0;
		for (int k = 0; k < n_classes; k++)
		{
			n1 += o[k];
			n2 += o[n_classes + k];
		}
		if (n1 > 0 && n2 > 0)
		{
			splits = true;
			t_imp = imp(o, n_classes, criterion);
		}
		else
		{
			t_imp = std::numeric_limits<double>::max();
			splits = false;
		}
	}
	else
	{
		if (o[2] > 0)
		{
			t_imp = (o[1] - 2 * o[0] * (o[0] / o[2]) + o[2] * (o[0] / o[2]) * (o[0] / o[2]));
			if (o[5] > 0)
			{
				t_imp += (o[4] - 2 * o[3] * (o[3] / o[5]) + o[5] * (o[3] / o[5]) * (o[3] / o[5]));
				splits = true;
			}
		}
		else if (o[5] > 0)
		{
			t_imp = (o[4] - 2 * o[3] * (o[3] / o[5]) + o[5] * (o[3] / o[5]) * (o[3] / o[5]));
			splits = true;
		}
		else
		{
			t_imp = std::numeric_limits<double>::max();
			splits = false;
		}
	}
	if (t_imp < UB && splits)
	{
		UB = t_imp;
		a = t_a;
		b = t_b;
	}
}
void BBNode::test_all()
{
	double t_a = 0;
	double t_b = 0;
	double p1_x = X[0];
	double p1_y = X[1];
	double p2_x = X[0];
	double p2_y = X[1];
	double m_x = 0;
	double m_y = 0;
	double t = 0;
	double a_u = amin;
	double b_u = bmax;
	double a_l = amin;
	double b_l = bmin;
	for (int i = 0; i < n_P; i++)
	{
		if (X[2 * i] > 0.0)
		{
			t_a = (X[2 * i + 1] - bmax) / X[2 * i];
			t = 1.0;
			while (t_a + eps / t > amax && t <= 128.0)
			{
				t *= 2.0;
			}
			t_a = t_a + eps / t;
			if (t_a > a_u)
			{
				a_u = t_a;
			}
		}
		if (X[2 * i] < 0.0)
		{
			t_a = (X[2 * i + 1] - bmin) / X[2 * i];
			if (t_a > a_l)
			{
				a_l = t_a;
			}
		}
	}
	bool above = true;
	bool below = true;
	for (int i = 0; i < n_P; i++)
	{
		if (X[2 * i + 1] >= a_u * X[2 * i] + b_u)
		{
			above = false;
			break;
		}
		if (X[2 * i + 1] < a_l * X[2 * i] + b_l)
		{
			below = false;
			break;
		}
	}
	if (above)
	{
		test(a_u, b_u);
	}
	if (below)
	{
		test(a_l, b_l);
	}
	bool test_borders = false;
	for (int i = 1; i < n_P; i++)
	{
		p2_x = X[2 * i];
		p2_y = X[2 * i + 1];
		if (p2_x != p1_x || p2_y != p1_y)
		{
			m_x = (p1_x + p2_x) / 2.0;
			m_y = (p1_y + p2_y) / 2.0;
			bool inside = false;
			if (m_x >= 0)
			{
				if (m_y < amax * m_x + bmax && m_y >= amin * m_x + bmin)
				{
					inside = true;
				}
			}
			else
			{
				if (m_y < amin * m_x + bmax && m_y >= amax * m_x + bmin)
				{
					inside = true;
				}
			}
			if (inside)
			{
				t_a = amin;
				t_b = m_y - t_a * m_x;
				if (t_b <= bmax && t_b >= bmin && t_a <= amax && t_a >= amin)
				{
					test(t_a, t_b);
				}
				t_a = amax;
				t_b = m_y - t_a * m_x;
				if (t_b <= bmax && t_b >= bmin && t_a <= amax && t_a >= amin)
				{
					test(t_a, t_b);
				}
				t_b = bmin;
				t_a = (m_y - t_b) / m_x;
				if (t_b <= bmax && t_b >= bmin && t_a <= amax && t_a >= amin)
				{
					test(t_a, t_b);
				}
				t_b = bmax;
				t_a = (m_y - t_b) / m_x;
				if (t_b <= bmax && t_b >= bmin && t_a <= amax && t_a >= amin)
				{
					test(t_a, t_b);
				}
			}
			else
			{
				test_borders = true;
			}
			p1_x = p2_x;
			p1_y = p2_y;
		}
	}
	if (test_borders)
	{
		t_a = amin;
		t_b = bmin;
		test(t_a, t_b);
		t_a = amin;
		t_b = bmax;
		test(t_a, t_b);
		t_a = amax;
		t_b = bmin;
		test(t_a, t_b);
		t_a = amax;
		t_b = bmax;
		test(t_a, t_b);
	}
}
void BBNode::test_4()
{
	double t_a = amin;
	double t_b = bmin;
	test(t_a, t_b);
	t_a = amax;
	t_b = bmax;
	test(t_a, t_b);
	t_a = amin;
	t_b = bmax;
	test(t_a, t_b);
	t_a = amax;
	t_b = bmin;
	test(t_a, t_b);
}
