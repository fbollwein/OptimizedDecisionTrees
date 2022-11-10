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
#include <mutex>
#include "CrossNode.h"
#include "PointValuePair.h"
#include "Result.h"
#include "CrossBranchAndBound.h"
using namespace std;
std::mutex mtx4;
class CompareNodes
{
public:
	bool operator()(CrossNode *a, CrossNode *b)
	{
		return a->LB > b->LB;
	}
};
void CrossBranchAndBound(int id, double X1[], double X2[], int ind1, int ind2, int n_classes, int total_points,
						 int **dataX1, int **dataX2, double **dataX1_reg, double **dataX2_reg, int cX1[], int cX2[], int countX, int countY, Result *res, Environment *env)
{
	stack<CrossNode *> que;
	int *Q = new int[12 * n_classes];
	Vector_reg **Q_reg = new Vector_reg *[12];
	if (env->regression)
	{
		for (int i = 0; i < 12; i++)
		{
			Q_reg[i] = new Vector_reg();
			if (env->criterion == Environment::mae)
			{
				Q_reg[i]->use_mae();
			}
		}
	}
	int *newP = new int[4];
	double *newN2;
	Vector_reg **newN_reg2;
	if (!env->regression)
	{
		newN2 = new double[4 * n_classes];
		for (int k = 0; k < 4; k++)
		{
			for (int m = 0; m < n_classes; m++)
			{
				newN2[k * n_classes + m] = 0;
			}
		}
	}
	else
	{
		newN_reg2 = new Vector_reg *[4];
		for (int i = 0; i < 4; i++)
		{
			newN_reg2[i] = new Vector_reg();
			if (env->criterion == Environment::mae)
			{
				newN_reg2[i]->use_mae();
			}
		}
	}
	newP[0] = 0;
	newP[1] = countX - 1;
	newP[2] = 0;
	newP[3] = countY - 1;
	CrossNode *root = new CrossNode(X1, X2, newN2, newN_reg2, newP, n_classes, total_points, dataX1, dataX2, dataX1_reg, dataX2_reg, cX1, cX2, env->regression, env->criterion);
	que.push(root);
	root->simpleBound();
	int maxNodes = 0;
	CrossNode v;
	bool flag;
	int f = 0;
	while (que.size() > 0)
	{
		CrossNode *vptr = que.top();
		v = *vptr;
		que.pop();
		v.simpleBound();
		maxNodes++;
		if (v.LB < res->bestImp)
		{
			if (v.P[1] - v.P[0] <= 1)
			{
				if (v.P[3] - v.P[2] <= 1)
				{
					if (v.P[1] != countX - 1 || v.P[3] != countY - 1)
					{
						v.getPoint(countX, countY);
						mtx4.lock();
						if (v.UB < res->bestImp)
						{
							res->split_found = true;
							res->bestImp = v.UB;
							res->bestX1 = v.x;
							res->bestX2 = v.y;
							res->bestInd1 = ind1;
							res->bestInd2 = ind2;
							res->isCross = true;
						}
						mtx4.unlock();
					}
				}
				else
				{
					v.getBestYonLine(countX);
					mtx4.lock();
					if (v.UB < res->bestImp)
					{
						res->split_found = true;
						res->bestImp = v.UB;
						res->bestX1 = v.x;
						res->bestX2 = v.y;
						res->bestInd1 = ind1;
						res->bestInd2 = ind2;
						res->isCross = true;
					}
					mtx4.unlock();
				}
			}
			else
			{
				if (v.P[3] - v.P[2] <= 1)
				{
					v.getBestXonLine(countY);
					mtx4.lock();
					if (v.UB < res->bestImp)
					{
						res->split_found = true;
						res->bestImp = v.UB;
						res->bestX1 = v.x;
						res->bestX2 = v.y;
						res->bestInd1 = ind1;
						res->bestInd2 = ind2;
						res->isCross = true;
					}
					mtx4.unlock();
				}
				else
				{
					v.getSectors_opt(Q, Q_reg);
					v.upperBound();
					mtx4.lock();
					if (v.UB < res->bestImp)
					{
						res->split_found = true;
						res->bestImp = v.UB;
						res->bestX1 = v.x;
						res->bestX2 = v.y;
						res->bestInd1 = ind1;
						res->bestInd2 = ind2;
						res->isCross = true;
					}
					mtx4.unlock();
					if (v.LB <= v.UB)
					{
						newP[0] = v.P[0];
						newP[1] = (v.P[0] + v.P[1]) / 2;
						newP[2] = v.P[2];
						newP[3] = (v.P[2] + v.P[3]) / 2;
						double *newN;
						Vector_reg **newN_reg;
						if (env->regression)
						{
							newN_reg = new Vector_reg *[4];
							for (int i = 0; i < 4; i++)
							{
								newN_reg[i] = new Vector_reg();
								if (env->criterion == Environment::mae)
								{
									newN_reg[i]->use_mae();
								}
							}
						}
						else
						{
							newN = new double[4 * n_classes];
						}
						for (int k = 0; k < 4; k++)
						{
							if (k == 0)
							{
								if (!env->regression)
								{
									for (int m = 0; m < n_classes; m++)
									{
										newN[k * n_classes + m] = v.N[k * n_classes + m];
									}
								}
								else
								{
									newN_reg[k]->add(v.N_reg[k]);
								}
							}
							else if (k == 1)
							{
								if (!env->regression)
								{
									for (int m = 0; m < n_classes; m++)
									{
										newN[k * n_classes + m] = (v.N[k * n_classes + m] + v.Q[1 * n_classes + m]);
									}
								}
								else
								{
									newN_reg[k]->add(v.N_reg[k]);
									newN_reg[k]->add(v.Q_reg[1]);
								}
							}
							else if (k == 2)
							{
								if (!env->regression)
								{
									for (int m = 0; m < n_classes; m++)
									{
										newN[k * n_classes + m] = (v.N[k * n_classes + m] + v.Q[6 * n_classes + m]);
									}
								}
								else
								{
									newN_reg[k]->add(v.N_reg[k]);
									newN_reg[k]->add(v.Q_reg[6]);
								}
							}
							else if (k == 3)
							{
								if (!env->regression)
								{
									for (int m = 0; m < n_classes; m++)
									{
										newN[k * n_classes + m] = (v.N[k * n_classes + m] + v.Q[8 * n_classes + m] + v.Q[9 * n_classes + m] + v.Q[11 * n_classes + m]);
									}
								}
								else
								{
									newN_reg[k]->add(v.N_reg[k]);
									newN_reg[k]->add(v.Q_reg[8]);
									newN_reg[k]->add(v.Q_reg[9]);
									newN_reg[k]->add(v.Q_reg[11]);
								}
							}
						}
						CrossNode *u1 = new CrossNode(X1, X2, newN, newN_reg, newP, n_classes, total_points, dataX1, dataX2, dataX1_reg, dataX2_reg, cX1, cX2, env->regression, env->criterion);
						que.push(u1);
						newP[0] = (v.P[0] + v.P[1]) / 2;
						newP[1] = v.P[1];
						newP[2] = v.P[2];
						newP[3] = (v.P[2] + v.P[3]) / 2;
						if (env->regression)
						{
							newN_reg = new Vector_reg *[4];
							for (int i = 0; i < 4; i++)
							{
								newN_reg[i] = new Vector_reg();
								if (env->criterion == Environment::mae)
								{
									newN_reg[i]->use_mae();
								}
							}
						}
						else
						{
							newN = new double[4 * n_classes];
						}
						for (int k = 0; k < 4; k++)
						{
							if (k == 1)
							{
								if (!env->regression)
								{
									for (int m = 0; m < n_classes; m++)
									{
										newN[k * n_classes + m] = (v.N[k * n_classes + m]);
									}
								}
								else
								{
									newN_reg[k]->add(v.N_reg[k]);
								}
							}
							else if (k == 0)
							{
								if (!env->regression)
								{
									for (int m = 0; m < n_classes; m++)
									{
										newN[k * n_classes + m] = (v.N[k * n_classes + m] + v.Q[0 * n_classes + m]);
									}
								}
								else
								{
									newN_reg[k]->add(v.N_reg[k]);
									newN_reg[k]->add(v.Q_reg[0]);
								}
							}
							else if (k == 3)
							{
								if (!env->regression)
								{
									for (int m = 0; m < n_classes; m++)
									{
										newN[k * n_classes + m] = (v.N[k * n_classes + m] + v.Q[9 * n_classes + m]);
									}
								}
								else
								{
									newN_reg[k]->add(v.N_reg[k]);
									newN_reg[k]->add(v.Q_reg[9]);
								}
							}
							else if (k == 2)
							{
								if (!env->regression)
								{
									for (int m = 0; m < n_classes; m++)
									{
										newN[k * n_classes + m] = (v.N[k * n_classes + m] + v.Q[6 * n_classes + m] + v.Q[7 * n_classes + m] + v.Q[10 * n_classes + m]);
									}
								}
								else
								{
									newN_reg[k]->add(v.N_reg[k]);
									newN_reg[k]->add(v.Q_reg[6]);
									newN_reg[k]->add(v.Q_reg[7]);
									newN_reg[k]->add(v.Q_reg[10]);
								}
							}
						}
						CrossNode *u2 = new CrossNode(X1, X2, newN, newN_reg, newP, n_classes, total_points, dataX1, dataX2, dataX1_reg, dataX2_reg, cX1, cX2, env->regression, env->criterion);
						que.push(u2);
						newP[0] = v.P[0];
						newP[1] = (v.P[0] + v.P[1]) / 2;
						newP[2] = (v.P[2] + v.P[3]) / 2;
						newP[3] = v.P[3];
						if (env->regression)
						{
							newN_reg = new Vector_reg *[4];
							for (int i = 0; i < 4; i++)
							{
								newN_reg[i] = new Vector_reg();
								if (env->criterion == Environment::mae)
								{
									newN_reg[i]->use_mae();
								}
							}
						}
						else
						{
							newN = new double[4 * n_classes];
						}
						for (int k = 0; k < 4; k++)
						{
							if (k == 2)
							{
								if (!env->regression)
								{
									for (int m = 0; m < n_classes; m++)
									{
										newN[k * n_classes + m] = (v.N[k * n_classes + m]);
									}
								}
								else
								{
									newN_reg[k]->add(v.N_reg[k]);
								}
							}
							else if (k == 0)
							{
								if (!env->regression)
								{
									for (int m = 0; m < n_classes; m++)
									{
										newN[k * n_classes + m] = (v.N[k * n_classes + m] + v.Q[2 * n_classes + m]);
									}
								}
								else
								{
									newN_reg[k]->add(v.N_reg[k]);
									newN_reg[k]->add(v.Q_reg[2]);
								}
							}
							else if (k == 3)
							{
								if (!env->regression)
								{
									for (int m = 0; m < n_classes; m++)
									{
										newN[k * n_classes + m] = (v.N[k * n_classes + m] + v.Q[11 * n_classes + m]);
									}
								}
								else
								{
									newN_reg[k]->add(v.N_reg[k]);
									newN_reg[k]->add(v.Q_reg[11]);
								}
							}
							else if (k == 1)
							{
								if (!env->regression)
								{
									for (int m = 0; m < n_classes; m++)
									{
										newN[k * n_classes + m] = (v.N[k * n_classes + m] + v.Q[1 * n_classes + m] + v.Q[4 * n_classes + m] + v.Q[5 * n_classes + m]);
									}
								}
								else
								{
									newN_reg[k]->add(v.N_reg[k]);
									newN_reg[k]->add(v.Q_reg[1]);
									newN_reg[k]->add(v.Q_reg[4]);
									newN_reg[k]->add(v.Q_reg[5]);
								}
							}
						}
						CrossNode *u3 = new CrossNode(X1, X2, newN, newN_reg, newP, n_classes, total_points, dataX1, dataX2, dataX1_reg, dataX2_reg, cX1, cX2, env->regression, env->criterion);
						que.push(u3);
						newP[0] = (v.P[0] + v.P[1]) / 2;
						newP[1] = v.P[1];
						newP[2] = (v.P[2] + v.P[3]) / 2;
						newP[3] = v.P[3];
						if (env->regression)
						{
							newN_reg = new Vector_reg *[4];
							for (int i = 0; i < 4; i++)
							{
								newN_reg[i] = new Vector_reg();
								if (env->criterion == Environment::mae)
								{
									newN_reg[i]->use_mae();
								}
							}
						}
						else
						{
							newN = new double[4 * n_classes];
						}
						for (int k = 0; k < 4; k++)
						{
							if (k == 3)
							{
								if (!env->regression)
								{
									for (int m = 0; m < n_classes; m++)
									{
										newN[k * n_classes + m] = (v.N[k * n_classes + m]);
									}
								}
								else
								{
									newN_reg[k]->add(v.N_reg[k]);
								}
							}
							else if (k == 1)
							{
								if (!env->regression)
								{
									for (int m = 0; m < n_classes; m++)
									{
										newN[k * n_classes + m] = (v.N[k * n_classes + m] + v.Q[5 * n_classes + m]);
									}
								}
								else
								{
									newN_reg[k]->add(v.N_reg[k]);
									newN_reg[k]->add(v.Q_reg[5]);
								}
							}
							else if (k == 2)
							{
								if (!env->regression)
								{
									for (int m = 0; m < n_classes; m++)
									{
										newN[k * n_classes + m] = (v.N[k * n_classes + m] + v.Q[10 * n_classes + m]);
									}
								}
								else
								{
									newN_reg[k]->add(v.N_reg[k]);
									newN_reg[k]->add(v.Q_reg[10]);
								}
							}
							else if (k == 0)
							{
								if (!env->regression)
								{
									for (int m = 0; m < n_classes; m++)
									{
										newN[k * n_classes + m] = (v.N[k * n_classes + m] + v.Q[0 * n_classes + m] + v.Q[2 * n_classes + m] + v.Q[3 * n_classes + m]);
									}
								}
								else
								{
									newN_reg[k]->add(v.N_reg[k]);
									newN_reg[k]->add(v.Q_reg[0]);
									newN_reg[k]->add(v.Q_reg[2]);
									newN_reg[k]->add(v.Q_reg[3]);
								}
							}
						}
						CrossNode *u4 = new CrossNode(X1, X2, newN, newN_reg, newP, n_classes, total_points, dataX1, dataX2, dataX1_reg, dataX2_reg, cX1, cX2, env->regression, env->criterion);
						que.push(u4);
					}
				}
			}
		}
		if (env->regression)
		{
			for (int l = 0; l < 4; l++)
			{
				delete v.N_reg[l];
			}
			delete[] v.N_reg;
		}
		else
		{
			delete[] v.N;
		}
		delete vptr;
	}
	delete[] newP;
	delete[] Q;
	if (env->regression)
	{
		for (int i = 0; i < 12; i++)
		{
			delete Q_reg[i];
		}
	}
	delete[] Q_reg;
}
