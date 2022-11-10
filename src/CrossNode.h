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
#include "Vector_reg.h"
#ifndef CROSSNODE
#define CROSSNODE
class CrossNode
{
public:
	int n_classes;
	int total_points;
	double LB;
	double UB;
	double *N;
	Vector_reg **N_reg;
	int *Q;
	Vector_reg **Q_reg;
	int P[4];
	int *cX1;
	int *cX2;
	double x;
	double y;
	int **dataX1;
	int **dataX2;
	double **dataX1_reg;
	double **dataX2_reg;
	int pointer[4];
	double xmin;
	double xmax;
	double ymin;
	double ymax;
	double *X1;
	double *X2;
	int criterion;
	bool regression;
	CrossNode();
	CrossNode(double X1[], double X2[], double *N, Vector_reg **N_reg, int *P,
			  int n_classes, int total_points, int **dataX1, int **dataX2, double **dataX1_reg, double **dataX2_reg,
			  int *cX1, int *cX2, bool regression, int criterion);
	void getSectors_opt(int *Q, Vector_reg *Q_reg[12]);
	void upperBound();
	void simpleBound();
	void getBestXonLine(int countY);
	void getBestYonLine(int countX);
	void getPoint(int countX1, int countX2);
	double getImpN();
	bool checkN();
};
#endif
