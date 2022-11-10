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
#include "Result.h"
#include "General.h"
#include "Vector_reg.h"
#ifndef BBNODE
#define BBNODE
class BBNode
{
public:
	int n_classes;
	int total_points;
	double LB;
	double UB;
	double *O;
	int n_P;
	bool defined;
	double *o;
	doubleIntPair *dists;
	double amin;
	double amax;
	double bmin;
	double bmax;
	double *X;
	int *Y;
	bool regression;
	double *Y_reg;
	Vector_reg *reg1;
	Vector_reg *reg2;
	Vector_reg *reg3;
	Vector_reg *reg1_tmp;
	Vector_reg *reg2_tmp;
	bool same_as_parent;
	int status;
	double eps;
	double a;
	double b;
	double branch_a;
	double branch_b;
	int type;
	int criterion;
	BBNode();
	BBNode(double *X, int *Y, double *Y_reg, double amin, double amax, double bmin, double bmax,
		   double *O, double *o, Vector_reg *reg1, Vector_reg *reg2, Vector_reg *reg3, int n_P, bool same_as_parent, double parentLB, int n_classes, int total_points, doubleIntPair *dists, int type, bool regression, int criterion);
	bool singleClassTest();
	void upperBound();
	void upperBound_mid();
	void upperBound(double a, double b);
	void find_Bound(double *o1);
	void simpleBound();
	void lowerBound();
	void lowerBound_reg();
	void lowerBound2();
	void lowerBound3();
	void reduce();
	void process();
	void test_4();
	void test_all();
	void test(double a, double b);
	bool computeBranchingRule();
	bool computeBranchingRule_B(double *o1, double *o2);
	bool computeBranchingRule_A(double *o1, double *o2);
	void computeBranchingRule_B_Bound(bool switched);
	void computeBranchingRule_A_Bound(bool switched);
};
#endif
