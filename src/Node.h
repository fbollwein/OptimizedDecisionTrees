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

#include <vector>
#include "Condition.h"
#include "Rule.h"
#include "ctpl_stl.h"
#include "Environment.h"
using namespace std;
#ifndef NODE
#define NODE
class Node
{
public:
	int *Y;
	double *Y_reg;
	double **X;
	int ***dataY;
	double ***dataY_reg;
	double *Ni;
	double opt_gap = 0;
	double time = 0;
	bool regression = false;
	double **Xnominal;
	vector<int> nominal_domain;
	vector<vector<bool>> is_nominal_val_in_dataset;
	int total_numeric_features = 0;
	int total_nominal_features = 0;
	Environment *env;
	int total_features;
	int total_points;
	int n_classes;
	bool isLeaf;
	bool hasData;
	int id;
	string str_id = "";
	int maxVote;
	string maxVote_str = "";
	double impurity;
	bool shouldSplit;
	bool shouldSplit_numeric;
	bool shouldSplit_nominal;
	int random_state;
	vector<double> a;
	double mean;
	double mse;
	double median;
	double mae;
	double LB;
	int bestInd1;
	int bestInd2;
	double bestA;
	double bestB;
	double bestX1;
	double bestX2;
	int type;
	int split;
	double errorC;
	double classification_error;
	double sq_error;
	double a_error;
	int depth;
	Condition cond;
	double impS1;
	double impS2;
	double impS;
	ctpl::thread_pool *p;
	bool isInTree;
	bool possesses_data = true;
	Node();
	Node(double **X, int *Y, double *Y_reg, int total_points, int total_numeric_features, int n_classes, ctpl::thread_pool *p, Environment *env);
	Node(double **X, double **Xnominal, int *Y, double *Y_reg, int total_points, int total_numeric_features, int total_nominal_features, vector<int> nominal_domain, vector<vector<bool>> is_nominal_val_in_dataset, int n_classes, ctpl::thread_pool *p, Environment *env);
	void initialize(double **X, double **Xnominal, int *Y, double *Y_reg, int total_points, int total_numeric_features, int total_nominal_features, vector<int> nominal_domain, vector<vector<bool>> is_nominal_val_in_dataset, int n_classes, ctpl::thread_pool *p, Environment *env);
	Node(double **X, int *Y, double *Y_reg, int ***dataY, double ***dataY_reg, int total_points, int total_numeric_features, int n_classes, ctpl::thread_pool *p, Environment *env);
	Node(double **X, double **Xnominal, int *Y, double *Y_reg, int ***dataY, double ***dataY_reg, int total_points, int total_numeric_features, int total_nominal_features, vector<int> nominal_domain, vector<vector<bool>> is_nominal_val_in_dataset, int n_classes, ctpl::thread_pool *p, Environment *env);
	void initialize(double **X, double **Xnominal, int *Y, double *Y_reg, int ***dataY, double ***dataY_reg, int total_points, int total_numeric_features, int total_nominal_features, vector<int> nominal_domain, vector<vector<bool>> is_nominal_val_in_dataset, int n_classes, ctpl::thread_pool *p, Environment *env);
	~Node();
	void getData();
	void getDataOC1();
	vector<Node *> branch();
	vector<Node *> branch_randomly();
	vector<Node *> divide(Rule rule);
	vector<Node *> distribute(vector<Node *> nodes);
	Rule branchNominal(vector<int> features, int max_features);
	Rule branchBivariate(Rule rule, int max_features);
	Rule branchBivariate(vector<int> features, double bound);
	Rule branchMultivariate(Rule rule, int max_features);
	Rule branchMultivariate(vector<int> features);
	Rule branchUnivariate(vector<int> features, int max_features);
	Rule branchCross(Rule rule, int max_features);
	Rule beautify(Rule &rule, bool max_margin);
	Rule chooseBetter(Rule &rule);
	bool evaluate(vector<double> &x, vector<double> &x_nominal);
	int classify();
	void setCondition(Condition cond);
	double intrinsicInformation(Rule &rule);
	vector<int> getLinearIndependentFeatures(Rule &rule, double precision, vector<int> features, bool use_order);
	Node *copy(int id, double **X, double **Xnominal, int *Y, double *Y_reg, Environment *env, bool take_foreign, bool data);
	void deleteData();
	void deleteDataOnly();
};
#endif
