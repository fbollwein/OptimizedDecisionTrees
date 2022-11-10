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

#include <string>
#include <random>
#include "FunctionParser.h"
#ifdef use_gurobi
#include "gurobi_c++.h"
extern GRBEnv grb_env;
#endif
using namespace std;
#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H
class Environment
{
public:
	static const int gini_imp = 0;
	static const int entropy_imp = 1;
	static const int twoing = -1;
	static const int mse = 2;
	static const int mae = 3;
	static const int geom = 0;
	static const int equi = 1;
	static const int no_scale = -1;
	int criterion;
	int criterion_reg;
	int max_depth;
	double min_samples_split;
	double min_cases;
	double min_weight_fraction_leaf;
	int max_features;
	int max_leaf_nodes;
	double min_impurity_decrease;
	double min_impurity_split;
	bool univariate;
	bool mult;
	bool cross;
	bool bi;
	bool lookahead;
	int random_state;
	double cf;
	int n_threads;
	int oblique_split_threshold;
	double nom_timelimit;
	double ga_timelimit;
	bool cross_entropy;
	bool simplex;
	bool branch_randomly;
	int nominal_partitions;
	bool regression;
	bool normalize;
	bool max_margin;
	int margin_norm;
	int cross_entropy_samples;
	string cross_entropy_samples_form;
	FunctionParser cross_entropy_samples_expression;
	double cross_entropy_alpha;
	double cross_entropy_rho;
	int cross_entropy_no_improvement;
	string cross_entropy_no_improvement_form;
	FunctionParser cross_entropy_no_improvement_expression;
	int spx_start_iterations;
	string spx_start_iterations_form;
	FunctionParser spx_start_iterations_expression;
	double spx_beta_T;
	string spx_no_improvement_form;
	int spx_no_improvement;
	FunctionParser spx_no_improvement_expression;
	double spx_deg_tol;
	double spx_feas_tol;
	double spx_dual_feas_tol;
	double spx_opt_requirement;
	double spx_zero_precision;
	double spx_comp_precision;
	double spx_comp_zero_precision;
	double spx_res_tol;
	double spx_piv_tol;
	int spx_scale;
	int spx_max_eta;
	bool deleteData = false;
	std::default_random_engine generator;
	std::uniform_real_distribution<double> uniform;
	std::uniform_int_distribution<> uniform_int;
	bool check_residual;
	Environment();
	void copy(Environment &env);
	void seed(int random_state);
	double Rand();
	double Rand(double lb, double ub);
	int rand();
};
#endif
