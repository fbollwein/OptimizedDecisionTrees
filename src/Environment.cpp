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

#include "Environment.h"
#include <limits>
#include <chrono>
#include <iostream>
#ifdef use_gurobi
GRBEnv grb_env=GRBEnv();
#endif
Environment::Environment()
{
  max_depth = std::numeric_limits<int>::max();
  criterion = gini_imp;
  criterion_reg = mse;
  min_samples_split = 2;
  min_cases = 1;
  min_weight_fraction_leaf = 0.0;
  max_features = std::numeric_limits<int>::max();
  max_leaf_nodes = std::numeric_limits<int>::max();
  min_impurity_decrease = -1;
  min_impurity_split = -1;
  random_state = std::chrono::system_clock::now().time_since_epoch().count();
  univariate = true;
  mult = false;
  bi = false;
  cross = false;
  lookahead = true;
  n_threads = 0;
  oblique_split_threshold = 0;
  nom_timelimit = std::numeric_limits<int>::max();
  ga_timelimit = std::numeric_limits<int>::max();
  cross_entropy = false;
  simplex = false;
  nominal_partitions = 2;
  regression = false;
  normalize = true;
  max_margin = false;
  margin_norm=1;
  cross_entropy_samples = -1;
  cross_entropy_samples_form = "";
  cross_entropy_alpha = 0.2;
  cross_entropy_rho = 0.1;
  cross_entropy_no_improvement = 3;
  cross_entropy_no_improvement_form="";
	cross_entropy_no_improvement_expression;
	spx_start_iterations=100;
	spx_start_iterations_form="";
	spx_start_iterations_expression;
	spx_beta_T=0.85;
  spx_no_improvement=-1;
	spx_no_improvement_form="";
	spx_no_improvement_expression;
  spx_deg_tol = 1e-6;
  spx_feas_tol = 1e-6;
  spx_dual_feas_tol = 1e-6;
  spx_opt_requirement = 1e-9;
  spx_zero_precision = 3e-10;
  spx_comp_precision = 1e-12;
  spx_comp_zero_precision = 1e-9;
  spx_res_tol = 5e-14;
  spx_piv_tol = 1e-6;
  spx_scale = 0;
  spx_max_eta = 100;
  generator = std::default_random_engine(random_state);
  uniform = std::uniform_real_distribution<double>(0.0, 1.0);
  uniform_int = std::uniform_int_distribution<int>(0, std::numeric_limits<int>::max());
  branch_randomly = false;
  check_residual = false;
}
void Environment::seed(int random_state)
{
  this->random_state = random_state;
  generator.seed(random_state);
}
double Environment::Rand()
{
  return uniform(generator);
}
double Environment::Rand(double lb, double ub)
{
  return lb + uniform(generator) * (ub - lb);
}
int Environment::rand()
{
  return uniform_int(generator);
}
void Environment::copy(Environment &env)
{
}
