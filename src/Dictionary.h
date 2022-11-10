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

using namespace std;
#include "Eta.h"
#include "Permutation.h"
#include "Environment.h"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Vector_reg.h"
class Dictionary
{
public:
    Environment *env;
    Eigen::VectorXd orientation;
    vector<bool> is_duplicate;
    Eigen::VectorXd b;
    Eigen::VectorXd c;
    Eigen::VectorXd cB;
    Eigen::VectorXd cN;
    Eigen::VectorXd red_c;
    Eigen::VectorXd y;
    Eigen::VectorXd w;
    Eigen::VectorXd xB;
    Eigen::VectorXd xN;
    Eigen::VectorXd x;
    Eigen::VectorXd a;
    Eigen::MatrixXd A;
    Eigen::MatrixXd W;
    Eigen::VectorXd lb_fixed;
    Eigen::VectorXd ub_fixed;
    Eigen::VectorXd lb;
    Eigen::VectorXd ub;
    vector<Eigen::Triplet<double>> coefficients;
    Eigen::SparseMatrix<double> sparse_AB;
    Eigen::SparseMatrix<double> sparse_AB_T;
    Eigen::MatrixXd AB;
    Eigen::MatrixXd AN_T;
    Eigen::SparseMatrix<double> sparse_AN_T;
    int d;
    int n;
    vector<int> B;
    vector<int> N;
    Permutation lex_order;
    Permutation lex_order_non_slack;
    vector<int> lex_order_nb;
    vector<int> lex_order_non_slack_nb;
    vector<int> orientation_non_slack;
    bool use_lex = true;
    vector<bool> is_basic;
    vector<bool> is_redundant;
    vector<bool> is_redundant_save;
    vector<bool> hasW;
    vector<bool> hasW_save;
    bool is_degenerate = false;
    bool is_shiftable = false;
    bool is_basis_feasible = false;
    bool bland = true;
    bool has_zero_basic_slack = false;
    vector<int> B_tmp;
    vector<int> N_tmp;
    int e_tmp;
    int l_tmp;
    Eigen::VectorXd orientation_tmp;
    std::map<int, vector<int>> duplicates;
    vector<int> position_B;
    vector<int> position_N;
    Eigen::VectorXd x_shifted;
    bool has_shifted = false;
    vector<Eta> eta;
    int max_eta = 100;
    double deg_tol = 1e-6;
    double feas_tol = 1e-6;
    double dual_feas_tol = 1e-6;
    double opt_requirement = 1e-9;
    double zero_precision = 3e-10;
    double comp_precision = 1e-12;
    double comp_zero_precision = 1e-9;
    double res_tol = 5e-14;
    double piv_tol = 1e-6;
    bool sparse = true;
    bool has_sparse = false;
    bool sparseN = true;
    bool pricing = false;
    bool has_LU = false;
    Eigen::SparseLU<Eigen::SparseMatrix<double>> dec_sp;
    Eigen::SparseLU<Eigen::SparseMatrix<double>> dec_sp_T;
    Eigen::PartialPivLU<Eigen::MatrixXd> dec;
    int n_classes;
    int *Y;
    double *Y_reg;
    vector<double> distribution;
    vector<Vector_reg> reg;
    bool scale = true;
    int scale_type = 0;
    vector<double> row_scale;
    vector<double> col_scale;
    Dictionary(Eigen::MatrixXd &A, Eigen::VectorXd &b, Eigen::VectorXd &c, Eigen::VectorXd &orientation, vector<int> &B, Environment *env);
    Dictionary(Eigen::MatrixXd &A, Eigen::VectorXd &b, Eigen::VectorXd &c, Eigen::VectorXd &orientation, vector<int> &B, int *Y, double *Y_reg, int n_classes, std::map<int, vector<int>> &duplicates, Environment *env);
    Dictionary(Dictionary &dict, bool computeBasis);
    void init(Eigen::MatrixXd &A, Eigen::VectorXd &b, Eigen::VectorXd &c, Eigen::VectorXd &orientation, vector<int> &B, Environment *env);
    void init_Bounds();
    bool enablePricing();
    bool disablePricing();
    int recompute_decomposition(bool update);
    Eigen::VectorXd reduced_costs();
    int pivot(int e, int type, int sign, bool use_order);
    int pivot_harris(int e, int type, int sign);
    vector<pair<int, pair<double, int>>> find_pivot(int e_pos, int type, int sign);
    int get_zero_slack(int method);
    int get_cobasic_slack(int method);
    int get_cobasic_var_rand(int method, bool exclude_last);
    int get_cobasic_var_price(int method, bool exclude_last);
    vector<int> get_Duplicates(int i);
    Eigen::VectorXd solve(Eigen::VectorXd rhs);
    Eigen::VectorXd solve_T(Eigen::VectorXd lhs);
    void setBasis(vector<int> &B2, Eigen::VectorXd &x2, const Permutation &order);
    void preparePhaseI(bool feasibility_check);
    void prepareFeasibilityCheck();
    void setOrientation(Eigen::VectorXd &orientation_tmp);
    void updateStates();
    void updateSol();
    bool is_inner(Eigen::VectorXd x_tmp);
    void reset_Redundancy();
    int compare(double a, double b, double prec);
    int lex_compare(const pair<int, pair<double, int>> &cand1, const pair<int, pair<double, int>> &cand2);
    int lex_compare_push(const pair<int, pair<double, int>> &cand1, const pair<int, pair<double, int>> &cand2, Eigen::VectorXd &w2);
    int lex_compare_push_co(int i_push, const pair<int, pair<double, int>> &cand, Eigen::VectorXd &w2);
    void change_orientation(int i, bool revert);
    int update_orientation(int i_push, int i_change);
    int push_slack(int i);
    void change_side(int i);
    bool is_lex_pos();
    Eigen::VectorXd getUnscaled();
    void shift();
    double getResidual();
};