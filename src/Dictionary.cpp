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

#include "Dictionary.h"
#include "General.h"
using namespace std;
#include <iostream>
Dictionary::Dictionary(Eigen::MatrixXd &A_tmp, Eigen::VectorXd &b, Eigen::VectorXd &c, Eigen::VectorXd &orientation, vector<int> &B, Environment *env)
{
    init(A_tmp, b, c, orientation, B, env);
}
void Dictionary::setBasis(vector<int> &B2, Eigen::VectorXd &x2, const Permutation &order)
{
    this->lex_order = Permutation(order);
    hasW = vector<bool>(d + n, false);
    this->B = vector<int>(B2);
    N.clear();
    is_basic = vector<bool>(A.cols(), false);
    position_B = vector<int>(A.cols());
    position_N = vector<int>(A.cols());
    x = Eigen::VectorXd::Zero(A.cols());
    xB = Eigen::VectorXd::Zero(B.size());
    xN = Eigen::VectorXd::Zero(A.cols() - B.size());
    for (int j_pos = 0; j_pos < B.size(); j_pos++)
    {
        is_basic[B[j_pos]] = true;
        position_B[B[j_pos]] = j_pos;
    }
    for (int k = 0; k < A.cols(); k++)
    {
        if (!is_basic[k])
        {
            N.push_back(k);
            position_N[k] = N.size() - 1;
        }
    }
    cB = Eigen::VectorXd::Zero(B.size());
    cN = Eigen::VectorXd::Zero(N.size());
    for (int j_pos = 0; j_pos < B.size(); j_pos++)
    {
        xB(j_pos) = x2(B[j_pos]);
        x(B[j_pos]) = x2(B[j_pos]);
        cB(j_pos) = c(B[j_pos]);
    }
    for (int j_pos = 0; j_pos < N.size(); j_pos++)
    {
        xN(j_pos) = x2(N[j_pos]);
        x(N[j_pos]) = x2(N[j_pos]);
        cN(j_pos) = c(N[j_pos]);
    }
    AB = Eigen::MatrixXd::Zero(n, n);
    sparse_AB = Eigen::SparseMatrix<double>(n, n);
    recompute_decomposition(true);
}
void Dictionary::init_Bounds()
{
    lb = Eigen::VectorXd::Zero(A.cols());
    ub = Eigen::VectorXd::Zero(A.cols());
    for (int j = 0; j < A.cols(); j++)
    {
        if (j < d)
        {
            lb(j) = std::numeric_limits<double>::lowest();
            ub(j) = std::numeric_limits<double>::max();
        }
        else
        {
            lb(j) = 0;
            ub(j) = std::numeric_limits<double>::max();
        }
    }
    lb_fixed = lb;
    ub_fixed = ub;
}
Dictionary::Dictionary(Dictionary &dict, bool computeBasis)
{
    this->n_classes = dict.n_classes;
    this->distribution = dict.distribution;
    this->reg = dict.reg;
    this->duplicates = dict.duplicates;
    this->env = dict.env;
    deg_tol = env->spx_deg_tol;
    feas_tol = env->spx_feas_tol;
    dual_feas_tol = env->spx_dual_feas_tol;
    opt_requirement = env->spx_opt_requirement;
    zero_precision = env->spx_zero_precision;
    comp_precision = env->spx_comp_precision;
    comp_zero_precision = env->spx_comp_zero_precision;
    res_tol = env->spx_res_tol;
    piv_tol = env->spx_piv_tol;
    scale_type = env->spx_scale;
    max_eta = env->spx_max_eta;
    d = dict.d;
    n = dict.n;
    Y = dict.Y;
    Y_reg = dict.Y_reg;
    A = dict.A;
    W = dict.W;
    this->w = dict.w;
    this->hasW = dict.hasW;
    this->lex_order = Permutation(dict.lex_order);
    this->lex_order_non_slack = Permutation(dict.lex_order_non_slack);
    this->b = dict.b;
    this->c = dict.c;
    this->cB = dict.cB;
    this->cN = dict.cN;
    this->orientation = dict.orientation;
    x = dict.x;
    lb = dict.lb;
    ub = dict.ub;
    if (computeBasis)
    {
        setBasis(dict.B, dict.x, lex_order);
    }
    is_redundant = dict.is_redundant;
}
pair<double, double> range(Eigen::MatrixXd &A_scale)
{
    double range_max = std::numeric_limits<double>::lowest();
    double range_min = std::numeric_limits<double>::max();
    for (int i = 0; i < A_scale.rows(); i++)
    {
        for (int j = 0; j < A_scale.cols(); j++)
        {
            if (A_scale(i, j) != 0)
            {
                if (abs(A_scale(i, j)) > range_max)
                {
                    range_max = abs(A_scale(i, j));
                }
                if (abs(A_scale(i, j)) < range_min)
                {
                    range_min = abs(A_scale(i, j));
                }
            }
        }
    }
    return pair<double, double>(range_min, range_max);
}
void eq_scale(Eigen::MatrixXd &A_scale, Eigen::VectorXd &b, vector<double> &row_scale, vector<double> &col_scale)
{
    int n = A_scale.rows();
    int d = A_scale.cols();
    for (int i = 0; i < n; i++)
    {
        Eigen::MatrixXd A_save = Eigen::MatrixXd(A_scale);
        double row_max = std::numeric_limits<double>::lowest();
        for (int j = 0; j < d; j++)
        {
            if (abs(A_scale(i, j)) != 0 && abs(A_scale(i, j)) > row_max)
            {
                row_max = abs(A_scale(i, j));
            }
        }
        if (row_max == 0)
        {
            row_max = 1;
        }
        double f = row_max;
        row_scale[i] *= f;
        A_scale.row(i) = A_scale.row(i) / f;
        b(i) = b(i) / f;
    }
    for (int j = 0; j < d; j++)
    {
        double col_max = std::numeric_limits<double>::lowest();
        for (int i = 0; i < n; i++)
        {
            if (abs(A_scale(i, j)) != 0 && abs(A_scale(i, j)) > col_max)
            {
                col_max = abs(A_scale(i, j));
            }
        }
        if (col_max == 0)
        {
            col_max = 1;
        }
        double f = col_max;
        col_scale[j] *= f;
        A_scale.col(j) = A_scale.col(j) / f;
    }
}
void geom_scale(Eigen::MatrixXd &A_scale, Eigen::VectorXd &b, vector<double> &row_scale, vector<double> &col_scale)
{
    int n = A_scale.rows();
    int d = A_scale.cols();
    pair<double, double> m_range;
    double range_min;
    double range_max;
    for (int k = 0; k < 25; k++)
    {
        m_range = range(A_scale);
        range_min = m_range.first;
        range_max = m_range.second;
        if (range_min >= 0.1 && range_max <= 10)
        {
            break;
        }
        for (int i = 0; i < n; i++)
        {
            double row_max = std::numeric_limits<double>::lowest();
            double row_min = std::numeric_limits<double>::max();
            for (int j = 0; j < d; j++)
            {
                if (abs(A_scale(i, j)) != 0 && abs(A_scale(i, j)) > row_max)
                {
                    row_max = abs(A_scale(i, j));
                }
                if (abs(A_scale(i, j)) != 0 && abs(A_scale(i, j)) < row_min)
                {
                    row_min = abs(A_scale(i, j));
                }
            }
            double f = 0;
            if (row_min > 0 && row_max > 0)
            {
                f = sqrt(row_min * row_max);
            }
            if (f == 0)
            {
                f = 1;
            }
            row_scale[i] *= f;
            A_scale.row(i) = A_scale.row(i) / f;
            b(i) = b(i) / f;
        }
        for (int j = 0; j < d; j++)
        {
            double col_max = std::numeric_limits<double>::lowest();
            double col_min = std::numeric_limits<double>::max();
            for (int i = 0; i < n; i++)
            {
                if (abs(A_scale(i, j)) != 0 && abs(A_scale(i, j)) > col_max)
                {
                    col_max = abs(A_scale(i, j));
                }
                if (abs(A_scale(i, j)) != 0 && abs(A_scale(i, j)) < col_min)
                {
                    col_min = abs(A_scale(i, j));
                }
            }
            double f = 0;
            if (col_min > 0 && col_max > 0)
            {
                f = sqrt(col_min * col_max);
            }
            if (f == 0)
            {
                f = 1;
            }
            col_scale[j] *= f;
            A_scale.col(j) = A_scale.col(j) / f;
        }
    }
}
void pow2scale(Eigen::MatrixXd &A_scale, Eigen::VectorXd &b, vector<double> &row_scale, vector<double> &col_scale)
{
    int n = A_scale.rows();
    int d = A_scale.cols();
    double tmp;
    for (int i = 0; i < n; i++)
    {
        tmp = pow(2, floor(log2((4. / 3.) * row_scale[i])));
        if (tmp == 0)
        {
            tmp = 1;
        }
        row_scale[i] = tmp;
        A_scale.row(i) = A_scale.row(i) / tmp;
    }
    for (int j = 0; j < d; j++)
    {
        tmp = pow(2, floor(log2((4. / 3.) * col_scale[j])));
        if (tmp == 0)
        {
            tmp = 1;
        }
        col_scale[j] = tmp;
        A_scale.col(j) = A_scale.col(j) / tmp;
    }
}
void full_scale(Eigen::MatrixXd &A_scale, Eigen::VectorXd &b_scale, vector<double> &row_scale, vector<double> &col_scale, int type)
{
    if (type == Environment::no_scale)
    {
        return;
    }
    pair<double, double> m_range = range(A_scale);
    double range_min = m_range.first;
    double range_max = m_range.second;
    if (range_min >= 0.1 && range_max <= 10)
    {
        return;
    }
    Eigen::MatrixXd A_tmp = Eigen::MatrixXd(A_scale);
    Eigen::VectorXd b_tmp = Eigen::VectorXd(b_scale);
    if (type == Environment::geom)
    {
        if (range_min < 0.1 || range_max > 10)
        {
            geom_scale(A_scale, b_scale, row_scale, col_scale);
        }
    }
    m_range = range(A_scale);
    range_min = m_range.first;
    range_max = m_range.second;
    if (type == Environment::equi)
    {
        if (range_min < 0.1 || range_max > 10)
        {
            eq_scale(A_scale, b_scale, row_scale, col_scale);
        }
    }
}
void Dictionary::init(Eigen::MatrixXd &A_tmp, Eigen::VectorXd &b, Eigen::VectorXd &c, Eigen::VectorXd &orientation, vector<int> &B, Environment *env)
{
    this->env = env;
    deg_tol = env->spx_deg_tol;
    feas_tol = env->spx_feas_tol;
    dual_feas_tol = env->spx_dual_feas_tol;
    opt_requirement = env->spx_opt_requirement;
    zero_precision = env->spx_zero_precision;
    comp_precision = env->spx_comp_precision;
    comp_zero_precision = env->spx_comp_zero_precision;
    res_tol = env->spx_res_tol;
    piv_tol = env->spx_piv_tol;
    scale_type = env->spx_scale;
    if (scale_type == -1)
    {
        scale = false;
    }
    max_eta = env->spx_max_eta;
    d = A_tmp.cols();
    n = A_tmp.rows();
    row_scale = vector<double>(n, 1);
    col_scale = vector<double>(d, 1);
    Eigen::VectorXd b_scale = Eigen::VectorXd(b);
    Eigen::MatrixXd A_scale = Eigen::MatrixXd(A_tmp);
    pair<double, double> m_range = range(A_scale);
    double range_min = m_range.first;
    double range_max = m_range.second;
    full_scale(A_scale, b_scale, row_scale, col_scale, scale_type);
    m_range = range(A_scale);
    range_min = m_range.first;
    range_max = m_range.second;
    A = Eigen::MatrixXd::Zero(n, d + n + 1);
    W = Eigen::MatrixXd::Zero(n, d + n);
    hasW = vector<bool>(d + n, false);
    vector<int> order = vector<int>(A.rows());
    this->b = Eigen::VectorXd(b_scale);
    this->c = Eigen::VectorXd(c);
    for (int i = 0; i < A.rows(); i++)
    {
        for (int j = 0; j < d; j++)
        {
            A(i, j) = A_scale(i, j);
        }
        A(i, d + i) = orientation(i);
        A(i, A.cols() - 1) = orientation(i);
        order[i] = i;
    }
    shuffle(order.begin(), order.end(), env->generator);
    lex_order = Permutation(order);
    this->cB = Eigen::VectorXd::Zero(B.size());
    this->cN = Eigen::VectorXd::Zero(N.size());
    this->orientation = Eigen::VectorXd(orientation);
    x = Eigen::VectorXd::Zero(A.cols());
    init_Bounds();
    reset_Redundancy();
    setBasis(B, x, lex_order);
}
Dictionary::Dictionary(Eigen::MatrixXd &A_tmp, Eigen::VectorXd &b, Eigen::VectorXd &c, Eigen::VectorXd &orientation, vector<int> &B, int *Y, double *Y_reg, int n_classes, std::map<int, vector<int>> &duplicates, Environment *env)
{
    this->duplicates = duplicates;
    init(A_tmp, b, c, orientation, B, env);
    this->n_classes = n_classes;
    this->Y = Y;
    this->Y_reg = Y_reg;
    if (!env->regression)
    {
        distribution = vector<double>(2 * n_classes, 0);
    }
    else
    {
        reg.push_back(Vector_reg());
        reg.push_back(Vector_reg());
        if (env->criterion == Environment::mse)
        {
            reg[0].keep_vec = false;
            reg[1].keep_vec = false;
        }
        else
        {
            reg[0].use_mae();
            reg[1].use_mae();
        }
        int reservation = 0;
        for (int i = 0; i < n; i++)
        {
            reservation += duplicates[i].size();
        }
        reg[0].reserve(reservation);
        reg[1].reserve(reservation);
    }
    if (!env->regression)
    {
        for (int i = 0; i < n; i++)
        {
            if (orientation[i] > 0)
            {
                for (int k = 0; k < duplicates[i].size(); k++)
                {
                    distribution[Y[duplicates[i][k]]] += 1;
                }
            }
            else
            {
                for (int k = 0; k < duplicates[i].size(); k++)
                {
                    distribution[n_classes + Y[duplicates[i][k]]] += 1;
                }
            }
        }
    }
    else
    {
        for (int i = 0; i < n; i++)
        {
            if (orientation[i] > 0)
            {
                for (int k = 0; k < duplicates[i].size(); k++)
                {
                    reg[0].push_back(Y_reg[duplicates[i][k]]);
                }
            }
            else
            {
                for (int k = 0; k < duplicates[i].size(); k++)
                {
                    reg[1].push_back(Y_reg[duplicates[i][k]]);
                }
            }
        }
        reg[0].sort();
        reg[1].sort();
    }
}
void Dictionary::updateStates()
{
    is_degenerate = false;
    is_basis_feasible = true;
    is_shiftable = true;
    for (int k = 0; k < B.size(); k++)
    {
        if (B[k] >= d && abs(xB(k)) < comp_zero_precision)
        {
            is_shiftable = false;
        }
        if (lb(B[k]) > std::numeric_limits<double>::lowest() && abs(xB(k) - lb(B[k])) < feas_tol)
        {
            is_degenerate = true;
        }
        if (ub(B[k]) < std::numeric_limits<double>::max() && abs(xB(k) - ub(B[k])) < feas_tol)
        {
            is_degenerate = true;
        }
        if (B[k] >= d && xB(k) < -feas_tol)
        {
            is_basis_feasible = false;
        }
    }
}
int Dictionary::recompute_decomposition(bool update)
{
    hasW = vector<bool>(d + n, false);
    eta.clear();
    has_sparse = false;
    if (sparse)
    {
        coefficients.clear();
        int i;
        for (int k = 0; k < B.size(); k++)
        {
            if (B[k] >= d && B[k] < A.cols() - 1)
            {
                i = B[k] - d;
                coefficients.push_back(Eigen::Triplet<double>(i, k, A(i, B[k])));
            }
            else
            {
                for (int i = 0; i < n; i++)
                {
                    coefficients.push_back(Eigen::Triplet<double>(i, k, A(i, B[k])));
                }
            }
        }
        sparse_AB.setFromTriplets(coefficients.begin(), coefficients.end());
        dec_sp.analyzePattern(sparse_AB);
        dec_sp.factorize(sparse_AB);
        if (dec_sp.info() != 0)
        {
            has_LU = false;
            return 0;
        }
        else
        {
            if (pricing)
            {
                has_sparse = true;
            }
            else
            {
                has_sparse = true;
            }
            if (has_sparse)
            {
                if (update)
                {
                    updateSol();
                }
            }
        }
    }
    if (!has_sparse)
    {
        AB = Eigen::MatrixXd::Zero(n, n);
        for (int j_pos = 0; j_pos < B.size(); j_pos++)
        {
            AB.col(j_pos) = A.col(B[j_pos]);
        }
        dec = Eigen::PartialPivLU<Eigen::MatrixXd>(AB);
        if (update)
        {
            updateSol();
        }
    }
    updateStates();
    has_LU = true;
    return 1;
}
void Dictionary::updateSol()
{
    Eigen::VectorXd tmp = b;
    for (int k = 0; k < N.size(); k++)
    {
        if (!is_basic[N[k]] && x(N[k]) != 0)
        {
            tmp -= A.col(N[k]) * x(N[k]);
        }
    }
    xB = solve(tmp);
    for (int k = 0; k < B.size(); k++)
    {
        x(B[k]) = xB[k];
    }
    updateStates();
}
bool Dictionary::enablePricing()
{
    pricing = true;
    if (eta.size() > 0)
    {
        int status = recompute_decomposition(true);
        if (status == 0)
        {
            return false;
        }
    }
    AN_T = Eigen::MatrixXd::Zero(A.cols() - n, n);
    for (int k = 0; k < N.size(); k++)
    {
        AN_T.row(k) = A.col(N[k]);
        cN(k) = c(N[k]);
    }
    for (int k = 0; k < N.size(); k++)
    {
        cB(k) = c(B[k]);
    }
    return true;
}
bool Dictionary::disablePricing()
{
    pricing = false;
    if (eta.size() > 0)
    {
        int status = recompute_decomposition(true);
        if (status == 0)
        {
            return false;
        }
    }
    return true;
}
Eigen::VectorXd Dictionary::reduced_costs()
{
    y = solve_T(cB);
    red_c = cN - AN_T * y;
    return red_c;
}
vector<int> Dictionary::get_Duplicates(int i)
{
    if (i > n)
    {
        return vector<int>();
    }
    else
    {
        return duplicates[i];
    }
}
int Dictionary::push_slack(int i)
{
    if (is_shiftable)
    {
        return i;
    }
    if (!hasW[d + i])
    {
        W.col(d + i) = solve(A.col(d + i));
    }
    Eigen::MatrixXd W2;
    Eigen::VectorXd w2 = W.col(d + i);
    vector<int> order2 = vector<int>(d);
    orientation_non_slack = vector<int>(d);
    for (int j = 0; j < d; j++)
    {
        order2[j] = j;
        orientation_non_slack[j] = 1;
        if (env->Rand() > 0.5)
        {
            int h = 0;
        }
    }
    shuffle(order2.begin(), order2.end(), env->generator);
    lex_order_non_slack = Permutation(order2);
    lex_order_nb.clear();
    for (int l = 0; l < lex_order.indices.size(); l++)
    {
        if (!is_basic[d + lex_order.indices[l]])
        {
            lex_order_nb.push_back(lex_order.indices[l]);
        }
    }
    lex_order_non_slack_nb.clear();
    for (int l = 0; l < lex_order_non_slack.indices.size(); l++)
    {
        if (!is_basic[lex_order_non_slack.indices[l]])
        {
            lex_order_non_slack_nb.push_back(lex_order_non_slack.indices[l]);
        }
    }
    double MAX = std::numeric_limits<double>::lowest();
    int j_pos2 = -1;
    pair<int, pair<double, int>> cand;
    bool has_cand = false;
    if (w2.isZero(piv_tol))
    {
        return i;
    }
    for (int j_pos = 0; j_pos < B.size(); j_pos++)
    {
        if (B[j_pos] < d || w2(j_pos) >= 0 || abs(w2(j_pos)) < piv_tol || xB(j_pos) > comp_zero_precision)
        {
            continue;
        }
        double tmp = xB(j_pos) / w2(j_pos);
        int cmp;
        if (!has_cand)
        {
            cmp = compare(tmp, MAX, comp_precision);
        }
        else
        {
            cmp = lex_compare_push(pair<int, pair<double, int>>(j_pos, pair<double, int>(tmp, -1)), cand, w2);
        }
        if (cmp == 1)
        {
            MAX = tmp;
            cand = pair<int, pair<double, int>>(j_pos, pair<double, int>(tmp, -1));
            has_cand = true;
        }
        else if (cmp == 0)
        {
            return -1;
        }
        else
        {
            continue;
        }
    }
    if (has_cand)
    {
        int cmp = lex_compare_push_co(i, cand, w2);
        if (cmp == -1)
        {
            return i;
        }
        else if (cmp == 1)
        {
            return B[cand.first] - d;
        }
        else
        {
            return -1;
        }
    }
    else
    {
        return i;
    }
}
void Dictionary::change_side(int i)
{
    if (i > n)
    {
        return;
    }
    vector<int> D = get_Duplicates(i);
    orientation(i) = -orientation(i);
    for (int k = 0; k < D.size(); k++)
    {
        int i2 = D[k];
        if (orientation(i) > 0)
        {
            if (!env->regression)
            {
                distribution[Y[i2]] += 1;
                distribution[n_classes + Y[i2]] -= 1;
            }
            else
            {
                reg[1].remove(Y_reg[i2]);
                reg[0].insert_inplace(Y_reg[i2]);
            }
        }
        else
        {
            if (!env->regression)
            {
                distribution[Y[i2]] -= 1;
                distribution[n_classes + Y[i2]] += 1;
            }
            else
            {
                reg[0].remove(Y_reg[i2]);
                reg[1].insert_inplace(Y_reg[i2]);
            }
        }
    }
}
bool Dictionary::is_lex_pos()
{
    return true;
    for (int k = 0; k < B.size(); k++)
    {
        if (B[k] < d || B[k] == d + n)
        {
            continue;
        }
        double val = x[B[k]];
        int i = -1;
        int l = -1;
        vector<double> vals;
        vals.push_back(val);
        if (abs(val) <= zero_precision)
        {
            for (l = 0; l < lex_order.indices.size(); l++)
            {
                i = lex_order.indices[l];
                if (!hasW[d + i])
                {
                    if (is_basic[d + i])
                    {
                        W.col(d + i) = Eigen::VectorXd::Zero(n);
                        W(position_B[d + i], d + i) = 1;
                    }
                    else
                    {
                        W.col(d + i) = solve(A.col(i + d));
                    }
                }
                vals.push_back(-W(k, d + i));
                if (abs(W(k, d + i)) > zero_precision)
                {
                    val = -W(k, d + i);
                    break;
                }
            }
        }
        if (val <= 0)
        {
            return false;
        }
    }
    return true;
}
int Dictionary::update_orientation(int i_push, int i_change)
{
    if (i_push != i_change)
    {
    }
    vector<int> order_basic;
    vector<int> order_non_basic;
    order_basic.reserve(n);
    order_non_basic.reserve(n);
    int i;
    for (int l = 0; l < lex_order.indices.size(); l++)
    {
        i = lex_order.indices[l];
        if (i == i_push || i == i_change)
        {
            continue;
        }
        if (is_basic[d + i])
        {
            order_basic.push_back(i);
        }
        else
        {
            order_non_basic.push_back(i);
        }
    }
    if (i_push == i_change)
    {
        order_non_basic.push_back(i_push);
    }
    else
    {
        order_non_basic.push_back(i_change);
        order_basic.push_back(i_push);
    }
    shuffle(order_basic.begin(), order_basic.end(), env->generator);
    order_non_basic.insert(order_non_basic.end(), order_basic.begin(), order_basic.end());
    lex_order = Permutation(order_non_basic);
    if (i_push == i_change)
    {
        change_orientation(i_change, false);
    }
    else
    {
        Eigen::VectorXd w2;
        if (!hasW[d + i_push])
        {
            w2 = solve(A.col(d + i_push));
        }
        else
        {
            w2 = W.col(d + i_push);
        }
        Eigen::VectorXd x_save(x);
        int l = d + i_change;
        int e = d + i_push;
        int l_pos = position_B[l];
        int e_pos = position_N[e];
        double t = xB(l_pos) / w2(l_pos);
        bool recomp = false;
        if (t < -feas_tol)
        {
            recomp = true;
        }
        for (int k = 0; k < B.size(); k++)
        {
            if (k != l_pos && B[k] >= d && xB(k) - t * w2(k) < -feas_tol)
            {
                recomp = true;
            }
        }
        if (recomp)
        {
            return -1;
        }
        xB(l_pos) = xN(e_pos) + t;
        x(e) = xB(l_pos);
        for (int k = 0; k < B.size(); k++)
        {
            if (k != l_pos)
            {
                xB(k) = xB(k) - t * w2(k);
                x(B[k]) = xB(k);
            }
        }
        xN(e_pos) = 0;
        x(l) = 0;
        B[l_pos] = e;
        N[e_pos] = l;
        position_B[e] = l_pos;
        position_N[l] = e_pos;
        is_basic[e] = true;
        is_basic[l] = false;
        if (pricing)
        {
            AN_T.row(e_pos) = A.col(l);
            double tmp = cB(l_pos);
            cB(l_pos) = cN(e_pos);
            cN(e_pos) = tmp;
        }
        hasW = vector<bool>(d + n, false);
        eta.push_back(Eta(w2, l_pos));
        change_orientation(i_change, false);
        updateStates();
    }
    return 1;
}
void Dictionary::change_orientation(int i, bool revert)
{
    if (i > n)
    {
        return;
    }
    if (hasW[d + i])
    {
        W.col(d + i) = -W.col(d + i);
    }
    has_shifted = false;
    orientation(i) = -orientation(i);
    A(i, A.cols() - 1) = orientation(i);
    A(i, d + i) = orientation(i);
    if (is_basic[d + i])
    {
        if (revert)
        {
            eta.pop_back();
            is_redundant = is_redundant_save;
            hasW = hasW_save;
        }
        else
        {
            is_redundant_save = is_redundant;
            is_redundant = vector<bool>(A.rows(), false);
            eta.push_back(Eta(position_B[d + i]));
            hasW_save = hasW;
        }
    }
    else
    {
        if (revert)
        {
            is_redundant = is_redundant_save;
            hasW = hasW_save;
        }
        else
        {
            is_redundant_save = is_redundant;
            is_redundant = vector<bool>(A.rows(), false);
            hasW_save = hasW;
        }
        if (pricing)
        {
            AN_T(position_N[d + i], i) = orientation(i);
        }
    }
    vector<int> D = get_Duplicates(i);
    for (int k = 0; k < D.size(); k++)
    {
        int i2 = D[k];
        if (orientation(i) > 0)
        {
            if (!env->regression)
            {
                distribution[Y[i2]] += 1;
                distribution[n_classes + Y[i2]] -= 1;
            }
            else
            {
                reg[1].remove(Y_reg[i2]);
                reg[0].insert_inplace(Y_reg[i2]);
            }
        }
        else
        {
            if (!env->regression)
            {
                distribution[Y[i2]] -= 1;
                distribution[n_classes + Y[i2]] += 1;
            }
            else
            {
                reg[0].remove(Y_reg[i2]);
                reg[1].insert_inplace(Y_reg[i2]);
            }
        }
    }
}
bool approximatelyEqual(float a, float b, float epsilon)
{
    return fabs(a - b) <= ((fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}
bool essentiallyEqual(float a, float b, float epsilon)
{
    return fabs(a - b) <= ((fabs(a) > fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}
bool definitelyGreaterThan(float a, float b, float epsilon)
{
    return (a - b) > ((fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}
bool definitelyLessThan(float a, float b, float epsilon)
{
    return (b - a) > ((fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}
int Dictionary::compare(double a, double b, double prec)
{
    if (a == std::numeric_limits<double>::max() && b == std::numeric_limits<double>::max())
    {
        return 0;
    }
    if (a == std::numeric_limits<double>::lowest() && b == std::numeric_limits<double>::lowest())
    {
        return 0;
    }
    if (a == std::numeric_limits<double>::max())
    {
        return 1;
    }
    if (a == std::numeric_limits<double>::lowest())
    {
        return -1;
    }
    if (b == std::numeric_limits<double>::max())
    {
        return -1;
    }
    if (b == std::numeric_limits<double>::lowest())
    {
        return 1;
    }
    if (definitelyLessThan(a, b, prec))
    {
        return -1;
    }
    else if (definitelyGreaterThan(a, b, prec))
    {
        return 1;
    }
    else
    {
        return 0;
    }
}
int Dictionary::lex_compare_push_co(int i_push, const pair<int, pair<double, int>> &cand2, Eigen::VectorXd &w2)
{
    int j_pos2 = cand2.first;
    double v1 = 0;
    double v2 = xB(j_pos2);
    if (abs(v1) > comp_zero_precision || abs(v2) > comp_zero_precision)
    {
        int cmp = compare(v1, v2, comp_precision);
        if (cmp != 0)
        {
            return cmp;
        }
    }
    else
    {
        for (int l = 0; l < lex_order_nb.size(); l++)
        {
            int i = lex_order_nb[l];
            if (is_basic[d + i])
            {
                continue;
            }
            if (!hasW[d + i])
            {
                if (is_basic[d + i])
                {
                    W.col(d + i) = Eigen::VectorXd::Zero(n);
                    W(position_B[d + i], d + i) = 1;
                }
                else
                {
                    W.col(d + i) = solve(A.col(d + i));
                }
            }
            hasW[d + i] = true;
            double v1 = 0;
            if (i == i_push)
            {
                v1 = -w2(j_pos2);
            }
            double v2 = -W(j_pos2, d + i);
            if (abs(v1) > comp_zero_precision || abs(v2) > comp_zero_precision)
            {
                int cmp = compare(v1, v2, comp_precision);
                if (cmp != 0)
                {
                    return cmp;
                }
            }
        }
    }
    return 0;
}
void Dictionary::shift()
{
    x_shifted = x;
    double tmp;
    for (int l = 0; l < lex_order.indices.size(); l++)
    {
        double t = std::numeric_limits<double>::max();
        int i = lex_order.indices[l];
        if (is_basic[d + i])
        {
            continue;
        }
        if (!hasW[d + i])
        {
            w = solve(A.col(d + i));
        }
        else
        {
            w = W.col(d + i);
        }
        for (int j_pos = 0; j_pos < B.size(); j_pos++)
        {
            if (B[j_pos] < d)
            {
                continue;
            }
            if (abs(w(j_pos)) < piv_tol || w(j_pos) < 0)
            {
                continue;
            }
            tmp = x_shifted(B[j_pos]) / w(j_pos);
            if (tmp < t)
            {
                t = tmp;
            }
        }
        t = t / 2.;
        t = min(1., t);
        x_shifted(d + i) = t;
        for (int j_pos = 0; j_pos < B.size(); j_pos++)
        {
            x_shifted(B[j_pos]) -= t * w(j_pos);
        }
    }
    has_shifted = true;
}
int Dictionary::lex_compare_push(const pair<int, pair<double, int>> &cand1, const pair<int, pair<double, int>> &cand2, Eigen::VectorXd &w2)
{
    int j_pos1 = cand1.first;
    int j_pos2 = cand2.first;
    double v1 = cand1.second.first;
    double v2 = cand2.second.first;
    int bound1 = cand1.second.second;
    int bound2 = cand2.second.second;
    v1 = xB(j_pos1);
    v2 = xB(j_pos2);
    if (abs(xB(j_pos1)) > comp_zero_precision || abs(xB(j_pos2)) > comp_zero_precision)
    {
        int cmp = compare(v1, v2, comp_precision);
        if (cmp != 0)
        {
            return cmp;
        }
    }
    for (int l = 0; l < lex_order_nb.size(); l++)
    {
        int i = lex_order_nb[l];
        if (is_basic[d + i])
        {
            continue;
        }
        if (!hasW[d + i])
        {
            if (is_basic[d + i])
            {
                W.col(d + i) = Eigen::VectorXd::Zero(n);
                W(position_B[d + i], d + i) = 1;
            }
            else
            {
                W.col(d + i) = solve(A.col(d + i));
            }
        }
        hasW[d + i] = true;
        double v1 = W(j_pos1, d + i);
        double v2 = W(j_pos2, d + i);
        if (abs(v1) > comp_zero_precision || abs(v2) > comp_zero_precision)
        {
            v1 = -v1 / w2(j_pos1);
            v2 = -v2 / w2(j_pos2);
            int cmp = compare(v1, v2, comp_precision);
            if (cmp != 0)
            {
                return cmp;
            }
        }
    }
    return 0;
}
int Dictionary::lex_compare(const pair<int, pair<double, int>> &cand1, const pair<int, pair<double, int>> &cand2)
{
    int j_pos1 = cand1.first;
    int j_pos2 = cand2.first;
    double v1 = cand1.second.first;
    double v2 = cand2.second.first;
    int bound1 = cand1.second.second;
    int bound2 = cand2.second.second;
    if (abs(xB(j_pos1)) > comp_zero_precision || abs(xB(j_pos2)) > comp_zero_precision)
    {
        int cmp = compare(v1, v2, comp_precision);
        if (cmp != 0)
        {
            return cmp;
        }
    }
    for (int l = 0; l < lex_order.indices.size(); l++)
    {
        int i = lex_order.indices[l];
        if (!hasW[d + i])
        {
            if (is_basic[d + i])
            {
                W.col(d + i) = Eigen::VectorXd::Zero(n);
                W(position_B[d + i], d + i) = 1;
            }
            else
            {
                W.col(d + i) = solve(A.col(d + i));
            }
        }
        hasW[d + i] = true;
        double v1 = W(j_pos1, d + i);
        double v2 = W(j_pos2, d + i);
        if (abs(v1) > comp_zero_precision || abs(v2) > comp_zero_precision)
        {
            v1 = -v1 / w(j_pos1);
            v2 = -v2 / w(j_pos2);
            int cmp = compare(v1, v2, comp_precision);
            if (cmp != 0)
            {
                return cmp;
            }
        }
    }
    return 0;
}
vector<pair<int, pair<double, int>>> Dictionary::find_pivot(int e, int type, int sign)
{
    double upperbound = std::numeric_limits<double>::max();
    double lowerbound = std::numeric_limits<double>::lowest();
    double MAX = std::numeric_limits<double>::max();
    double LOWEST = std::numeric_limits<double>::lowest();
    vector<pair<int, pair<double, int>>> cands_u;
    vector<pair<int, pair<double, int>>> cands_l;
    cands_u.reserve(B.size());
    cands_l.reserve(B.size());
    double eps = 0;
    double prec = comp_precision;
    if (w.isZero(piv_tol))
    {
        return cands_u;
    }
    for (int j_pos = 0; j_pos < B.size(); j_pos++)
    {
        if (abs(w(j_pos)) < piv_tol)
        {
            continue;
        }
        else
        {
            double tmp;
            if (lb(B[j_pos]) > LOWEST)
            {
                tmp = (xB(j_pos) - (lb(B[j_pos]) - eps)) / w(j_pos);
                if (w(j_pos) > 0)
                {
                    if (sign >= 0)
                    {
                        int comp = compare(tmp, upperbound, prec);
                        if (cands_u.size() > 0 && use_lex)
                        {
                            comp = lex_compare(pair<int, pair<double, int>>(j_pos, pair<double, int>(tmp, -1)), cands_u[0]);
                        }
                        if (comp == 0)
                        {
                            cands_u.push_back(pair<int, pair<double, int>>(j_pos, pair<double, int>(tmp, -1)));
                        }
                        else if (comp == -1)
                        {
                            cands_u.clear();
                            cands_u.push_back(pair<int, pair<double, int>>(j_pos, pair<double, int>(tmp, -1)));
                            upperbound = tmp;
                        }
                    }
                }
                else
                {
                    if (sign <= 0)
                    {
                        int comp = compare(tmp, lowerbound, prec);
                        if (cands_l.size() > 0 && use_lex)
                        {
                            comp = lex_compare(pair<int, pair<double, int>>(j_pos, pair<double, int>(tmp, -1)), cands_l[0]);
                        }
                        if (comp == 0)
                        {
                            cands_l.push_back(pair<int, pair<double, int>>(j_pos, pair<double, int>(tmp, -1)));
                        }
                        else if (comp == 1)
                        {
                            cands_l.clear();
                            cands_l.push_back(pair<int, pair<double, int>>(j_pos, pair<double, int>(tmp, -1)));
                            lowerbound = tmp;
                        }
                    }
                }
            }
            if (ub(B[j_pos]) < MAX)
            {
                tmp = (xB(j_pos) - (ub(B[j_pos]) + eps)) / w(j_pos);
                if (w(j_pos) > 0)
                {
                    if (sign <= 0)
                    {
                        int comp = compare(tmp, lowerbound, prec);
                        if (cands_l.size() > 0 && use_lex)
                        {
                            comp = lex_compare(pair<int, pair<double, int>>(j_pos, pair<double, int>(tmp, 1)), cands_l[0]);
                        }
                        if (comp == 0)
                        {
                            cands_l.push_back(pair<int, pair<double, int>>(j_pos, pair<double, int>(tmp, 1)));
                        }
                        else if (comp == 1)
                        {
                            cands_l.clear();
                            cands_l.push_back(pair<int, pair<double, int>>(j_pos, pair<double, int>(tmp, 1)));
                            lowerbound = tmp;
                        }
                    }
                }
                else
                {
                    if (sign >= 0)
                    {
                        int comp = compare(tmp, upperbound, prec);
                        if (cands_u.size() > 0 && use_lex)
                        {
                            comp = lex_compare(pair<int, pair<double, int>>(j_pos, pair<double, int>(tmp, 1)), cands_u[0]);
                        }
                        if (comp == 0)
                        {
                            cands_u.push_back(pair<int, pair<double, int>>(j_pos, pair<double, int>(tmp, 1)));
                        }
                        else if (comp == -1)
                        {
                            cands_u.clear();
                            cands_u.push_back(pair<int, pair<double, int>>(j_pos, pair<double, int>(tmp, 1)));
                            upperbound = tmp;
                        }
                    }
                }
            }
        }
    }
    if (sign == 1)
    {
        if (cands_u.size() > 0)
        {
        }
        return cands_u;
    }
    else if (sign == -1)
    {
        return cands_l;
    }
    else
    {
        for (int k = 0; k < cands_l.size(); k++)
        {
            cands_u.push_back(cands_l[k]);
        }
        return cands_u;
    }
}
int Dictionary::pivot(int e_pos, int type, int sign, bool use_order)
{
    vector<int> save_B(B);
    Eigen::VectorXd save_x(x);
    int e = N[e_pos];
    vector<pair<int, pair<double, int>>> cands;
    cands.reserve(B.size());
    w = solve(A.col(e));
    if (e < d + n)
    {
        W.col(e) = w;
        hasW[e] = true;
    }
    if (type == 1 && !is_basic[A.cols() - 1] && N[e_pos] == A.cols() - 1)
    {
        double min = std::numeric_limits<double>::max();
        for (int k = 0; k < B.size(); k++)
        {
            if (xB(k) >= 0)
            {
                continue;
            }
            if (xB(k) == min && B[k] >= d)
            {
                cands.push_back(pair<int, pair<double, int>>(k, pair<double, int>(-xB(k), -1)));
            }
            else if (xB(k) < min && B[k] >= d)
            {
                min = xB(k);
                cands.clear();
                cands.push_back(pair<int, pair<double, int>>(k, pair<double, int>(-xB(k), -1)));
            }
        }
        if (cands.size() == 0)
        {
            return -1;
        }
        if (use_lex)
        {
            pair<int, pair<double, int>> cand = cands[0];
            int i2 = B[cand.first] - d;
            for (int k = 1; k < cands.size(); k++)
            {
                int i = B[cands[k].first] - d;
                if (lex_order.indexOf[i] < lex_order.indexOf[i2])
                {
                    cand = cands[k];
                    i2 = i;
                }
            }
            cands.clear();
            cands.push_back(cand);
        }
    }
    else
    {
        cands = find_pivot(e, type, sign);
    }
    int l_pos = -1;
    double t;
    int lOu;
    if ((cands.size() > 1 && e >= d && use_lex) || (cands.size() > 2 && e < d && use_lex))
    {
        if (eta.size() > 0)
        {
            return -2;
        }
        else
        {
            return -3;
        }
    }
    if (cands.size() > 0)
    {
        Eigen::VectorXd x_old = x;
        if (type == 0)
        {
            int rand_ind = env->rand() % cands.size();
            l_pos = cands[rand_ind].first;
            t = cands[rand_ind].second.first;
            lOu = cands[rand_ind].second.second;
        }
        else
        {
            bool flag = false;
            if (type == 1)
            {
                if (is_basic[A.cols() - 1])
                {
                    for (int k = 0; k < cands.size(); k++)
                    {
                        if (B[cands[k].first] == A.cols() - 1)
                        {
                            l_pos = cands[k].first;
                            t = cands[k].second.first;
                            lOu = cands[k].second.second;
                            flag = true;
                        }
                    }
                }
            }
            if (!flag)
            {
                int rand_ind = env->rand() % cands.size();
                l_pos = cands[rand_ind].first;
                t = cands[rand_ind].second.first;
                lOu = cands[rand_ind].second.second;
                if (bland)
                {
                    for (int k = 0; k < cands.size(); k++)
                    {
                        if (B[cands[k].first] < B[l_pos])
                        {
                            l_pos = cands[k].first;
                            t = cands[k].second.first;
                            lOu = cands[k].second.second;
                        }
                    }
                }
            }
        }
        int l = B[l_pos];
        double save = xB(l_pos);
        double save2 = x(B[l_pos]);
        updateStates();
        if (!is_basis_feasible && !type == 1)
        {
            recompute_decomposition(true);
            return -3;
        }
        bool recomp = false;
        if (e >= d && x(e) + t < lb(e) - feas_tol)
        {
            recomp = true;
        }
        for (int k = 0; k < B.size(); k++)
        {
            if (k != l_pos && B[k] >= d && xB(k) - t * w(k) < -feas_tol)
            {
                recomp = true;
            }
        }
        if (type == 1)
        {
            recomp = false;
        }
        if (recomp)
        {
            return -2;
        }
        for (int k = 0; k < B.size(); k++)
        {
            if (k != l_pos)
            {
                xB(k) = xB(k) - t * w(k);
                x(B[k]) = xB(k);
            }
        }
        xB(l_pos) = xN(e_pos) + t;
        x(e) = xB(l_pos);
        if (lOu == -1)
        {
            x(l) = lb(B[l_pos]);
        }
        if (lOu == 0)
        {
            x(l) = 0;
        }
        else if (lOu == 1)
        {
            x(l) = ub(B[l_pos]);
        }
        xN(e_pos) = x(l);
        B[l_pos] = e;
        N[e_pos] = l;
        position_B[e] = l_pos;
        position_N[l] = e_pos;
        is_basic[e] = true;
        is_basic[l] = false;
        hasW = vector<bool>(d + n, false);
        if (pricing)
        {
            AN_T.row(e_pos) = A.col(l);
            double tmp = cB(l_pos);
            cB(l_pos) = cN(e_pos);
            cN(e_pos) = tmp;
        }
        updateStates();
        if (!is_basis_feasible)
        {
            if (eta.size() == 0)
            {
                setBasis(save_B, save_x, Permutation(lex_order));
                updateStates();
                return -1;
            }
            else
            {
                setBasis(save_B, save_x, Permutation(lex_order));
                e_pos = -1;
                for (int k = 0; k < N.size(); k++)
                {
                    if (N[k] == e)
                    {
                        e_pos = k;
                    }
                }
                return pivot(e_pos, type, sign, use_order);
            }
        }
        vector<int> order_basic;
        vector<int> order_non_basic;
        order_basic.reserve(n);
        order_non_basic.reserve(n);
        int i;
        for (int l = 0; l < lex_order.indices.size(); l++)
        {
            i = lex_order.indices[l];
            if (is_basic[d + i])
            {
                order_basic.push_back(i);
            }
            else
            {
                order_non_basic.push_back(i);
            }
        }
        shuffle(order_basic.begin(), order_basic.end(), env->generator);
        order_non_basic.insert(order_non_basic.end(), order_basic.begin(), order_basic.end());
        lex_order = Permutation(order_non_basic);
        eta.push_back(Eta(w, l_pos));
        return 1;
    }
    return 0;
}
int Dictionary::get_cobasic_var_rand(int method, bool exclude_last)
{
    int e_pos = env->rand() % (N.size());
    while (exclude_last && N[e_pos] == A.cols() - 1)
    {
        e_pos = env->rand() % (N.size());
    }
    return e_pos;
}
int Dictionary::get_cobasic_var_price(int method, bool exclude_last)
{
    reduced_costs();
    vector<int> cands;
    bool opt = true;
    for (int k = 0; k < N.size(); k++)
    {
        if (N[k] == A.cols() - 1 && exclude_last)
        {
            continue;
        }
        if (red_c(k) >= dual_feas_tol && x(N[k]) < ub(N[k]) - zero_precision)
        {
            opt = false;
            cands.push_back(k);
        }
        if (red_c(k) <= -dual_feas_tol && x(N[k]) > lb(N[k]) + zero_precision)
        {
            opt = false;
            cands.push_back(k);
        }
    }
    if (opt)
    {
        return -1;
    }
    if (cands.size() == 0)
    {
        return -1;
    }
    if (method == 0)
    {
        return cands[env->rand() % (cands.size())];
    }
    else if (method > 0)
    {
        int e_pos = cands[env->rand() % (cands.size())];
        if (bland)
        {
            for (int k = 0; k < cands.size(); k++)
            {
                if (N[cands[k]] < N[e_pos])
                {
                    e_pos = cands[k];
                }
            }
        }
        return e_pos;
    }
    return cands[0];
}
int Dictionary::get_zero_slack(int method)
{
    vector<int> N2;
    if (!is_shiftable)
    {
        for (int k = 0; k < B.size(); k++)
        {
            if (B[k] >= d && B[k] < A.cols() - 1 && abs(x(B[k])) < deg_tol && !is_redundant[B[k] - d])
            {
                N2.push_back(B[k]);
            }
        }
    }
    for (int k = 0; k < N.size(); k++)
    {
        if (N[k] >= d && N[k] < A.cols() - 1 && !is_redundant[N[k] - d])
        {
            N2.push_back(N[k]);
        }
    }
    if (N2.size() == 0)
    {
        return -1;
    }
    return N2[env->rand() % N2.size()] - d;
}
bool compareMoves(const pair<pair<int, int>, double> &a, const pair<pair<int, int>, double> &b)
{
    return a.second < b.second;
}
int Dictionary::get_cobasic_slack(int method)
{
    if (method == 0)
    {
        vector<int> N2;
        for (int k = 0; k < N.size(); k++)
        {
            if (N[k] >= d && N[k] < A.cols() - 1 && !is_redundant[N[k] - d])
            {
                N2.push_back(N[k]);
            }
        }
        if (N2.size() == 0)
        {
            return -1;
        }
        return N2[env->rand() % N2.size()] - d;
    }
    else
    {
        vector<pair<pair<int, int>, double>> cands;
        for (int k = 0; k < N.size(); k++)
        {
            if (N[k] < d || N[k] >= A.cols() - 1)
            {
                continue;
            }
            int i_push2 = N[k] - d;
            int i_change2 = push_slack(i_push2);
            change_side(i_change2);
            double imp_t;
            ;
            if (!env->regression)
            {
                imp_t = imp(distribution, n_classes, env->criterion);
            }
            else
            {
                imp_t = reg[0].error(env->criterion) + reg[1].error(env->criterion);
            }
            cands.push_back(pair<pair<int, int>, double>(pair<int, int>(i_push2, i_change2), imp_t));
            change_side(i_change2);
        }
        std::shuffle(cands.begin(), cands.end(), env->generator);
        std::sort(cands.begin(), cands.end(), compareMoves);
        if (method == 2)
        {
            return cands[0].first.first;
        }
        vector<double> ranks(cands.size());
        double sum = 0;
        for (int k = 0; k < cands.size(); k++)
        {
            ranks[k] = cands.size() - k;
            sum += ranks[k];
        }
        sum = 0;
        double p = env->Rand();
        int i_push = -1;
        for (int k = 0; k < cands.size(); k++)
        {
            i_push = cands[k].first.first;
            if (sum + ranks[k] > p)
            {
                return i_push;
            }
            sum += ranks[k];
        }
        return cands[0].first.first;
    }
}
Eigen::VectorXd Dictionary::solve(Eigen::VectorXd rhs)
{
    Eigen::VectorXd w = rhs;
    if (has_sparse)
    {
        w = dec_sp.solve(w);
    }
    else
    {
        w = dec.solve(w);
    }
    for (int k = 0; k < eta.size(); k++)
    {
        eta[k].solve(w);
    }
    return w;
}
Eigen::VectorXd Dictionary::solve_T(Eigen::VectorXd lhs)
{
    Eigen::VectorXd y = lhs;
    for (int k = eta.size() - 1; k >= 0; k--)
    {
        eta[k].solve_transpose(y);
    }
    if (has_sparse)
    {
        y = dec_sp.transpose().solve(y);
    }
    else
    {
        y = dec.transpose().solve(y);
    }
    return y;
}
bool Dictionary::is_inner(Eigen::VectorXd x_tmp)
{
    has_shifted = true;
    for (int i = 0; i < n; i++)
    {
        double sum = 0;
        for (int j = 0; j < d; j++)
        {
            sum += x_tmp(j) * A(i, j);
        }
        if (orientation(i) > 0 && sum >= b(i) - opt_requirement)
        {
            has_shifted = false;
        }
        if (orientation(i) < 0 && sum <= b(i) + opt_requirement)
        {
            has_shifted = false;
        }
    }
    if (has_shifted)
    {
        x_shifted = x_tmp;
    }
    return has_shifted;
}
void Dictionary::setOrientation(Eigen::VectorXd &orientation_tmp)
{
    reset_Redundancy();
    orientation = orientation_tmp;
    for (int i = 0; i < orientation.size(); i++)
    {
        orientation(i) = orientation_tmp(i);
    }
    for (int i = 0; i < n; i++)
    {
        A(i, d + i) = orientation(i);
    }
    init_Bounds();
    if (!env->regression)
    {
        distribution = vector<double>(2 * n_classes, 0);
        for (int k = 0; k < 2 * n_classes; k++)
        {
            distribution[k] = 0;
        }
    }
    else
    {
        reg[0].clear();
        reg[1].clear();
    }
    if (!env->regression)
    {
        for (int i = 0; i < n; i++)
        {
            if (orientation[i] > 0)
            {
                for (int k = 0; k < duplicates[i].size(); k++)
                {
                    distribution[Y[duplicates[i][k]]] += 1;
                }
            }
            else
            {
                for (int k = 0; k < duplicates[i].size(); k++)
                {
                    distribution[n_classes + Y[duplicates[i][k]]] += 1;
                }
            }
        }
    }
    else
    {
        for (int i = 0; i < n; i++)
        {
            if (orientation[i] > 0)
            {
                for (int k = 0; k < duplicates[i].size(); k++)
                {
                    reg[0].push_back(Y_reg[duplicates[i][k]]);
                }
            }
            else
            {
                for (int k = 0; k < duplicates[i].size(); k++)
                {
                    reg[1].push_back(Y_reg[duplicates[i][k]]);
                }
            }
        }
        reg[0].sort();
        reg[1].sort();
    }
}
void Dictionary::reset_Redundancy()
{
    is_redundant = vector<bool>(A.rows(), false);
    is_redundant_save = vector<bool>(A.rows(), false);
}
void Dictionary::preparePhaseI(bool feasibility_check)
{
    init_Bounds();
    Eigen::VectorXd tmp = b;
    for (int i = 0; i < n; i++)
    {
        if (orientation(i) * tmp(i) < 0)
        {
            A(i, A.cols() - 1) = -orientation(i);
        }
        else
        {
            A(i, A.cols() - 1) = 0;
        }
    }
    c = Eigen::VectorXd::Zero(A.cols());
    c(A.cols() - 1) = -1;
    cB = Eigen::VectorXd::Zero(B.size());
    cN = Eigen::VectorXd::Zero(N.size());
    cN(position_N[A.cols() - 1]) = -1;
    if (feasibility_check)
    {
        B.clear();
        N.clear();
        for (int i = 0; i < d + n + 1; i++)
        {
            if (i < d || i == d + n)
            {
                N.push_back(i);
            }
            else
            {
                B.push_back(i);
            }
        }
        x = Eigen::VectorXd::Zero(A.cols());
        setBasis(B, x, lex_order);
    }
}
void Dictionary::prepareFeasibilityCheck()
{
    init_Bounds();
    for (int i = 0; i < n; i++)
    {
        A(i, A.cols() - 1) = orientation(i);
    }
    c(A.cols() - 1) = 1;
    cB = Eigen::VectorXd::Zero(B.size());
    cN = Eigen::VectorXd::Zero(N.size());
    if (!is_basic[A.cols() - 1])
    {
        cN(position_N[A.cols() - 1]) = 1;
    }
    else
    {
        cB(position_B[A.cols() - 1]) = 1;
    }
}
Eigen::VectorXd Dictionary::getUnscaled()
{
    Eigen::VectorXd a_min(d + 1);
    for (int j = 0; j < d; j++)
    {
        a_min(j) = x_shifted(j) / col_scale[j];
    }
    a_min(d) = -1;
    return a_min;
}
double Dictionary::getResidual()
{
    if (!env->check_residual)
    {
        return 0;
    }
    Eigen::VectorXd prod = Eigen::VectorXd::Zero(b.size());
    for (int k = 0; k < B.size(); k++)
    {
        if (B[k] < d || B[k] == A.cols() - 1)
        {
            prod += x(B[k]) * A.col(B[k]);
        }
        else
        {
            int i2 = B[k] - d;
            prod(i2) += orientation(i2) * x(d + i2);
        }
    }
    return (prod - b).lpNorm<Eigen::Infinity>();
}