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

#define EIGEN_NO_DEBUG
#include "Heuristics.h"
#include "Dictionary.h"
#include "Environment.h"
#include "Eigen/Sparse"
#include "vMF.h"
#include "General.h"
#include "Normalize.h"
/*#ifdef use_gurobi
#include "gurobi_c++.h"
#endif*/
#ifndef use_gurobi
#include "ClpSimplex.hpp"
#include "CoinBuild.hpp"
#endif
#include "linear.h"
using namespace std;
typedef std::chrono::high_resolution_clock Clock;
double evaluate(Eigen::VectorXd &a, double **X, int *Y, double *Y_reg, int total_points, int n_classes, bool reg, int criterion)
{
    if (!reg)
    {
        vector<double> distr(2 * n_classes, 0);
        int d = a.size();
        for (int i = 0; i < total_points; i++)
        {
            double sum = 0;
            for (int j = 0; j < d; j++)
            {
                sum += a(j) * X[i][j];
            }
            if (sum >= 0)
            {
                distr[n_classes + Y[i]] += 1;
            }
            else
            {
                distr[Y[i]] += 1;
            }
        }
        return imp(distr, n_classes, criterion);
    }
    else
    {
        Vector_reg left;
        Vector_reg right;
        if (criterion == Environment::mse)
        {
            left.keep_vec = false;
            right.keep_vec = false;
        }
        else
        {
            left.use_mae();
            right.use_mae();
        }
        left.reserve(total_points);
        right.reserve(total_points);
        int d = a.size();
        double sumL = 0;
        double sumL2 = 0;
        double n_L = 0;
        double sumR = 0;
        double sumR2 = 0;
        double n_R = 0;
        for (int i = 0; i < total_points; i++)
        {
            double sum = 0;
            for (int j = 0; j < d; j++)
            {
                sum += a(j) * X[i][j];
            }
            if (sum >= 0)
            {
                sumL += Y_reg[i];
                sumL2 += Y_reg[i] * Y_reg[i];
                n_L += 1;
                left.push_back(Y_reg[i]);
            }
            else
            {
                sumR += Y_reg[i];
                sumR2 += Y_reg[i] * Y_reg[i];
                n_R += 1;
                right.push_back(Y_reg[i]);
            }
        }
        double error = 0;
        if (n_L > 0)
        {
            error = left.error(criterion);
            if (n_R > 0)
            {
                error += right.error(criterion);
            }
        }
        else if (n_R > 0)
        {
            error = right.error(criterion);
        }
        else
        {
            error = std::numeric_limits<double>::max();
        }
        return error;
    }
}
double evaluate_p(int id, Eigen::VectorXd &a, double **X, int *Y, double *Y_reg, int total_points, int n_classes, bool reg, int criterion)
{
    return evaluate(a, X, Y, Y_reg, total_points, n_classes, reg, criterion);
}
bool metropolis(double old_imp, double new_imp, double t, double rand, bool preAccept)
{
    if (preAccept || new_imp <= old_imp)
    {
        return true;
    }
    else
    {
        double alpha = exp(-(new_imp - old_imp) / t);
        if (rand <= alpha)
        {
            return true;
        }
    }
    return false;
}
bool PhaseI(Dictionary &dict, bool feasibility_check)
{
    dict.preparePhaseI(feasibility_check);
    if (!dict.enablePricing())
    {
        dict.disablePricing();
        return false;
    }
    int e_pos = dict.position_N[dict.A.cols() - 1];
    int status = dict.pivot(e_pos, 1, 1, false);
    if (status == -1)
    {
        dict.disablePricing();
        return true;
    }
    vector<int> B_recov = vector<int>(dict.B);
    VectorXd orientation_recov = dict.orientation;
    bool has_shifted_recov = dict.has_shifted;
    VectorXd x_recov = dict.x;
    VectorXd x_recov_shifted = dict.x_shifted;
    Permutation lex_order_recov = Permutation(dict.lex_order);
    int it = 0;
    while (true)
    {
        it++;
        if (it % 10000 == 0)
        {
            dict.disablePricing();
            return false;
        }
        if (!dict.is_basic[dict.A.cols() - 1])
        {
            if (!dict.disablePricing())
            {
                return false;
            }
            return true;
        }
        bool is_opt = true;
        for (int i = 0; i < dict.n; i++)
        {
            if (abs(dict.x(dict.A.cols() - 1) + dict.x(dict.d + i)) < dict.opt_requirement)
            {
                is_opt = false;
            }
        }
        double residual = 0;
        residual = dict.getResidual();
        if (dict.eta.size() > dict.max_eta || residual > dict.res_tol)
        {
            dict.recompute_decomposition(true);
            if (!dict.has_LU)
            {
                dict.setOrientation(orientation_recov);
                dict.setBasis(B_recov, x_recov, lex_order_recov);
                if (!dict.has_LU)
                {
                    dict.disablePricing();
                    return false;
                }
            }
            else
            {
                dict.updateStates();
                if (dict.is_basis_feasible)
                {
                    B_recov = vector<int>(dict.B);
                    orientation_recov = dict.orientation;
                    has_shifted_recov = dict.has_shifted;
                    x_recov = dict.x;
                    x_recov_shifted = dict.x_shifted;
                    lex_order_recov = Permutation(dict.lex_order);
                }
            }
        }
        e_pos = dict.get_cobasic_var_price(1, false);
        if (e_pos < 0)
        {
            for (int i = 0; i < dict.n; i++)
            {
                double min = dict.x(dict.A.cols() - 1);
                if (abs(dict.x(dict.A.cols() - 1) + dict.x(dict.d + i)) < dict.opt_requirement)
                {
                    is_opt = false;
                }
            }
            dict.disablePricing();
            return false;
        }
        int sign = 1;
        if (dict.red_c(e_pos) < 0 && dict.x(dict.N[e_pos]) > dict.lb_fixed(dict.N[e_pos]))
        {
            sign = -1;
        }
        status = dict.pivot(e_pos, 0, sign, false);
        if (status == -2)
        {
            dict.recompute_decomposition(true);
            if (!dict.has_LU)
            {
                dict.setOrientation(orientation_recov);
                dict.setBasis(B_recov, x_recov, lex_order_recov);
                if (!dict.has_LU)
                {
                    dict.disablePricing();
                    return false;
                }
            }
        }
        if (status == -3)
        {
            dict.setOrientation(orientation_recov);
            dict.setBasis(B_recov, x_recov, lex_order_recov);
            if (!dict.has_LU)
            {
                dict.disablePricing();
                return false;
            }
        }
        if (status == 0)
        {
            dict.disablePricing();
            return false;
        }
    }
}
bool isCell(Dictionary &dict)
{
    int n = dict.n;
    int d = dict.d;
    dict.prepareFeasibilityCheck();
    if (!dict.enablePricing())
    {
        dict.disablePricing();
        return false;
    }
    int e_pos = dict.position_N[dict.A.cols() - 1];
    int status;
    vector<int> B_recov = vector<int>(dict.B);
    VectorXd orientation_recov = dict.orientation;
    bool has_shifted_recov = dict.has_shifted;
    VectorXd x_recov = dict.x;
    VectorXd x_recov_shifted = dict.x_shifted;
    Permutation lex_order_recov = Permutation(dict.lex_order);
    int it = 0;
    while (true)
    {
        it++;
        if (it % 10000 == 0)
        {
            dict.disablePricing();
            return false;
        }
        double residual = 0;
        residual = dict.getResidual();
        if (dict.eta.size() > dict.max_eta || residual > dict.res_tol)
        {
            dict.recompute_decomposition(true);
            if (!dict.has_LU)
            {
                dict.setOrientation(orientation_recov);
                dict.setBasis(B_recov, x_recov, lex_order_recov);
                if (!dict.has_LU)
                {
                    dict.disablePricing();
                    return false;
                }
            }
            if (dict.is_basis_feasible)
            {
                B_recov = vector<int>(dict.B);
                orientation_recov = dict.orientation;
                has_shifted_recov = dict.has_shifted;
                x_recov = dict.x;
                x_recov_shifted = dict.x_shifted;
                lex_order_recov = Permutation(dict.lex_order);
            }
        }
        e_pos = dict.get_cobasic_var_price(3, false);
        if (e_pos < 0 || dict.x(dict.A.cols() - 1) > 1)
        {
            dict.disablePricing();
            double min = dict.x(dict.A.cols() - 1);
            if (min > dict.opt_requirement)
            {
                dict.x_shifted = dict.x;
                dict.has_shifted = true;
                return true;
            }
            for (int i = 0; i < dict.n; i++)
            {
                if (dict.x(dict.A.cols() - 1) + dict.x(dict.d + i) < min)
                {
                    min = dict.x(dict.A.cols() - 1) + dict.x(dict.d + i);
                }
            }
            return false;
        }
        int sign = 1;
        if (dict.red_c(e_pos) < 0 && dict.x(dict.N[e_pos]) > dict.lb(dict.N[e_pos]))
        {
            sign = -1;
        }
        double save = dict.xN(e_pos);
        status = dict.pivot(e_pos, 0, sign, false);
        if (status == -2)
        {
            dict.recompute_decomposition(true);
            if (!dict.has_LU)
            {
                dict.setOrientation(orientation_recov);
                dict.setBasis(B_recov, x_recov, lex_order_recov);
                if (!dict.has_LU)
                {
                    dict.disablePricing();
                    return false;
                }
            }
        }
        if (status == -3)
        {
            dict.setOrientation(orientation_recov);
            dict.setBasis(B_recov, x_recov, lex_order_recov);
            if (!dict.has_LU)
            {
                dict.disablePricing();
                return false;
            }
        }
        if (status == 0)
        {
            double obj = dict.x(dict.A.cols() - 1);
            double max_obj = 1;
            if (obj > max_obj)
            {
                max_obj = obj;
            }
            double t = (max_obj - obj) / dict.red_c(e_pos);
            dict.x(dict.N[e_pos]) = dict.x(dict.N[e_pos]) + t;
            dict.xN(e_pos) = dict.x(dict.N[e_pos]);
            for (int k = 0; k < dict.B.size(); k++)
            {
                dict.xB(k) -= t * dict.w(k);
                dict.x(dict.B[k]) = dict.xB(k);
            }
            double min = dict.x(dict.A.cols() - 1);
            if (min > dict.opt_requirement)
            {
                dict.x_shifted = dict.x;
                dict.has_shifted = true;
                dict.disablePricing();
                return true;
            }
            for (int i = 0; i < dict.n; i++)
            {
                if (dict.x(dict.A.cols() - 1) + dict.x(dict.d + i) < min)
                {
                    min = dict.x(dict.A.cols() - 1) + dict.x(dict.d + i);
                }
            }
            dict.disablePricing();
            return false;
        }
    }
    dict.disablePricing();
    return false;
}
void print_null2(const char *s) {}
#define Malloc2(type, n) (type *)malloc((n) * sizeof(type))
Eigen::VectorXd max_margin(double **X, vector<int> &orientation, int total_points, int d, Eigen::VectorXd &a_start, bool scale)
{
    vector<double> denom(d, 1);
    if (scale)
    {
        for (int j = 0; j < d; j++)
        {
            denom[j] = 0;
            double mean = 0;
            for (int i = 0; i < total_points; i++)
            {
                mean += X[i][j];
            }
            mean /= total_points;
            double dev = 0;
            for (int i = 0; i < total_points; i++)
            {
                dev += (X[i][j] - mean) * (X[i][j] - mean);
            }
            dev = sqrt(dev / (total_points - 1));
            if (dev != 0)
            {
                denom[j] = dev;
            }
            else
            {
                denom[j] = 1;
            }
        }
    }
    struct parameter param;
    struct problem prob;
    struct model *model;
    struct feature_node *x_space;
    int elements = d * total_points;
    prob.l = total_points;
    prob.bias = 1;
    prob.y = Malloc2(double, prob.l);
    prob.x = Malloc2(struct feature_node *, prob.l);
    x_space = Malloc2(struct feature_node, elements + prob.l);
    int j2 = 0;
    int max_index = 0;
    for (int i = 0; i < total_points; i++)
    {
        prob.x[i] = &x_space[j2];
        prob.y[i] = orientation[i];
        for (int j = 0; j < d - 1; j++)
        {
            x_space[j2].index = j + 1;
            x_space[j2].value = X[i][j] / denom[j];
            if (x_space[j2].index > max_index)
            {
                max_index = x_space[j2].index;
            }
            j2++;
        }
        if (prob.bias >= 0)
        {
            x_space[j2++].value = prob.bias;
        }
        x_space[j2++].index = -1;
    }
    if (prob.bias >= 0)
    {
        prob.n = max_index + 1;
        for (int i = 1; i < prob.l; i++)
        {
            (prob.x[i] - 2)->index = prob.n;
        }
        x_space[j2 - 2].index = prob.n;
    }
    else
    {
        prob.n = max_index;
    }
    double mindiff = std::numeric_limits<double>::max();
    for (int i = 0; i < total_points; i++)
    {
        double sum = 0;
        for (int j = 0; j < d; j++)
        {
            sum += a_start(j) * X[i][j];
        }
        if (abs(sum) < mindiff)
        {
            mindiff = abs(sum);
        }
    }
    if (mindiff > 0)
    {
        mindiff = 0.99 * mindiff;
    }
    else
    {
        mindiff = 1;
    }
    param.solver_type = L2R_L2LOSS_SVC;
    param.C = 1e10;
    param.p = 0.1;
    param.nu = 0.;
    param.eps = 1e-9;
    param.nr_weight = 0;
    param.regularize_bias = 0;
    param.weight_label = NULL;
    param.weight = NULL;
    param.init_sol = NULL;
    param.init_sol = Malloc2(double, d);
    for (int j = 0; j < d; j++)
    {
        param.init_sol[j] = a_start(j) / mindiff;
    }
    set_print_string_function(print_null2);
    model = train(&prob, &param);
    Eigen::VectorXd a_min(d);
    for (int j = 0; j < d; j++)
    {
        a_min(j) = model->w[j] / denom[j];
    }
    int errors = 0;
    for (int i = 0; i < total_points; i++)
    {
        double sum = 0;
        for (int j = 0; j < d; j++)
        {
            sum += a_min(j) * X[i][j];
        }
        if (sum * orientation[i] < 0)
        {
            errors++;
        }
    }
    free_and_destroy_model(&model);
    destroy_param(&param);
    free(prob.y);
    free(prob.x);
    free(x_space);
    return a_min;
}
Eigen::VectorXd max_margin_solver(int ***data, double ***data_reg, double **X2, double **X, vector<int> &orientation, int total_points, int d, Eigen::VectorXd &a_start, bool scale, bool l1, Environment *env)
{
    vector<double> denom(d, 1);
    vector<double> offset(d, 0);
    for (int j = 0; j < d - 1; j++)
    {
        double min = std::numeric_limits<double>::max();
        double max = std::numeric_limits<double>::min();
        double median;
        double uq;
        double lq;
        double dev = 0;
        if (!env->regression)
        {
            double *X1 = X2[j];
            int **dataX1 = data[j];
            if (total_points % 2 == 0)
            {
                median = 0.5 * (X1[dataX1[(total_points) / 2 - 1][0]] + X1[dataX1[(total_points) / 2][0]]);
            }
            else
            {
                median = X1[dataX1[(total_points - 1) / 2][0]];
            }
            if (total_points % 2 == 0)
            {
                if ((total_points / 2) % 2 == 0)
                {
                    int ind = (total_points) / 4;
                    lq = 0.5 * (X1[dataX1[ind - 1][0]] + X1[dataX1[ind][0]]);
                    uq = 0.5 * (X1[dataX1[total_points / 2 + ind - 1][0]] + X1[dataX1[total_points / 2 + ind][0]]);
                }
                else
                {
                    lq = X1[dataX1[(total_points / 2 - 1) / 2][0]];
                    uq = X1[dataX1[(total_points / 2 - 1) / 2 + total_points / 2][0]];
                }
            }
            else
            {
                if (((total_points + 1) / 2) % 2 == 0)
                {
                    int ind = ((total_points + 1) / 2) / 2;
                    lq = 0.5 * (X1[dataX1[ind - 1][0]] + X1[dataX1[ind][0]]);
                    uq = 0.5 * (X1[dataX1[ind - 1][0]] + X1[dataX1[ind][0]]);
                }
                else
                {
                    int ind = ((total_points + 1) / 2 - 1) / 2;
                    lq = X1[dataX1[ind][0]];
                    uq = X1[dataX1[ind][0]];
                }
            }
        }
        else
        {
            double *X1 = X2[j];
            double **dataX1 = data_reg[j];
            if (total_points % 2 == 0)
            {
                median = 0.5 * (X1[(int)dataX1[(total_points) / 2 - 1][0]] + X1[(int)dataX1[(total_points) / 2][0]]);
            }
            else
            {
                median = X1[(int)dataX1[(total_points - 1) / 2][0]];
            }
            if (total_points % 2 == 0)
            {
                if ((total_points / 2) % 2 == 0)
                {
                    int ind = (total_points) / 4;
                    lq = 0.5 * (X1[(int)dataX1[ind - 1][0]] + X1[(int)dataX1[ind][0]]);
                    uq = 0.5 * (X1[(int)dataX1[total_points / 2 + ind - 1][0]] + X1[(int)dataX1[total_points / 2 + ind][0]]);
                }
                else
                {
                    lq = X1[(int)dataX1[(total_points / 2 - 1) / 2][0]];
                    uq = X1[(int)dataX1[(total_points / 2 - 1) / 2 + total_points / 2][0]];
                }
            }
            else
            {
                if (((total_points + 1) / 2) % 2 == 0)
                {
                    int ind = ((total_points + 1) / 2) / 2;
                    lq = 0.5 * (X1[(int)dataX1[ind - 1][0]] + X1[(int)dataX1[ind][0]]);
                    uq = 0.5 * (X1[(int)dataX1[ind - 1][0]] + X1[(int)dataX1[ind][0]]);
                }
                else
                {
                    int ind = ((total_points + 1) / 2 - 1) / 2;
                    lq = X1[(int)dataX1[ind][0]];
                    uq = X1[(int)dataX1[ind][0]];
                }
            }
        }
        denom[j] = 0;
        double mean = 0;
        double maxim = 0;
        double minim = std::numeric_limits<double>::max();
        for (int i = 0; i < total_points; i++)
        {
            mean += X[i][j];
            if (abs(X[i][j]) > maxim)
            {
                maxim = abs(X[i][j]);
            }
            if (minim != 0 && abs(minim) < minim)
            {
                minim = abs(X[i][j]);
            }
        }
        mean /= total_points;
        for (int i = 0; i < total_points; i++)
        {
            dev += (X[i][j] - mean) * (X[i][j] - mean);
        }
        dev = sqrt(dev / (total_points - 1));
        double den = dev;
        if (den != 0)
        {
            denom[j] = den;
        }
        else
        {
            denom[j] = 1;
        }
    }
    #ifndef use_gurobi
        if(!l1){
            return max_margin(X, orientation, total_points, d, a_start, scale);
        }
        ClpSimplex clpmodel;
        clpmodel.resize(0, d+d-1);
        vector<int> a_ind(d);
        for(int j=0;j<d;j++){
            a_ind[j]=j;
        }
        vector<int> abs_ind(d-1);
        for(int j=0;j<d-1;j++){
            abs_ind[j]=d+j;
        }
        for(int j=0;j<d;j++){
            clpmodel.setColumnLower(a_ind[j],-COIN_DBL_MAX);
            clpmodel.setColumnUpper(a_ind[j],COIN_DBL_MAX);
            if (a_start(j) == 0 && j < d - 1)
            {
                clpmodel.setColumnLower(a_ind[j],0);
                clpmodel.setColumnUpper(a_ind[j],0);
            }
            if(j<d-1){
                clpmodel.setColumnLower(abs_ind[j],0);
                clpmodel.setColumnUpper(abs_ind[j],COIN_DBL_MAX);
            }
        }
        for(int j=0;j<d-1;j++){
            if(!scale){
                clpmodel.setObjectiveCoefficient(abs_ind[j],1/(denom[j] * denom[j]));
            }
            else{
                clpmodel.setObjectiveCoefficient(abs_ind[j],1);
            }
        }
        CoinBuild buildObject;
        int** rowInds=new int*[total_points];
        double** rowVals=new double*[total_points];
        for(int i=0;i<total_points;i++){
            double sum = 0;
            for (int j = 0; j < d; j++)
            {
                sum += (X[i][j]) * a_start(j);
            }
            rowInds[i]=new int[d];
            rowVals[i]=new double[d];
            for(int j=0;j<d;j++){
                rowInds[i][j]=a_ind[j];
                rowVals[i][j]=((X[i][j] - offset[j]) / denom[j]);
            }
            if(sum>=0){
                buildObject.addRow(d, rowInds[i], rowVals[i],1.0, COIN_DBL_MAX);
                //clpmodel.addRow(d, rowInds[i], rowVals[i],1.0, COIN_DBL_MAX);
            }
            else{
                buildObject.addRow(d, rowInds[i], rowVals[i],-COIN_DBL_MAX, -1.0);
                //clpmodel.addRow(d, rowInds[i], rowVals[i],-COIN_DBL_MAX, -1.0);
            }
        }
        int** rowInds2=new int*[d-1];
        double** rowVals2=new double*[d-1];
        for(int j=0;j<d-1;j++){
            rowInds2[j]=new int[2];
            rowVals2[j]=new double[2];
            rowInds2[j][0]=abs_ind[j];
            rowInds2[j][1]=a_ind[j];
            rowVals2[j][0]=1;
            rowVals2[j][1]=-1;
            buildObject.addRow(2, rowInds2[j], rowVals2[j],0, COIN_DBL_MAX);
            //clpmodel.addRow(2, rowInds2[j], rowVals2[j],0, COIN_DBL_MAX);
        }
        int** rowInds3=new int*[d-1];
        double** rowVals3=new double*[d-1];
        for(int j=0;j<d-1;j++){
            rowInds3[j]=new int[2];
            rowVals3[j]=new double[2];
            rowInds3[j][0]=abs_ind[j];
            rowInds3[j][1]=a_ind[j];
            rowVals3[j][0]=1;
            rowVals3[j][1]=1;
            buildObject.addRow(2, rowInds3[j], rowVals3[j],0, COIN_DBL_MAX);
            //clpmodel.addRow(2, rowInds3[j], rowVals3[j],0, COIN_DBL_MAX);
        }
        clpmodel.setLogLevel(0);
        clpmodel.addRows(buildObject);
        bool primal=false;
        if(primal){
            clpmodel.primal();
            if(!clpmodel.isProvenOptimal()){
                clpmodel.dual();
            }
        }
        else{
            clpmodel.dual();
            if(!clpmodel.isProvenOptimal()){
                clpmodel.primal();
            }
        }
        bool has_solution=clpmodel.isProvenOptimal();
        Eigen::VectorXd a_min = a_start;
        if(has_solution){
            double * columnPrimal = clpmodel.primalColumnSolution();
            a_min(d - 1) = columnPrimal[a_ind[d-1]] / denom[d - 1];
            for (int j = 0; j < d - 1; j++)
            {
                a_min(j) = columnPrimal[a_ind[j]]/ denom[j];
                a_min(d - 1) -= a_min(j) * (offset[j]);
            }
        }
        for(int i=0;i<total_points;i++){
            delete[] rowInds[i];
            delete[] rowVals[i];
        }
        delete[] rowInds;
        delete[] rowVals;
        for(int j=0;j<d-1;j++){
            delete[] rowInds2[j];
            delete[] rowVals2[j];
        }
        delete[] rowInds2;
        delete[] rowVals2;
        for(int j=0;j<d-1;j++){
            delete[] rowInds3[j];
            delete[] rowVals3[j];
        }
        delete[] rowInds3;
        delete[] rowVals3;
        if(has_solution){
            return a_min;
        }
        else{
            return max_margin(X, orientation, total_points, d, a_start, scale);
        }
    #endif
    #ifdef use_gurobi
        grb_env.start();
        GRBModel model = GRBModel(grb_env);
        model.set(GRB_IntParam_NumericFocus, 3);
        model.set(GRB_IntParam_OutputFlag, 0);
        model.set(GRB_IntParam_LogToConsole, 0);
        vector<GRBVar> x(d);
        for (int j = 0; j < d; j++)
        {
            if (a_start(j) == 0 && j < d - 1)
            {
                x[j] = model.addVar(-0, 0, 0.0, GRB_CONTINUOUS, "x" + j);
            }
            else
            {
                x[j] = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "x" + j);
            }
        }
        if (l1)
        {
            vector<GRBVar> abs(d);
            for (int j = 0; j < d; j++)
            {
                abs[j] = model.addVar(0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "abs" + j);
                model.addConstr(abs[j] >= x[j]);
                model.addConstr(abs[j] >= -x[j]);
            }
            GRBLinExpr obj = GRBLinExpr();
            for (int j = 0; j < d - 1; j++)
            {
                if (!scale)
                {
                    obj += abs[j] / (denom[j] * denom[j]);
                }
                else
                {
                    obj += abs[j];
                }
            }
            model.setObjective(obj, GRB_MINIMIZE);
        }
        else
        {
            GRBQuadExpr obj = GRBQuadExpr();
            for (int j = 0; j < d - 1; j++)
            {
                if (!scale)
                {
                    obj += (x[j] * x[j]) / (denom[j] * denom[j]);
                }
                else
                {
                    obj += x[j] * x[j];
                }
            }
        }
        for (int i = 0; i < total_points; i++)
        {
            double sum = 0;
            for (int j = 0; j < d; j++)
            {
                sum += (X[i][j]) * a_start(j);
            }
            GRBLinExpr e = GRBLinExpr();
            for (int j = 0; j < d; j++)
            {
                e += x[j] * ((X[i][j] - offset[j]) / denom[j]);
            }
            if (sum >= 0)
            {
                model.addConstr(e >= 1);
            }
            else
            {
                model.addConstr(e <= -1);
            }
        }
        model.optimize();
        Eigen::VectorXd a_min = a_start;
        int status = model.get(GRB_IntAttr_Status);
        if (status == GRB_OPTIMAL || status == GRB_SUBOPTIMAL)
        {
            a_min(d - 1) = x[d - 1].get(GRB_DoubleAttr_X) / denom[d - 1];
            for (int j = 0; j < d - 1; j++)
            {
                a_min(j) = x[j].get(GRB_DoubleAttr_X) / denom[j];
                a_min(d - 1) -= a_min(j) * (offset[j]);
            }
            return a_min;
        }
        else
        {
            return max_margin(X, orientation, total_points, d, a_start, scale);
        }
    #endif
}
Eigen::VectorXd reduce_features(double **X, vector<int> orientation, int total_points, int d, Eigen::VectorXd a_start, bool scale, Environment *env)
{
    VectorXd a_min(a_start);
    vector<double> denom(d, 1);
    if (scale)
    {
        for (int j = 0; j < d - 1; j++)
        {
            denom[j] = 0;
            double mean = 0;
            double maxim = 0;
            double minim = std::numeric_limits<double>::max();
            double dev = 0;
            for (int i = 0; i < total_points; i++)
            {
                mean += X[i][j];
                if (abs(X[i][j]) > maxim)
                {
                    maxim = abs(X[i][j]);
                }
                if (minim != 0 && abs(minim) < minim)
                {
                    minim = abs(X[i][j]);
                }
            }
            mean /= total_points;
            for (int i = 0; i < total_points; i++)
            {
                dev += (X[i][j] - mean) * (X[i][j] - mean);
            }
            dev = sqrt(dev / (total_points - 1));
            double den = dev;
            if (den != 0)
            {
                denom[j] = den;
            }
            else
            {
                denom[j] = 1;
            }
            a_min(j) = a_min(j) * denom[j];
        }
    }
    int count = 0;
    for (int j = 0; j < d; j++)
    {
        if (abs(a_min(j)) > 0)
        {
            count++;
        }
    }
    vector<int> order;
    for (int j = 0; j < d; j++)
    {
        order.push_back(j);
    }
    bool reduced = true;
    for (int l = 0; l < 2 * d && count > 1 && reduced; l++)
    {
        reduced = false;
        std::shuffle(order.begin(), order.end(), env->generator);
        for (int k = 0; k < d; k++)
        {
            int j = order[k];
            if (abs(a_min(j)) == 0)
            {
                continue;
            }
            double min = std::numeric_limits<double>::lowest();
            double max = std::numeric_limits<double>::max();
            double sum;
            for (int i = 0; i < total_points; i++)
            {
                sum = 0;
                if (X[i][j] == 0)
                {
                    continue;
                }
                for (int j2 = 0; j2 < d; j2++)
                {
                    sum += a_min(j2) * (X[i][j2] / denom[j2]);
                }
                bool case1 = true;
                if (sum < 0)
                {
                    case1 = false;
                }
                sum -= a_min(j) * (X[i][j] / denom[j]);
                double tmp = -sum / (X[i][j] / denom[j]);
                if (case1)
                {
                    if ((X[i][j] / denom[j]) > 0)
                    {
                        if (tmp > min)
                        {
                            min = tmp;
                        }
                    }
                    else
                    {
                        if (tmp < max)
                        {
                            max = tmp;
                        }
                    }
                }
                else
                {
                    if ((X[i][j] / denom[j]) < 0)
                    {
                        if (tmp > min)
                        {
                            min = tmp;
                        }
                    }
                    else
                    {
                        if (tmp < max)
                        {
                            max = tmp;
                        }
                    }
                }
            }
            if (min < 0 && max > 0)
            {
                reduced = true;
                a_min(j) = 0;
            }
            else
            {
                if (min > std::numeric_limits<double>::lowest() && max < std::numeric_limits<double>::max() && (max - min) > 1e-6)
                {
                    a_min(j) = min + (max - min) / 2.;
                }
                else if (min > std::numeric_limits<double>::lowest() && max == std::numeric_limits<double>::max())
                {
                    double tmp = env->Rand() * (1 - 0);
                    while (tmp < 1e-6)
                    {
                        tmp = env->Rand() * (1 - 0);
                    }
                    a_min(j) = min + tmp;
                }
                else if (min == std::numeric_limits<double>::lowest() && max < std::numeric_limits<double>::max())
                {
                    double tmp = env->Rand() * (1 - 0);
                    while (tmp < 1e-6)
                    {
                        tmp = env->Rand() * (1 - 0);
                    }
                    a_min(j) = max - tmp;
                }
            }
        }
        count = 0;
        for (int j = 0; j < d; j++)
        {
            if (abs(a_min(j)) > 0)
            {
                count++;
            }
        }
    }
    if (scale)
    {
        for (int j = 0; j < d; j++)
        {
            a_min(j) = a_min(j) / denom[j];
        }
    }
    return a_min;
}
void BrightSide(ctpl::thread_pool *p, int ***data, double ***data_reg, double **X2, vector<int> order, double **X, int *Y, double *Y_reg, int total_points, int d, int n_classes, double *distr_S_G, Result *res, Environment *env, int norm, double start_time)
{
    int id = res->generateId();
    double *offset = new double[d];
    double *denominator = new double[d];
    bool normal = false;
    double **X3 = new double *[total_points];
    bool sparse = true;
    double eps = 1e-9;
    bool start_warm = res->split_found;
    for (int i = 0; i < total_points; i++)
    {
        X3[i] = new double[d];
        for (int j = 0; j < d; j++)
        {
            X3[i][j] = X[i][j];
        }
    }
    if (!env->regression)
    {
        if (normal)
        {
            if (norm == 0)
            {
                normalize(data, X2, X, total_points, d, offset, denominator, 0);
            }
            else if (norm == 1)
            {
                normalize(data, X2, X, total_points, d, offset, denominator, 1);
            }
        }
    }
    else
    {
        if (normal)
        {
            if (norm == 0)
            {
                normalize_reg(data_reg, X2, X, total_points, d, offset, denominator, 0);
            }
            else if (norm == 1)
            {
                normalize_reg(data_reg, X2, X, total_points, d, offset, denominator, 1);
            }
        }
    }
    vector<bool> is_duplicated(total_points, false);
    int n = 0;
    int n2 = 0;
    std::map<int, vector<int>> duplicates_tmp;
    for (int i = 0; i < total_points; i++)
    {
        if (is_duplicated[i])
        {
            continue;
        }
        duplicates_tmp[i] = vector<int>();
        duplicates_tmp[i].push_back(i);
        n++;
        for (int i2 = i + 1; i2 < total_points; i2++)
        {
            bool flag = true;
            for (int j = 0; j < d; j++)
            {
                if (X[i][j] != X[i2][j])
                {
                    flag = false;
                    break;
                }
            }
            if (flag)
            {
                duplicates_tmp[i].push_back(i2);
                is_duplicated[i2] = true;
            }
        }
    }
    int opt_try = -1;
    bool resume = false;
    for (int restart = 0; restart < 1; restart++)
    {
        VectorXd a_start(d);
        VectorXd a_min(d);
        if (restart > 0)
        {
            start_warm = res->split_found;
        }
        double min_tmp = std::numeric_limits<double>::max();
        if (normal && start_warm)
        {
            a_start(d - 1) = -res->bestB;
            for (int j = 0; j < d - 1; j++)
            {
                a_start(j) = res->a[j] * denominator[j];
                a_start(d - 1) += res->a[j] * offset[j];
            }
        }
        else if (start_warm)
        {
            a_start(d - 1) = -res->bestB;
            for (int j = 0; j < d - 1; j++)
            {
                a_start(j) = res->a[j];
            }
        }
        else
        {
            a_start(d - 1) = env->Rand(eps, 1);
            for (int j = 0; j < d - 1; j++)
            {
                a_start(j) = env->Rand(-1, 1);
            }
        }
        bool has_max_diff = false;
        double max_diff = std::numeric_limits<double>::min();
        if (abs(a_start(d - 1)) < 1e-8)
        {
            double min_t = std::numeric_limits<double>::lowest();
            double max_t = std::numeric_limits<double>::max();
            for (int i = 0; i < total_points; i++)
            {
                double sum = 0;
                for (int j = 0; j < d; j++)
                {
                    sum += a_start(j) * X[i][j];
                }
                if (sum >= 0)
                {
                    sum -= a_start(d - 1) * X[i][d - 1];
                    if (-sum > min_t)
                    {
                        min_t = -sum;
                    }
                }
                else
                {
                    sum -= a_start(d - 1) * X[i][d - 1];
                    if (-sum < max_t)
                    {
                        max_t = -sum;
                    }
                }
            }
            a_start(d - 1) = min_t + (max_t - min_t) / 2.;
            if (abs(a_start(d - 1)) < 1e-10)
            {
                a_start(d - 1) = min_t + (max_t - min_t) / 4.;
            }
        }
        if (a_start(d - 1) > 0)
        {
            a_start = -a_start;
        }
        double tmp = abs(a_start(d - 1));
        a_start = a_start / tmp;
        min_tmp = evaluate(a_start, X, Y, Y_reg, total_points, n_classes, env->regression, env->criterion);
        if (min_tmp < res->bestImp)
        {
            if (normal)
            {
                for (int j = 0; j < d - 1; j++)
                {
                    a_start(d - 1) -= (a_start(j) * offset[j]) / denominator[j];
                    a_start(j) = a_start(j) / denominator[j];
                }
            }
            res->split_found = true;
            res->id = id;
            res->bestImp = min_tmp;
            res->bestB = -a_start(d - 1);
            for (int j = 0; j < d - 1; j++)
            {
                res->a[j] = a_start(j);
            }
        }
        d = d - 1;
        MatrixXd A(n, d);
        int i_tmp = 0;
        std::map<int, vector<int>> duplicates;
        for (int i = 0; i < total_points; i++)
        {
            if (is_duplicated[i])
            {
                continue;
            }
            duplicates[i_tmp] = duplicates_tmp[i];
            for (int j = 0; j < d; j++)
            {
                A(i_tmp, j) = X[i][j];
            }
            i_tmp++;
        }
        vector<int> B(n);
        for (int i = 0; i < n; i++)
        {
            B[i] = d + i;
        }
        VectorXd b = VectorXd::Ones(n);
        VectorXd orientation(n);
        for (int i = 0; i < n; i++)
        {
            double sum = 0;
            for (int j = 0; j < d; j++)
            {
                sum += a_start(j) * A(i, j);
            }
            if (sum < b(i))
            {
                orientation(i) = 1;
            }
            else
            {
                orientation(i) = -1;
            }
        }
        VectorXd c = VectorXd::Zero(n + d + 1);
        Dictionary dict(A, b, c, orientation, B, Y, Y_reg, n_classes, duplicates, env);
        double imp_start;
        double imp_a;
        double imp_a_n;
        double imp_min = std::numeric_limits<double>::max();
        bool phaseI = false;
        int phaseI_tries = 0;
        VectorXd x_start = dict.x;
        vector<int> B_start(dict.B);
        Permutation lex_order_start = Permutation(dict.lex_order);
        while (!phaseI && phaseI_tries < 2)
        {
            if (phaseI_tries != 0)
            {
                a_start(d) = env->Rand(eps, 1);
                for (int j = 0; j < d; j++)
                {
                    a_start(j) = env->Rand(-1, 1);
                }
                for (int i = 0; i < n; i++)
                {
                    double sum = 0;
                    for (int j = 0; j < d; j++)
                    {
                        sum += a_start(j) * A(i, j);
                    }
                    if (sum < b(i))
                    {
                        orientation(i) = 1;
                    }
                    else
                    {
                        orientation(i) = -1;
                    }
                }
                dict.setOrientation(orientation);
                dict.setBasis(B_start, x_start, lex_order_start);
            }
            phaseI = PhaseI(dict, false);
            if (phaseI)
            {
                if (!env->regression)
                {
                    imp_a = imp(dict.distribution, n_classes, env->criterion);
                }
                else
                {
                    imp_a = dict.reg[0].error(env->criterion) + dict.reg[1].error(env->criterion);
                }
                imp_min = imp_a;
                resume = true;
            }
            else
            {
                resume = false;
                imp_min = std::numeric_limits<double>::max();
                imp_a = imp_min;
            }
            phaseI_tries++;
        }
        if (!dict.has_LU)
        {
            resume = false;
            imp_min = std::numeric_limits<double>::max();
        }
        int start_iterations = 100;
        if(env->spx_start_iterations>0){
            start_iterations=env->spx_start_iterations;
        }
        if (env->spx_start_iterations_expression.is_compiled)
        {
            int tmp = round(env->spx_start_iterations_expression.evaluate(d, total_points));
            if (tmp > 0)
            {
                start_iterations = tmp;
            }
        }
        double T;
        double T0;
        double beta_T = env->spx_beta_T;
        double min_T;
        vector<double> start_diffs;
        double cl = 0;
        double c97 = 0;
        int start_orientation_changes = 0;
        bool start = true;
        bool set_opt = false;
        ;
        bool accept = false;
        int no_improvements = 0;
        int max_no_improvements = max(d, 100);
        if(env->spx_no_improvement>0){
            max_no_improvements=env->spx_no_improvement;
        }
        if (env->spx_no_improvement_expression.is_compiled)
        {
            int tmp = round(env->spx_no_improvement_expression.evaluate(d, total_points));
            if (tmp > 0)
            {
                max_no_improvements = tmp;
            }
        }
        int not_shiftable_counter = 0;
        bool force_pivots = false;
        bool is_congested = false;
        if (has_max_diff)
        {
            T = max_diff / (-log(0.85));
            T0 = T;
            start = false;
        }
        if (resume)
        {
            for (int j = 0; j < dict.d; j++)
            {
                if (!dict.is_basic[j])
                {
                    for (int e_pos = 0; e_pos < dict.N.size(); e_pos++)
                    {
                        if (dict.N[e_pos] == j)
                        {
                            dict.pivot(e_pos, 0, 0, false);
                        }
                    }
                }
            }
        }
        vector<int> B_opt = vector<int>(dict.B);
        VectorXd orientation_opt = dict.orientation;
        bool has_shifted_opt = dict.has_shifted;
        VectorXd x_opt = dict.x;
        VectorXd x_opt_shifted = dict.x_shifted;
        Permutation lex_order_opt = Permutation(dict.lex_order);
        vector<int> B_recov = vector<int>(dict.B);
        VectorXd orientation_recov = dict.orientation;
        bool has_shifted_recov = dict.has_shifted;
        VectorXd x_recov = dict.x;
        VectorXd x_recov_shifted = dict.x_shifted;
        Permutation lex_order_recov = Permutation(dict.lex_order);
        dict.cB = Eigen::VectorXd::Zero(dict.B.size());
        dict.cN = Eigen::VectorXd::Zero(dict.N.size());
        int it = 0;
        int opt_it = 0;
        int best_it = 0;
        bool neighborhood = false;
        if (start && start_diffs.size() >= start_iterations)
        {
            start = false;
            no_improvements = 0;
            sort(start_diffs.begin(), start_diffs.end());
            for (int k = 0; k < start_diffs.size(); k++)
            {
                if (start_diffs[k] != 0)
                {
                    cl = start_diffs[k];
                }
            }
            double ind = 0.97 * start_diffs.size();
            if (round(ind) == ind)
            {
                c97 = 0.5 * (start_diffs[(int)ind] + start_diffs[(int)ind + 1]);
            }
            else
            {
                c97 = start_diffs[(int)floor(ind) + 1];
            }
            T = c97;
            min_T = c97 / (-log(0.05));
            T = T / (-log(0.85));
            T0 = T;
        }
        while (imp_min > 0 && resume)
        {
            if (no_improvements >= max_no_improvements)
            {
                neighborhood = true;
                dict.setOrientation(orientation_opt);
                dict.setBasis(B_opt, x_opt, lex_order_opt);
                if (!dict.has_LU)
                {
                    resume = false;
                    imp_min = std::numeric_limits<double>::max();
                    break;
                }
                imp_a = imp_min;
                set_opt = false;
            }
            it++;
            double residual = 0;
            residual = dict.getResidual();
            if (dict.eta.size() > dict.max_eta || residual > dict.res_tol)
            {
                dict.recompute_decomposition(true);
                if (!dict.has_LU)
                {
                    dict.setOrientation(orientation_recov);
                    dict.setBasis(B_recov, x_recov, lex_order_recov);
                    if (!dict.has_LU)
                    {
                        resume = false;
                        imp_min = std::numeric_limits<double>::max();
                        break;
                    }
                }
                else
                {
                    dict.updateStates();
                    if (dict.is_basis_feasible)
                    {
                        B_recov = vector<int>(dict.B);
                        orientation_recov = dict.orientation;
                        has_shifted_recov = dict.has_shifted;
                        x_recov = dict.x;
                        x_recov_shifted = dict.x_shifted;
                        lex_order_recov = Permutation(dict.lex_order);
                    }
                }
            }
            if (env->Rand() < 0.5)
            {
                int e_pos = dict.get_cobasic_var_rand(0, true);
                int sign = 0;
                if (dict.N[e_pos] >= d)
                {
                    sign = 1;
                }
                int status = dict.pivot(e_pos, 0, sign, false);
                if (status == -2)
                {
                    dict.recompute_decomposition(true);
                    if (!dict.has_LU)
                    {
                        dict.setOrientation(orientation_recov);
                        dict.setBasis(B_recov, x_recov, lex_order_recov);
                        if (!dict.has_LU)
                        {
                            resume = false;
                            imp_min = std::numeric_limits<double>::max();
                            break;
                        }
                    }
                    continue;
                }
                if (status == -3)
                {
                    dict.setOrientation(orientation_recov);
                    dict.setBasis(B_recov, x_recov, lex_order_recov);
                    if (!dict.has_LU)
                    {
                        resume = false;
                        imp_min = std::numeric_limits<double>::max();
                        break;
                    }
                    if (!env->regression)
                    {
                        imp_a = imp(dict.distribution, n_classes, env->criterion);
                    }
                    else
                    {
                        imp_a = dict.reg[0].error(env->criterion) + dict.reg[1].error(env->criterion);
                    }
                    if (!dict.is_basis_feasible)
                    {
                        exit(0);
                    }
                }
                if (status == 0)
                {
                }
            }
            else
            {
                if (!start)
                {
                    opt_it++;
                }
                VectorXd orient = Eigen::VectorXd(dict.orientation);
                int i_push;
                if (neighborhood)
                {
                    i_push = dict.get_cobasic_slack(2);
                }
                else
                {
                    if (d >= 10 || false)
                    {
                        i_push = dict.get_cobasic_slack(1);
                    }
                    else
                    {
                        i_push = dict.get_cobasic_slack(0);
                    }
                }
                if (i_push == -1)
                {
                    continue;
                }
                int i_change = dict.push_slack(i_push);
                if (i_change < 0 && neighborhood)
                {
                    break;
                }
                if (i_change < 0)
                {
                    dict.recompute_decomposition(true);
                    if (!dict.has_LU)
                    {
                        dict.setOrientation(orientation_recov);
                        dict.setBasis(B_recov, x_recov, lex_order_recov);
                        if (!dict.has_LU)
                        {
                            resume = false;
                            imp_min = std::numeric_limits<double>::max();
                            break;
                        }
                    }
                    i_change = dict.push_slack(i_push);
                    if (i_change < 0)
                    {
                        break;
                    }
                }
                dict.change_side(i_change);
                if (!env->regression)
                {
                    imp_a_n = imp(dict.distribution, n_classes, env->criterion);
                }
                else
                {
                    imp_a_n = dict.reg[0].error(env->criterion) + dict.reg[1].error(env->criterion);
                }
                dict.change_side(i_change);
                if (neighborhood && imp_a_n >= imp_a)
                {
                    break;
                }
                double difference = imp_a - imp_a_n;
                if (start && start_diffs.size() < start_iterations)
                {
                    no_improvements = 0;
                    accept = true;
                    start_diffs.push_back(abs(difference));
                }
                else if (start && start_diffs.size() >= start_iterations)
                {
                    start = false;
                    set_opt = true;
                    no_improvements = 0;
                    sort(start_diffs.begin(), start_diffs.end());
                    for (int k = 0; k < start_diffs.size(); k++)
                    {
                        if (start_diffs[k] != 0)
                        {
                            cl = start_diffs[k];
                        }
                    }
                    double ind = 0.97 * start_diffs.size();
                    if (round(ind) == ind)
                    {
                        c97 = 0.5 * (start_diffs[(int)ind] + start_diffs[(int)ind + 1]);
                    }
                    else
                    {
                        c97 = start_diffs[(int)floor(ind) + 1];
                    }
                    T = c97;
                    min_T = c97 / (-log(0.05));
                    T = T / (-log(0.85));
                    T0 = T;
                }
                else
                {
                    accept = metropolis(imp_a, imp_a_n, T, env->Rand(), false);
                }
                bool test = false;
                bool shiftable = true;
                double T_save = T;
                if (shiftable && !start)
                {
                    T = max(beta_T * T, min_T);
                    not_shiftable_counter = 0;
                }
                if (!shiftable || !accept)
                {
                    if (shiftable && !start)
                    {
                        no_improvements++;
                    }
                    if (set_opt)
                    {
                        dict.setOrientation(orientation_opt);
                        dict.setBasis(B_opt, x_opt, lex_order_opt);
                        if (!dict.has_LU)
                        {
                            resume = false;
                            imp_min = std::numeric_limits<double>::max();
                            break;
                        }
                        imp_a = imp_min;
                        set_opt = false;
                        T = T0;
                    }
                    continue;
                }
                else
                {
                    int status = dict.update_orientation(i_push, i_change);
                    if (status == -1)
                    {
                        T = T_save;
                        continue;
                    }
                    if (!env->regression)
                    {
                        imp_a_n = imp(dict.distribution, n_classes, env->criterion);
                    }
                    else
                    {
                        imp_a_n = dict.reg[0].error(env->criterion) + dict.reg[1].error(env->criterion);
                    }
                    if (!dict.is_basis_feasible)
                    {
                        dict.setOrientation(orientation_recov);
                        dict.setBasis(B_recov, x_recov, lex_order_recov);
                        if (!dict.has_LU)
                        {
                            resume = false;
                            imp_min = std::numeric_limits<double>::max();
                            break;
                        }
                        if (!env->regression)
                        {
                            imp_a = imp(dict.distribution, n_classes, env->criterion);
                        }
                        else
                        {
                            imp_a = dict.reg[0].error(env->criterion) + dict.reg[1].error(env->criterion);
                        }
                        imp_a_n = imp_a;
                    }
                }
                imp_a = imp_a_n;
                if (imp_a < imp_min)
                {
                    best_it = opt_it;
                    no_improvements = 0;
                    B_opt = vector<int>(dict.B);
                    lex_order_opt = Permutation(dict.lex_order);
                    orientation_opt = dict.orientation;
                    has_shifted_opt = dict.has_shifted;
                    x_opt = dict.x;
                    x_opt_shifted = dict.x_shifted;
                    imp_min = imp_a;
                }
                else
                {
                    if (!start)
                    {
                        no_improvements++;
                    }
                }
                if (set_opt)
                {
                    dict.setOrientation(orientation_opt);
                    dict.setBasis(B_opt, x_opt, lex_order_opt);
                    if (!dict.has_LU)
                    {
                        resume = false;
                        imp_min = std::numeric_limits<double>::max();
                        break;
                    }
                    imp_a = imp_min;
                    set_opt = false;
                    T = T0;
                }
            }
        }
        if (imp_min < res->bestImp && resume)
        {
            VectorXd a_best = a_start;
            double imp_best;
            if (!has_shifted_opt)
            {
                dict.shift();
                a_best = dict.getUnscaled();
                imp_best = evaluate(a_best, X, Y, Y_reg, total_points, n_classes, env->regression, env->criterion);
            }
            else
            {
                a_best = dict.getUnscaled();
                imp_best = evaluate(a_best, X, Y, Y_reg, total_points, n_classes, env->regression, env->criterion);
            }
            if (abs(imp_best - imp_min) > 1e-5)
            {
                dict.setOrientation(orientation_opt);
                dict.setBasis(B_opt, x_opt, lex_order_opt);
                if (dict.has_LU)
                {
                    bool status = isCell(dict);
                    if (status)
                    {
                        VectorXd a_tmp = dict.getUnscaled();
                        double imp_tmp = evaluate(a_tmp, X, Y, Y_reg, total_points, n_classes, env->regression, env->criterion);
                        if (imp_tmp < imp_best)
                        {
                            a_best = a_tmp;
                            imp_best = imp_tmp;
                        }
                    }
                    else
                    {
                        imp_min = std::numeric_limits<double>::max();
                    }
                }
            }
            a_min = a_best;
            imp_min = imp_best;
            d = d + 1;
            bool should_reduce = true;
            if (normal)
            {
                for (int j = 0; j < d - 1; j++)
                {
                    a_min(d - 1) -= (a_min(j) * offset[j]) / denominator[j];
                    a_min(j) = a_min(j) / denominator[j];
                }
            }
            if (should_reduce)
            {
                vector<int> orientation(total_points);
                int i2;
                for (int i = 0; i < dict.n; i++)
                {
                    vector<int> D = dict.get_Duplicates(i);
                    for (int k = 0; k < D.size(); k++)
                    {
                        i2 = D[k];
                        orientation[i2] = -dict.orientation(i);
                    }
                }
                VectorXd a_best_r = reduce_features(X3, orientation, total_points, d, a_min, false, env);
                double imp_best_r = evaluate(a_best_r, X3, Y, Y_reg, total_points, n_classes, env->regression, env->criterion);
                if (imp_best_r == imp_min)
                {
                    a_min = a_best_r;
                }
                a_best = max_margin_solver(data, data_reg, X2, X3, orientation, total_points, d, a_min, true, true, env);
            }
            a_min = a_best;
            imp_best = evaluate(a_min, X3, Y, Y_reg, total_points, n_classes, env->regression, env->criterion);
            imp_min = imp_best;
            if (imp_min < res->bestImp)
            {
                opt_try = restart;
                res->split_found = true;
                res->id = id;
                res->bestImp = imp_min;
                res->bestB = -a_min(d - 1);
                for (int j = 0; j < d - 1; j++)
                {
                    res->a[j] = a_min(j);
                }
                double sum = res->bestB * res->bestB;
                double min = abs(res->bestB);
                if (min == 0)
                {
                    min = 1;
                }
                for (int j = 0; j < d - 1; j++)
                {
                    sum += res->a[j] * res->a[j];
                    if (abs(res->a[j]) < min && abs(res->a[j]) != 0)
                    {
                        min = abs(res->a[j]);
                    }
                }
                sum = sqrt(sum);
                double denom = sum;
                res->bestB = res->bestB / denom;
                for (int j = 0; j < d - 1; j++)
                {
                    res->a[j] = res->a[j] / denom;
                }
                VectorXd a_t(d);
                for (int j = 0; j < d - 1; j++)
                {
                    a_t(j) = res->a[j];
                }
                a_t(d - 1) = -res->bestB;
            }
        }
    }
    delete[] offset;
    delete[] denominator;
    for (int i = 0; i < total_points; i++)
    {
        for (int j = 0; j < d; j++)
        {
            X[i][j] = X3[i][j];
        }
        delete[] X3[i];
    }
    delete[] X3;
}
struct Sample
{
    double imp;
    Eigen::VectorXd a;
};
bool compareSamples(const Sample &a, const Sample &b)
{
    return a.imp < b.imp;
}
void SimpleCrossEntropy(ctpl::thread_pool *p, int ***data, double ***data_reg, double **X2, vector<int> order, double **X, int *Y, double *Y_reg, int total_points, int d, int d_orig, int n_classes, double *distr_S_G, Result *res, Environment *env, int norm, double start_time)
{
    double **X3 = new double *[total_points];
    for (int i = 0; i < total_points; i++)
    {
        X3[i] = new double[d];
        for (int j = 0; j < d; j++)
        {
            X3[i][j] = X[i][j];
        }
    }
    double lowerbound = 0;
    if (!env->regression && env->criterion == Environment::gini_imp)
    {
        vector<double> distribution(2 * n_classes, 0);
        for (int i = 0; i < total_points; i++)
        {
            distribution[Y[i]] += 1;
        }
        int k_max = 0;
        for (int k = 0; k < n_classes; k++)
        {
            if (distribution[k] > distribution[k_max])
            {
                k_max = k;
            }
        }
        distribution[k_max] = 0;
        lowerbound = imp(distribution, n_classes, env->criterion);
    }
    int id = res->generateId();
    double *offset = new double[d];
    double *denominator = new double[d];
    bool normal = env->normalize;
    if (!env->regression)
    {
        if (normal)
        {
            if (norm == 0)
            {
                normalize(data, X2, X, total_points, d, offset, denominator, 0);
            }
            else if (norm == 1)
            {
                normalize(data, X2, X, total_points, d, offset, denominator, 1);
            }
        }
    }
    else
    {
        if (normal)
        {
            if (norm == 0)
            {
                normalize_reg(data_reg, X2, X, total_points, d, offset, denominator, 0);
            }
            else if (norm == 1)
            {
                normalize_reg(data_reg, X2, X, total_points, d, offset, denominator, 1);
            }
        }
    }
    vector<vMF> vMFsamplers;
    vMF vMFsampler(env, d);
    vMFsamplers.push_back(vMFsampler);
    if (env->n_threads > 1)
    {
        for (int i = 1; i < env->n_threads; i++)
        {
            vMF vMFsampler(env, d);
            vMFsamplers.push_back(vMFsampler);
        }
    }
    Eigen::VectorXd mu(d);
    for (int i = 0; i < d; i++)
    {
        mu(i) = 0;
    }
    double UB;
    if (!env->regression)
    {
        vector<double> distro(2 * n_classes, 0);
        for (int i = 0; i < total_points; i++)
        {
            distro[Y[i]] += 1;
        }
        UB = imp(distro, n_classes, env->criterion);
    }
    else
    {
        VectorXd a_tmp(d);
        for (int j = 0; j < d; j++)
        {
            a_tmp(j) = 0;
        }
        UB = evaluate(a_tmp, X, Y, Y_reg, total_points, n_classes, env->regression, env->criterion);
    }
    vector<double> a(d, 0);
    VectorXd a_best(d);
    for (int j = 0; j < d; j++)
    {
        a_best(j) = 0.5;
    }
    double min_imp = std::numeric_limits<double>::max();
    double alpha = 0.2;
    if (env->cross_entropy_alpha != -1)
    {
        alpha = env->cross_entropy_alpha;
    }
    double frac = 0.1;
    if (env->cross_entropy_rho != -1)
    {
        frac = env->cross_entropy_rho;
    }
    int N = round(2 * (d_orig - 1) * log2(total_points));
    if (N < 100)
    {
        N = 100;
    }
    if (env->cross_entropy_samples != -1)
    {
        N = env->cross_entropy_samples;
    }
    if (env->cross_entropy_samples_expression.is_compiled)
    {
        int N_tmp = round(env->cross_entropy_samples_expression.evaluate(d_orig - 1, total_points));
        if (N_tmp > 0)
        {
            N = N_tmp;
        }
    }
    int evaluations = 0;
    int no_improvement_counter = 0;
    int no_gamma_improvement_counter = 0;
    int max_no_gamma_improvement = 3;
    if (env->cross_entropy_no_improvement != -1)
    {
        max_no_gamma_improvement = env->cross_entropy_no_improvement;
    }
    if (env->cross_entropy_no_improvement_expression.is_compiled)
    {
        int max_no_gamma_improvement_tmp = round(env->cross_entropy_no_improvement_expression.evaluate(d_orig - 1, total_points));
        if (max_no_gamma_improvement_tmp > 0)
        {
            max_no_gamma_improvement = max_no_gamma_improvement_tmp;
        }
    }
    double iterations = 0;
    vector<Sample> samples;
    double gamma = std::numeric_limits<double>::max();
    double gamma_min = std::numeric_limits<double>::max();
    double kappa = 0;
    VectorXd theta(mu.size());
    for (int j = 0; j < d; j++)
    {
        theta(j) = 0;
    }
    while (true)
    {
        if (no_gamma_improvement_counter > max_no_gamma_improvement)
        {
            break;
        }
        theta = kappa * mu;
        Eigen::MatrixXd sampleM;
        samples.clear();
        if (env->n_threads > 1)
        {
            sampleM = vMFsampler.rvMF(N, theta);
            vector<future<double>> fs;
            for (int g = 0; g < N; g++)
            {
                Sample s;
                s.a = sampleM.row(g);
                fs.push_back(p->push(evaluate_p, s.a, X, Y, Y_reg, total_points, n_classes, env->regression, env->criterion));
                samples.push_back(s);
                evaluations += 1;
            }
            for (int k = 0; k < fs.size(); k++)
            {
                fs[k].wait();
            }
            for (int g = 0; g < N; g++)
            {
                samples[g].imp = fs[g].get();
            }
        }
        else
        {
            sampleM = vMFsampler.rvMF(N, theta);
            for (int g = 0; g < N; g++)
            {
                Sample s;
                s.a = sampleM.row(g);
                s.imp = evaluate(s.a, X, Y, Y_reg, total_points, n_classes, env->regression, env->criterion);
                evaluations += 1;
                samples.push_back(s);
            }
        }
        std::sort(samples.begin(), samples.end(), compareSamples);
        if (samples[0].imp < min_imp)
        {
            min_imp = samples[0].imp;
            for (int j = 0; j < d; j++)
            {
                a_best(j) = samples[0].a(j);
            }
            no_improvement_counter = 0;
        }
        else
        {
            no_improvement_counter++;
        }
        gamma = samples[ceil((frac)*samples.size()) - 1].imp;
        if (gamma < gamma_min)
        {
            gamma_min = gamma;
            no_gamma_improvement_counter = 0;
        }
        else
        {
            no_gamma_improvement_counter++;
        }
        bool weighted_update = false;
        if (!weighted_update)
        {
            double total_weight = 0;
            VectorXd numerator(d);
            VectorXd mu_w(d);
            double kappa_w = 0;
            for (int j = 0; j < d; j++)
            {
                numerator(j) = 0;
            }
            double R;
            for (int g = 0; g < samples.size(); g++)
            {
                if (samples[g].imp > gamma)
                {
                    break;
                }
                numerator += samples[g].a;
                total_weight += 1;
            }
            R = numerator.lpNorm<2>() / total_weight;
            kappa_w = (R * ((double)d - R * R)) / (1 - R * R);
            mu_w = numerator / numerator.lpNorm<2>();
            mu = alpha * mu + (1 - alpha) * mu_w;
            mu.normalize();
            kappa = alpha * kappa + (1 - alpha) * kappa_w;
        }
        else
        {
            double m = 0;
            double b = 1;
            if (samples[0].imp < samples[samples.size() - 1].imp)
            {
                m = 1 / (samples[0].imp - samples[samples.size() - 1].imp);
            }
            b = 1 - m * samples[0].imp;
            double total_weight = 0;
            VectorXd numerator(d);
            for (int j = 0; j < d; j++)
            {
                numerator(j) = 0;
            }
            for (int g = 0; g < samples.size(); g++)
            {
                double w = m * samples[g].imp + b;
                numerator += w * samples[g].a;
                total_weight += w;
            }
            double R = numerator.lpNorm<2>() / total_weight;
            double kappa_w = (R * ((double)d - R * R)) / (1 - R * R);
            VectorXd mu_w(d);
            mu_w = numerator / numerator.lpNorm<2>();
            mu = alpha * mu + (1 - alpha) * mu_w;
            mu.normalize();
            kappa = alpha * kappa + (1 - alpha) * kappa_w;
        }
        iterations++;
    }
    if (normal)
    {
        for (int j = 0; j < d - 1; j++)
        {
            a_best(d - 1) -= (a_best(j) * offset[j]) / denominator[j];
            a_best(j) = a_best(j) / denominator[j];
        }
    }
    double imp_best = min_imp;
    VectorXd a_min = a_best;
    double imp_min = imp_best;
    bool should_reduce = true;
    if (should_reduce)
    {
        vector<int> orientation(total_points);
        int i2;
        for (int i = 0; i < total_points; i++)
        {
            double tmp = 0;
            for (int j = 0; j < d; j++)
            {
                tmp += a_min(j) * X3[i][j];
            }
            if (tmp < 0)
            {
                orientation[i] = 1;
            }
            else
            {
                orientation[i] = -1;
            }
        }
        VectorXd a_best_r = reduce_features(X3, orientation, total_points, d, a_min, false, env);
        double imp_best_r = evaluate(a_best_r, X3, Y, Y_reg, total_points, n_classes, env->regression, env->criterion);
        if (imp_best_r == imp_min)
        {
            a_min = a_best_r;
        }
        a_best = max_margin_solver(data, data_reg, X2, X3, orientation, total_points, d, a_min, true, true, env);
    }
    a_min = a_best;
    imp_best = evaluate(a_min, X3, Y, Y_reg, total_points, n_classes, env->regression, env->criterion);
    imp_min = imp_best;
    if (imp_min < res->bestImp)
    {
        res->split_found = true;
        res->id = id;
        res->bestImp = imp_min;
        res->bestB = -a_min(d - 1);
        for (int j = 0; j < d - 1; j++)
        {
            res->a[j] = a_min(j);
        }
        double sum = res->bestB * res->bestB;
        double min = abs(res->bestB);
        if (min == 0)
        {
            min = 1;
        }
        for (int j = 0; j < d - 1; j++)
        {
            sum += res->a[j] * res->a[j];
            if (abs(res->a[j]) < min && abs(res->a[j]) != 0)
            {
                min = abs(res->a[j]);
            }
        }
        sum = sqrt(sum);
        double denom = sum;
        res->bestB = res->bestB / denom;
        for (int j = 0; j < d - 1; j++)
        {
            res->a[j] = res->a[j] / denom;
        }
        VectorXd a_t(d);
        for (int j = 0; j < d - 1; j++)
        {
            a_t(j) = res->a[j];
        }
        a_t(d - 1) = -res->bestB;
    }
    delete[] offset;
    delete[] denominator;
    for (int i = 0; i < total_points; i++)
    {
        for (int j = 0; j < d; j++)
        {
            X[i][j] = X3[i][j];
        }
        delete[] X3[i];
    }
    delete[] X3;
}