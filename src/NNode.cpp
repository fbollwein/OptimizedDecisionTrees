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

#include "NNode.h"
#include <limits>
#include <algorithm>
#include "Environment.h"
#ifdef use_gurobi
#include "gurobi_c++.h"
#endif
NNode::NNode(int **A, vector<vector<int>> fixed, vector<int> &inner, vector<vector<double>> &distr, vector<double> &distr_i, int N, int L, int K, vector<double> hyp, int criterion, int depth, vector<double> minL, vector<double> maxL, bool is_binary, int sortclass, Environment *env)
{
    this->A = A;
    this->fixed = fixed;
    this->inner = inner;
    this->N = N;
    this->L = L;
    this->K = K;
    this->depth = depth;
    this->hyp = hyp;
    this->distr = distr;
    this->distr_i = distr_i;
    this->criterion = criterion;
    this->maxL = maxL;
    this->minL = minL;
    this->is_binary = is_binary;
    this->sortclass = sortclass;
    this->env = env;
    LB = 0;
    UB = std::numeric_limits<double>::max();
}
vector<NNode> NNode::branch()
{
    vector<NNode> nodes;
    if (inner.size() == 0)
    {
        return nodes;
    }
    bool is_sep = true;
    int i = inner[inner.size() - 1];
    bool choose_max_inc = false;
    if (choose_max_inc)
    {
        double max_inc = std::numeric_limits<double>::lowest();
        int max_i = i;
        int max_total = 0;
        double max_imp = 0;
        for (int i2 = 0; i2 < inner.size(); i2++)
        {
            double inc = 0;
            int total = 0;
            double impu = 0;
            vector<double> tmp(K, 0);
            for (int k = 0; k < K; k++)
            {
                total += A[inner[i2]][k];
                tmp[k] = A[inner[i2]][k];
            }
            impu = imp(tmp, K, criterion);
            double max_inc_l = std::numeric_limits<double>::lowest();
            ;
            double min_inc_l = std::numeric_limits<double>::max();
            ;
            int max_inc_l_ind = 0;
            int min_inc_l_ind = 0;
            for (int l = 0; l < L; l++)
            {
                for (int k = 0; k < K; k++)
                {
                    tmp[k] = distr[l][k] + A[inner[i2]][k];
                }
                double inc_tmp = imp(tmp, K, criterion) - impu;
                if (inc_tmp > max_inc_l)
                {
                    max_inc_l = inc_tmp;
                }
                if (inc_tmp < min_inc_l)
                {
                    min_inc_l = inc_tmp;
                }
                inc += inc_tmp;
            }
            if (impu != 0 && max_inc_l > max_inc)
            {
                i = inner[i2];
                max_inc = max_inc_l;
                max_total = total;
                max_imp = impu;
            }
        }
    }
    bool choose_max_class = false;
    if (choose_max_class)
    {
        double max_non = 0;
        int max_k = 0;
        vector<double> tmp(K, 0);
        for (int i2 = 0; i2 < inner.size(); i2++)
        {
            for (int k = 0; k < K; k++)
            {
                tmp[k] += A[inner[i2]][k];
            }
        }
        for (int k = 0; k < K; k++)
        {
            if (tmp[k] > max_non)
            {
                max_non = tmp[k];
                max_k = k;
            }
        }
        double maxim = 0;
        for (int i2 = 0; i2 < inner.size(); i2++)
        {
            if (A[inner[i2]][max_k] > maxim)
            {
                maxim = A[inner[i2]][max_k];
                i = inner[i2];
            }
        }
    }
    bool farest_away = false;
    if (env->Rand() > 0.5)
    {
        farest_away = true;
    }
    if (farest_away)
    {
        double maxim = 0;
        vector<double> tmp(K, 0);
        double sum = 0;
        for (int l = 0; l < L; l++)
        {
            for (int k = 0; k < K; k++)
            {
                tmp[k] += distr[l][k];
                sum += distr[l][k];
            }
        }
        for (int k = 0; k < K; k++)
        {
            tmp[k] /= sum;
        }
        for (int i2 = 0; i2 < inner.size(); i2++)
        {
            vector<double> tmp2(K, 0);
            sum = 0;
            for (int k = 0; k < K; k++)
            {
                tmp2[k] += A[inner[i2]][k];
                sum += A[inner[i2]][k];
            }
            if (sum != 0)
            {
                for (int k = 0; k < K; k++)
                {
                    tmp2[k] /= sum;
                }
            }
            double distance = 0;
            for (int k = 0; k < K; k++)
            {
                distance += (tmp[k] - tmp2[k]) * (tmp[k] - tmp2[k]);
            }
            if (sum * distance > maxim)
            {
                maxim = sum * distance;
                i = inner[i2];
            }
        }
    }
    int fixed_l = -1;
    double p_sort = 0;
    if (is_binary && L == 2)
    {
        double sum = 0;
        if (A[i][sortclass] > 0)
        {
            for (int k = 0; k < K; k++)
            {
                sum += A[i][k];
            }
            p_sort = A[i][sortclass] / sum;
            for (int l = 0; l < L; l++)
            {
                if (p_sort < maxL[l] && p_sort > minL[l])
                {
                    fixed_l = l;
                }
            }
        }
    }
    for (int l = 0; l < L; l++)
    {
        if (fixed_l != -1)
        {
            l = fixed_l;
        }
        bool flag = false;
        if (fixed[l].size() == 0)
        {
            flag = true;
        }
        vector<int> inner2;
        for (int i2 = 0; i2 < inner.size(); i2++)
        {
            if (inner[i2] != i)
            {
                inner2.push_back(inner[i2]);
            }
        }
        vector<vector<int>> fixed2;
        for (int l2 = 0; l2 < L; l2++)
        {
            vector<int> f;
            fixed2.push_back(f);
            for (int i2 = 0; i2 < fixed[l2].size(); i2++)
            {
                fixed2[l2].push_back(fixed[l2][i2]);
            }
        }
        fixed2[l].push_back(i);
        vector<double> hyp2(K, 0);
        for (int k = 0; k < K; k++)
        {
            hyp2[k] = hyp[k];
        }
        vector<double> distr2_i(K, 0);
        for (int k = 0; k < K; k++)
        {
            distr2_i[k] = distr_i[k] - A[i][k];
        }
        vector<vector<double>> distr2;
        for (int l2 = 0; l2 < L; l2++)
        {
            vector<double> d(K, 0);
            for (int k = 0; k < K; k++)
            {
                d[k] = distr[l2][k];
            }
            distr2.push_back(d);
        }
        for (int k = 0; k < K; k++)
        {
            distr2[l][k] += A[i][k];
        }
        vector<double> maxL2 = maxL;
        vector<double> minL2 = minL;
        if (is_binary)
        {
            if (p_sort > maxL2[l])
            {
                maxL2[l] = p_sort;
            }
            if (p_sort < minL2[l])
            {
                minL2[l] = p_sort;
            }
        }
        NNode node(A, fixed2, inner2, distr2, distr2_i, N, L, K, hyp2, criterion, depth + 1, minL2, maxL2, is_binary, sortclass, env);
        nodes.push_back(node);
        if (fixed_l != -1)
        {
            break;
        }
        if (flag)
        {
            break;
        }
    }
    if (nodes.size() > 0)
    {
    }
    return nodes;
}
#ifdef use_gurobi
bool NNode::isLinearSeparable()
{
    vector<int> Ls(L, 0);
    for (int l = 0; l < L; l++)
    {
        Ls[l] = l;
    }
    shuffle(Ls.begin(), Ls.end(), env->generator);
    int l1 = Ls[0];
    int l2 = Ls[1];
    vector<double> max1(K, 0);
    vector<double> min1(K, 1);
    vector<double> max2(K, 0);
    vector<double> min2(K, 1);
    for (int i2 = 0; i2 < fixed[l1].size(); i2++)
    {
        int i = fixed[l1][i2];
        double sum = 0;
        for (int k = 0; k < K; k++)
        {
            sum += A[i][k];
        }
        for (int k = 0; k < K; k++)
        {
            double p = A[i][k] / sum;
            if (p > max1[k])
            {
                max1[k] = p;
            }
            if (p < min1[k])
            {
                min1[k] = p;
            }
        }
    }
    for (int i2 = 0; i2 < fixed[l2].size(); i2++)
    {
        int i = fixed[l2][i2];
        double sum = 0;
        for (int k = 0; k < K; k++)
        {
            sum += A[i][k];
        }
        for (int k = 0; k < K; k++)
        {
            double p = A[i][k] / sum;
            if (p > max2[k])
            {
                max2[k] = p;
            }
            if (p < min2[k])
            {
                min2[k] = p;
            }
        }
    }
    bool overlap = true;
    for (int k = 0; k < K; k++)
    {
        if (min2[k] < min1[k] && max2[k] < min1[k])
        {
            overlap = false;
            break;
        }
        if (min2[k] > max1[k] && max2[k] > max1[k])
        {
            overlap = false;
            break;
        }
        if (min1[k] < min2[k] && max1[k] < min2[k])
        {
            overlap = false;
            break;
        }
        if (min1[k] > max2[k] && max1[k] > max2[k])
        {
            overlap = false;
            break;
        }
    }
    if (!overlap)
    {
        return true;
    }
    grb_env.start();
    GRBModel model = GRBModel(grb_env);
    model.set(GRB_IntParam_OutputFlag, 0);
    model.set(GRB_DoubleParam_BestObjStop, 1e-6);
    vector<GRBVar> d(K);
    for (int k = 0; k < K; k++)
    {
        d[k] = model.addVar(-1, 1, 0.0, GRB_CONTINUOUS, "x" + k);
    }
    GRBVar eps = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, "eps");
    for (int i2 = 0; i2 < fixed[l1].size(); i2++)
    {
        int i = fixed[l1][i2];
        GRBLinExpr e = GRBLinExpr();
        double sum = 0;
        for (int k = 0; k < K; k++)
        {
            sum += A[i][k];
        }
        for (int k = 0; k < K; k++)
        {
            e += d[k] * (A[i][k] / sum);
        }
        e += eps;
        model.addConstr(e <= 0, "c" + i);
    }
    for (int i2 = 0; i2 < fixed[l2].size(); i2++)
    {
        int i = fixed[l2][i2];
        GRBLinExpr e = GRBLinExpr();
        double sum = 0;
        for (int k = 0; k < K; k++)
        {
            sum += A[i][k];
        }
        for (int k = 0; k < K; k++)
        {
            e += d[k] * (A[i][k] / sum);
        }
        e -= eps;
        model.addConstr(e >= 0, "c" + i);
    }
    GRBLinExpr e = GRBLinExpr();
    model.setObjective(eps + 0, GRB_MAXIMIZE);
    model.optimize();
    if ((model.get(GRB_IntAttr_Status) == 2 && model.getObjective().getValue() == 0))
    {
        if (model.get(GRB_IntAttr_Status) == 2)
        {
        }
        return false;
    }
    return true;
}
#endif
void NNode::lowerbound(vector<double> impurities)
{
    for (int l = 0; l < L; l++)
    {
        if (fixed[l].size() != 0)
        {
            LB += imp(distr[l], K, criterion);
        }
    }
    double LB2 = LB;
    int max_i = 0;
    double maxim = 0;
    for (int i2 = 0; i2 < inner.size(); i2++)
    {
        LB += impurities[inner[i2]];
    }
    vector<double> tmp(K, 0);
    for (int i2 = 0; i2 < inner.size(); i2++)
    {
        for (int k = 0; k < K; k++)
        {
            tmp[k] += A[inner[i2]][k];
        }
    }
    if (criterion == Environment::gini_imp)
    {
        for (int l = 0; l < L - 1; l++)
        {
            int max_k = 0;
            double max = 0;
            for (int k = 0; k < K; k++)
            {
                if (tmp[k] > max)
                {
                    max = tmp[k];
                    max_k = k;
                }
            }
            tmp[max_k] = 0;
        }
        LB2 += imp(tmp, K, criterion);
        LB = max(LB, LB2);
    }
    if (criterion == Environment::entropy_imp)
    {
        double sum = 0;
        for (int k = 0; k < K; k++)
        {
            sum += tmp[k];
        }
        double additional = imp(tmp, K, criterion) - sum * log2(L);
        LB2 += additional;
        LB = max(LB, LB2);
    }
}
pair<double, vector<vector<int>>> NNode::heuristic()
{
    vector<vector<int>> part(L);
    for (int l = 0; l < L; l++)
    {
        for (int k = 0; k < K; k++)
        {
            part[l] = fixed[l];
        }
    }
    vector<int> inner2 = inner;
    vector<vector<double>> distr_p(L);
    for (int l = 0; l < L; l++)
    {
        distr_p[l] = distr[l];
    }
    vector<double> impurities_p;
    for (int l = 0; l < L; l++)
    {
        impurities_p.push_back(imp(distr_p[l], K, criterion));
    }
    for (int i2 = inner2.size() - 1; i2 >= 0; i2--)
    {
        int i = inner[i2];
        double min = std::numeric_limits<double>::max();
        int min_l = 0;
        for (int l = 0; l < L; l++)
        {
            for (int k = 0; k < K; k++)
            {
                distr_p[l][k] += A[i][k];
            }
            double diff = imp(distr_p[l], K, criterion) - impurities_p[l];
            if (diff < min)
            {
                min = diff;
                min_l = l;
            }
            for (int k = 0; k < K; k++)
            {
                distr_p[l][k] -= A[i][k];
            }
        }
        for (int k = 0; k < K; k++)
        {
            distr_p[min_l][k] += A[i][k];
        }
        part[min_l].push_back(i);
        impurities_p[min_l] = imp(distr_p[min_l], K, criterion);
    }
    vector<pair<int, int>> order;
    for (int l = 0; l < L; l++)
    {
        for (int i2 = 0; i2 < part[l].size(); i2++)
        {
            order.push_back(pair<int, int>(l, part[l][i2]));
        }
    }
    bool cont = true;
    double minimp = 0;
    for (int l = 0; l < L; l++)
    {
        minimp += impurities_p[l];
    }
    int counter = -1;
    while (cont)
    {
        counter++;
        if (counter > 2 * order.size())
        {
            break;
        }
        cont = false;
        std::shuffle(order.begin(), order.end(), env->generator);
        for (int i2 = 0; i2 < order.size(); i2++)
        {
            int i = order[i2].second;
            int l = order[i2].first;
            int l2;
            double min_diff = std::numeric_limits<double>::max();
            int min_l = 0;
            for (l2 = 0; l2 < L; l2++)
            {
                if (l2 == order[i2].first)
                {
                    continue;
                }
                for (int k = 0; k < K; k++)
                {
                    distr_p[l][k] -= A[i][k];
                    distr_p[l2][k] += A[i][k];
                }
                double diff = imp(distr_p[l], K, criterion) + imp(distr_p[l2], K, criterion) - impurities_p[l] - impurities_p[l2];
                if (diff < min_diff)
                {
                    min_diff = diff;
                    min_l = l2;
                }
                for (int k = 0; k < K; k++)
                {
                    distr_p[l][k] += A[i][k];
                    distr_p[l2][k] -= A[i][k];
                }
            }
            if (min_diff < -1e-6)
            {
                cont = true;
                for (int k = 0; k < K; k++)
                {
                    distr_p[min_l][k] += A[i][k];
                    distr_p[l][k] -= A[i][k];
                }
                int index = 0;
                for (int ind = 0; ind < part[l].size(); ind++)
                {
                    if (part[l][ind] == i)
                    {
                        index = ind;
                        break;
                    }
                }
                part[min_l].push_back(i);
                part[l].erase(part[l].begin() + index);
                order[i2].first = min_l;
                impurities_p[l] = imp(distr_p[l], K, criterion);
                impurities_p[min_l] = imp(distr_p[min_l], K, criterion);
                minimp = 0;
                for (int l3 = 0; l3 < L; l3++)
                {
                    minimp += impurities_p[l3];
                }
            }
        }
    }
    pair<double, vector<vector<int>>> partpair = pair<double, vector<vector<int>>>(minimp, part);
    return partpair;
}
void NNode::upperbound(vector<double> impurities)
{
    UB = std::numeric_limits<double>::max();
    if (inner.size() == 0)
    {
        UB = 0;
        for (int l = 0; l < L; l++)
        {
            if (fixed[l].size() != 0)
            {
                UB += imp(distr[l], K, criterion);
            }
        }
    }
}
int NNode::partition_size()
{
    int size = 0;
    for (int l = 0; l < L; l++)
    {
        if (fixed[l].size() > 0)
        {
            size += 1;
        }
    }
    return size;
}