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

#include "NBranchAndBound.h"
#include <mutex>
#include <chrono>
using namespace std;
bool compsecond_l(const pair<int, double> &a, const pair<int, double> &b)
{
    return a.second < b.second;
}
bool compsecond_g(const pair<int, double> &a, const pair<int, double> &b)
{
    return a.second > b.second;
}
std::mutex mtx14;
int partition_size(vector<vector<int>> fixed)
{
    int size = 0;
    for (int l = 0; l < fixed.size(); l++)
    {
        if (fixed[l].size() > 0)
        {
            size += 1;
        }
    }
    return size;
}
void NBranchAndBound_reg(double **Xnominal, double *Y_reg, double ***data_reg, int total_points, int domain, int which_feature, vector<vector<bool>> is_nominal_val_in_dataset, Rule *rule, Environment *env)
{
    int N = domain;
    int L = env->nominal_partitions;
    vector<double> sum;
    vector<double> sum2;
    vector<Vector_reg> A;
    for (int i = 0; i <= N; i++)
    {
        A.push_back(Vector_reg());
        if (env->criterion == Environment::mae)
        {
            A[i].use_mae();
        }
    }
    vector<bool> has_val(N + 1, false);
    int n_values = 0;
    for (int i2 = 0; i2 < total_points; i2++)
    {
        int i = (int)Xnominal[(int)data_reg[0][i2][0]][which_feature];
        double k = Y_reg[(int)data_reg[0][i2][0]];
        A[i].push_back(k);
        if (!has_val[i])
        {
            n_values++;
            has_val[i] = true;
        }
    }
    if (n_values < 2)
    {
        return;
    }
    vector<pair<int, double>> impurity_order;
    vector<pair<int, double>> largest_order;
    for (int i = 0; i <= N; i++)
    {
        double impurity = A[i].error(env->criterion);
        impurity_order.push_back(pair<int, double>(i, impurity));
        largest_order.push_back(pair<int, double>(i, A[i].size()));
    }
    std::sort(largest_order.begin(), largest_order.end(), compsecond_l);
    vector<Vector_reg> fixed_y;
    vector<vector<int>> fixed;
    for (int l = 0; l < L; l++)
    {
        vector<int> f;
        fixed.push_back(f);
        fixed_y.push_back(Vector_reg());
        if (env->criterion == Environment::mae)
        {
            fixed_y[l].use_mae();
        }
    }
    vector<int> inner;
    for (int i2 = 0; i2 <= N; i2++)
    {
        int i = largest_order[i2].first;
        if (A[i].size() > 0)
        {
            inner.push_back(i);
        }
    }
    vector<double> maxL(L, std::numeric_limits<double>::lowest());
    vector<double> minL(L, std::numeric_limits<double>::max());
    NNode_reg root(A, fixed, inner, fixed_y, maxL, minL, N, L, 0, env->criterion);
    root.lowerBound();
    vector<NNode_reg> que;
    que.push_back(root);
    while (que.size() > 0)
    {
        NNode_reg node = que[que.size() - 1];
        que.pop_back();
        if (node.LB < rule->imp || (node.LB == rule->imp && partition_size(node.fixed) < partition_size(rule->partition)))
        {
            node.upperBound();
            if (node.UB < std::numeric_limits<double>::max() && partition_size(node.fixed) > 1)
            {
                if (node.UB < rule->imp || (node.UB == rule->imp && partition_size(node.fixed) < partition_size(rule->partition)))
                {
                    vector<int> used;
                    for (int l = 0; l < node.fixed_best.size(); l++)
                    {
                        if (node.fixed_best[l].size() > 0)
                        {
                            used.push_back(l);
                        }
                        for (int i = 0; i <= N; i++)
                        {
                            if (!has_val[i])
                            {
                                int l2 = rand() % used.size();
                                node.fixed_best[l2].push_back(i);
                            }
                        }
                    }
                    rule->imp = node.UB;
                    rule->partition = node.fixed_best;
                    rule->split_found = true;
                    rule->nominal_ind = which_feature;
                }
            }
            vector<NNode_reg> nodes = node.branch();
            vector<pair<int, double>> nodes_i(nodes.size());
            for (int l = 0; l < nodes.size(); l++)
            {
                nodes[l].lowerBound();
                nodes_i[l].first = l;
                nodes_i[l].second = nodes[l].LB;
            }
            std::sort(nodes_i.begin(), nodes_i.end(), compsecond_g);
            for (int l2 = 0; l2 < nodes.size(); l2++)
            {
                int l = nodes_i[l2].first;
                if (nodes[l].LB < rule->imp || (nodes[l].LB == rule->imp && partition_size(nodes[l].fixed) < partition_size(rule->partition)))
                {
                    que.push_back(nodes[l]);
                }
            }
        }
    }
}
vector<NNode> process_class(int id, NNode &node, int N, vector<double> impurities, vector<bool> has_val, int which_feature, Rule *rule, Environment *env)
{
    vector<NNode> nodes2;
    bool enter = false;
    mtx14.lock();
    if (node.LB < rule->imp || (node.LB == rule->imp && partition_size(node.fixed) < partition_size(rule->partition)))
    {
        enter = true;
    }
    mtx14.unlock();
    if (enter)
    {
        node.upperbound(impurities);
        double imp = std::numeric_limits<double>::max();
        vector<vector<int>> part;
        if (node.UB < std::numeric_limits<double>::max())
        {
            imp = node.UB;
            part = node.fixed;
        }
        if (env->Rand() < 0.01)
        {
            pair<double, vector<vector<int>>> partpair = node.heuristic();
            if (partpair.first < imp)
            {
                imp = partpair.first;
                part = partpair.second;
            }
        }
        if (imp < std::numeric_limits<double>::max() && partition_size(part) > 1)
        {
            enter = false;
            mtx14.lock();
            if (imp < rule->imp || (imp == rule->imp && partition_size(part) < partition_size(rule->partition)))
            {
                enter = true;
            }
            mtx14.unlock();
            if (enter)
            {
                vector<int> used;
                for (int l = 0; l < part.size(); l++)
                {
                    if (part[l].size() > 0)
                    {
                        used.push_back(l);
                    }
                    for (int i = 0; i <= N; i++)
                    {
                        if (!has_val[i])
                        {
                            int l2 = rand() % used.size();
                            part[l2].push_back(i);
                        }
                    }
                }
                mtx14.lock();
                if (imp < rule->imp || (imp == rule->imp && partition_size(part) < partition_size(rule->partition)))
                {
                    rule->split_found = true;
                    rule->imp = imp;
                    rule->partition = part;
                    rule->nominal_ind = which_feature;
                }
                mtx14.unlock();
            }
        }
        vector<NNode> nodes;
        nodes = node.branch();
        vector<pair<int, double>> nodes_i(nodes.size());
        bool test_sep = false;
        for (int l = 0; l < nodes.size(); l++)
        {
            nodes[l].lowerbound(impurities);
            if (nodes[l].LB < rule->imp)
            {
                test_sep = true;
            }
            nodes_i[l].first = l;
            nodes_i[l].second = nodes[l].LB;
        }
        bool add = true;
        #ifdef use_gurobi
        if (false && test_sep && node.depth >= 2 * node.K && node.depth % node.K)
        {
            add = node.isLinearSeparable();
        }
        #endif
        if (add)
        {
            for (int l2 = 0; l2 < nodes.size(); l2++)
            {
                int l = nodes_i[l2].first;
                mtx14.lock();
                if (nodes[l].LB >= rule->imp || (nodes[l].LB == rule->imp && partition_size(nodes[l].fixed) < partition_size(rule->partition)))
                {
                    add = false;
                }
                mtx14.unlock();
                if (add)
                {
                    nodes2.push_back(nodes[l]);
                }
            }
        }
    }
    return nodes2;
}
void NBranchAndBound_class(double **Xnominal, int *Y, int ***data, int total_points, int n_classes, int domain, int which_feature, vector<vector<bool>> is_nominal_val_in_dataset, Rule *rule, Environment *env, ctpl::thread_pool *p)
{
    int K = n_classes;
    int N = domain;
    int L = env->nominal_partitions;
    int **A = new int *[N + 1];
    for (int i = 0; i <= N; i++)
    {
        A[i] = new int[K];
        for (int k = 0; k < K; k++)
        {
            A[i][k] = 0;
        }
    }
    set<int> vals;
    vector<bool> has_val(N + 1, false);
    int n_values = 0;
    int sortclass = -1;
    for (int i2 = 0; i2 < total_points; i2++)
    {
        int i = (int)Xnominal[data[0][i2][0]][which_feature];
        int k = Y[data[0][i2][0]];
        vals.insert(k);
        if (sortclass == -1)
        {
            sortclass = k;
        }
        A[i][k]++;
        if (!has_val[i])
        {
            n_values++;
            has_val[i] = true;
        }
    }
    if (n_values < 2)
    {
        return;
    }
    bool is_binary = false;
    if (vals.size() <= 2)
    {
        is_binary = true;
    }
    vector<double> maxL(L, std::numeric_limits<double>::lowest());
    vector<double> minL(L, std::numeric_limits<double>::max());
    vector<pair<int, double>> impurity_order;
    vector<pair<int, double>> largest_order;
    int max_class = 0;
    int max_t = 0;
    int min_t = std::numeric_limits<int>::min();
    int min_class = 0;
    for (int k = 0; k < K; k++)
    {
        int sum = 0;
        for (int i = 0; i <= N; i++)
        {
            sum += A[i][k];
        }
        if (sum > max_t)
        {
            max_t = sum;
            max_class = k;
        }
        if (sum < min_t)
        {
            min_t = sum;
            min_class = k;
        }
    }
    vector<double> impurities;
    vector<double> beforedist(K, 0);
    for (int i = 0; i <= N; i++)
    {
        for (int k = 0; k < K; k++)
        {
            beforedist[k] += A[i][k];
        }
    }
    double before_imp = imp(beforedist, K, env->criterion);
    for (int i = 0; i <= N; i++)
    {
        vector<double> d(K, 0);
        double sum = 0;
        double pmax = 0;
        double pmin = std::numeric_limits<double>::max();
        for (int k = 0; k < K; k++)
        {
            d[k] = A[i][k];
            sum += A[i][k];
            if (d[k] > pmax)
            {
                pmax = d[k];
            }
            if (d[k] < pmin && d[k] != 0)
            {
                pmin = d[k];
            }
        }
        pmax = pmax / sum;
        pmin = pmin / sum;
        double impurity = imp(d, K, env->criterion);
        double purity;
        if (env->criterion == Environment::gini_imp)
        {
            purity = sum - impurity;
        }
        if (env->criterion == Environment::entropy_imp)
        {
            purity = sum * log2(K) - impurity;
        }
        impurities.push_back(impurity);
        impurity_order.push_back(pair<int, double>(i, impurity));
        double rank = 0;
        if (sum > 0)
        {
            vector<double> removedist(K, 0);
            for (int k = 0; k < K; k++)
            {
                beforedist[k] - A[i][k];
            }
            double rank = before_imp - imp(removedist, K, env->criterion);
        }
        if (env->criterion == Environment::gini_imp)
        {
            rank = (pmax - pmin) * sum;
        }
        else if (env->criterion == Environment::entropy_imp)
        {
            rank = (pmax - pmin) * sum;
            rank = sum;
        }
        else
        {
            rank = sum;
        }
        largest_order.push_back(pair<int, double>(i, rank));
    }
    std::sort(largest_order.begin(), largest_order.end(), compsecond_l);
    vector<vector<int>> fixed;
    for (int l = 0; l < L; l++)
    {
        vector<int> f;
        fixed.push_back(f);
    }
    vector<int> inner;
    for (int i2 = 0; i2 <= N; i2++)
    {
        int i = largest_order[i2].first;
        if (has_val[i])
        {
            inner.push_back(i);
        }
    }
    vector<double> hyp(K, 0);
    vector<vector<double>> distr;
    for (int l2 = 0; l2 < L; l2++)
    {
        vector<double> d(K, 0);
        distr.push_back(d);
    }
    vector<double> distr_i(K, 0);
    for (int i = 0; i < N; i++)
    {
        for (int k = 0; k < K; k++)
        {
            distr_i[k] += A[i][k];
        }
    }
    bool order_by_heuristic = false;
    if (env->criterion == Environment::entropy_imp)
    {
        order_by_heuristic = true;
    }
    if (order_by_heuristic)
    {
        NNode node(A, fixed, inner, distr, distr_i, N, L, K, hyp, env->criterion, 0, minL, maxL, is_binary, sortclass, env);
        pair<double, vector<vector<int>>> partpair = node.heuristic();
        vector<vector<int>> part = partpair.second;
        map<int, int> in_part;
        vector<vector<double>> distr_l;
        set<int> vals;
        vector<int> order;
        map<int, double> Ns;
        for (int i2 = 0; i2 < inner.size(); i2++)
        {
            Ns[inner[i2]] = 0;
        }
        for (int l = 0; l < part.size(); l++)
        {
            distr_l.push_back(vector<double>(K, 0));
            for (int i2 = 0; i2 < part[l].size(); i2++)
            {
                int i = part[l][i2];
                vals.insert(i);
                for (int k = 0; k < K; k++)
                {
                    distr_l[l][k] += A[i][k];
                    Ns[i] += A[i][k];
                }
                in_part[i] = l;
            }
        }
        while (vals.size() > 0)
        {
            int remove;
            double minim = std::numeric_limits<double>::max();
            double maxim = std::numeric_limits<double>::lowest();
            for (set<int>::iterator it = vals.begin(); it != vals.end(); it++)
            {
                int i = *it;
                int l2 = in_part[i];
                for (int k = 0; k < K; k++)
                {
                    distr_l[l2][k] -= A[i][k];
                }
                double LB = 0;
                for (int l = 0; l < distr_l.size(); l++)
                {
                    LB += imp(distr_l[l], K, env->criterion);
                }
                LB += impurities[i];
                for (int i2 = 0; i2 < order.size(); i2++)
                {
                    LB += impurities[order[i2]];
                }
                if (LB > maxim)
                {
                    maxim = LB;
                    remove = i;
                }
                if (LB == maxim && Ns[i] > Ns[remove])
                {
                    maxim = LB;
                    remove = i;
                }
                for (int k = 0; k < K; k++)
                {
                    distr_l[l2][k] += A[i][k];
                }
            }
            order.push_back(remove);
            vals.erase(remove);
        }
        inner.clear();
        for (int i = 0; i < order.size(); i++)
        {
            inner.push_back(order[i]);
        }
    }
    NNode root(A, fixed, inner, distr, distr_i, N, L, K, hyp, env->criterion, 0, minL, maxL, is_binary, sortclass, env);
    root.lowerbound(impurities);
    deque<NNode> que;
    que.push_back(root);
    int thread_count = 0;
    bool multi = true;
    int n_threads = p->size();
    int id = 0;
    int iter = 0;
    bool deep_dive = true;
    double min_imp = rule->imp;
    std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
    std::chrono::time_point<std::chrono::high_resolution_clock> end;
    if (n_threads >= 2)
    {
        int free_threads = n_threads;
        std::list<future<vector<NNode>>> fs;
        while (que.size() > 0 || fs.size() > 0)
        {
            end = std::chrono::high_resolution_clock::now();
            double diff = chrono::duration_cast<chrono::milliseconds>(end - start).count();
            if (diff > env->nom_timelimit * 1000)
            {
                break;
            }
            iter++;
            while (que.size() > 0 && free_threads > 0 && free_threads <= n_threads)
            {
                NNode node = que[que.size() - 1];
                if (deep_dive)
                {
                    if (iter % 5 < 1)
                    {
                        node = que[que.size() - 1];
                        que.pop_back();
                    }
                    else
                    {
                        node = que[0];
                        que.pop_front();
                    }
                }
                else
                {
                    node = que[0];
                    que.pop_front();
                }
                id++;
                free_threads--;
                fs.push_back(p->push(process_class, node, N, impurities, has_val, which_feature, rule, env));
            }
            if (fs.size() > 0)
            {
                std::list<future<vector<NNode>>>::iterator it = fs.begin();
                bool found = false;
                while (it != fs.end())
                {
                    if (it->wait_for(std::chrono::seconds(0)) == std::future_status::ready)
                    {
                        free_threads++;
                        vector<NNode> nodes = it->get();
                        for (int l = 0; l < nodes.size(); l++)
                        {
                            que.push_back(nodes[l]);
                        }
                        found = true;
                        break;
                    }
                    it++;
                }
                if (found)
                {
                    fs.erase(it);
                }
            }
        }
        for (std::list<future<vector<NNode>>>::iterator it = fs.begin(); it != fs.end(); it++)
        {
            it->wait();
        }
    }
    else
    {
        while (que.size() > 0)
        {
            end = std::chrono::high_resolution_clock::now();
            double diff = chrono::duration_cast<chrono::milliseconds>(end - start).count();
            if (diff > env->nom_timelimit * 1000)
            {
                break;
            }
            iter++;
            NNode node = que[que.size() - 1];
            if (deep_dive)
            {
                if (iter % 5 < 1)
                {
                    node = que[que.size() - 1];
                    que.pop_back();
                }
                else
                {
                    node = que[0];
                    que.pop_front();
                }
            }
            else
            {
                node = que[0];
                que.pop_front();
            }
            if (node.LB < rule->imp || (node.LB == rule->imp && partition_size(node.fixed) < partition_size(rule->partition)))
            {
                node.upperbound(impurities);
                double imp = std::numeric_limits<double>::max();
                vector<vector<int>> part;
                if (node.UB < std::numeric_limits<double>::max())
                {
                    imp = node.UB;
                    part = node.fixed;
                }
                if (env->Rand() < 0.01)
                {
                    pair<double, vector<vector<int>>> partpair = node.heuristic();
                    if (partpair.first < imp)
                    {
                        imp = partpair.first;
                        part = partpair.second;
                    }
                }
                if (imp < std::numeric_limits<double>::max() && partition_size(part) > 1)
                {
                    if (imp < rule->imp || (imp == rule->imp && partition_size(part) < partition_size(rule->partition)))
                    {
                        vector<int> used;
                        for (int l = 0; l < part.size(); l++)
                        {
                            if (part[l].size() > 0)
                            {
                                used.push_back(l);
                            }
                            for (int i = 0; i <= N; i++)
                            {
                                if (!has_val[i])
                                {
                                    int l2 = rand() % used.size();
                                    part[l2].push_back(i);
                                }
                            }
                        }
                        rule->split_found = true;
                        rule->imp = imp;
                        rule->partition = part;
                        rule->nominal_ind = which_feature;
                    }
                }
                vector<NNode> nodes = node.branch();
                vector<pair<int, double>> nodes_i(nodes.size());
                bool test_sep = false;
                for (int l = 0; l < nodes.size(); l++)
                {
                    nodes[l].lowerbound(impurities);
                    if (nodes[l].LB < rule->imp)
                    {
                        test_sep = true;
                    }
                    nodes_i[l].first = l;
                    nodes_i[l].second = nodes[l].LB;
                }
                bool add = true;
                #ifdef use_gurobi
                if (false && test_sep && node.depth >= 2 * K && node.depth % K)
                {
                    add = node.isLinearSeparable();
                }
                #endif
                if (add)
                {
                    std::sort(nodes_i.begin(), nodes_i.end(), compsecond_g);
                    for (int l2 = 0; l2 < nodes.size(); l2++)
                    {
                        int l = nodes_i[l2].first;
                        if (nodes[l].LB < rule->imp || (nodes[l].LB == rule->imp && partition_size(nodes[l].fixed) < partition_size(rule->partition)))
                        {
                            que.push_back(nodes[l]);
                        }
                    }
                }
            }
        }
    }
    end = std::chrono::high_resolution_clock::now();
    double diff = chrono::duration_cast<chrono::milliseconds>(end - start).count();
    rule->time = diff;
    if (que.size() > 0)
    {
        double minim = std::numeric_limits<double>::max();
        for (int i = 0; i < que.size(); i++)
        {
            if (que[i].LB < minim)
            {
                minim = que[i].LB;
            }
        }
        rule->opt_gap = 1 - minim / rule->imp;
    }
    for (int i = 0; i <= N; i++)
    {
        delete[] A[i];
    }
    delete[] A;
}
