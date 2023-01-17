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

#include "Node.h"
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sys/time.h>
#include <algorithm>
#include <limits>
#include <string>
#include "BranchAndBound.h"
#include "CrossBranchAndBound.h"
#include "NBranchAndBound.h"
#include "Heuristics.h"
#include <thread>
#include "linear.h"
#include "Eigen/Dense"
#ifndef use_gurobi
  #include "ClpSimplex.hpp"
  #include "CoinBuild.hpp"
#endif
using namespace std;
int compare_X(const void *arg1, const void *arg2)
{
  double const *lhs = *(const double **)(arg1);
  double const *rhs = *(const double **)(arg2);
  if (lhs[0] < rhs[0])
  {
    return -1;
  }
  if (rhs[0] < lhs[0])
  {
    return 1;
  }
  else
  {
    return 0;
  }
}
bool compareP(const pair<int, double> &a, const pair<int, double> &b)
{
  return a.second < b.second;
}
void Node::getData()
{
  shouldSplit = true;
  shouldSplit_numeric = true;
  shouldSplit_nominal = true;
  if (total_points < 2)
  {
    shouldSplit = false;
  }
  if (!regression)
  {
    Ni = new double[n_classes];
    for (int i = 0; i < n_classes; i++)
    {
      Ni[i] = 0;
    }
    if (total_points != 0)
    {
      maxVote = 0;
      int maxNi = 0;
      for (int i = 0; i < n_classes; i++)
      {
        Ni[i] = 0;
      }
      for (int i = 0; i < total_points; i++)
      {
        Ni[dataY[0][i][1]] += 1;
        if (Ni[dataY[0][i][1]] > maxNi)
        {
          maxVote = dataY[0][i][1];
          maxNi = Ni[dataY[0][i][1]];
        }
      }
      classification_error = total_points - Ni[maxVote];
      impurity = imp1(Ni, n_classes, env->criterion);
      int tmp = Ni[maxVote];
      if (total_points - tmp != 0 && env->criterion == Environment::gini_imp)
      {
        Ni[maxVote] = 0;
        LB = imp1(Ni, n_classes, env->criterion);
        Ni[maxVote] = tmp;
      }
      else
      {
        LB = 0;
      }
      int different_classes = 0;
      for (int i = 0; i < n_classes; i++)
      {
        if (Ni[i] > 0)
        {
          different_classes++;
        }
      }
      if (different_classes <= 1)
      {
        shouldSplit = false;
      }
    }
    else
    {
      impurity = 0;
      classification_error = 0;
      shouldSplit = false;
    }
  }
  else
  {
    LB = 0;
    mean = 0;
    median = 0;
    bool has_different = false;
    for (int i = 0; i < total_points; i++)
    {
      mean += dataY_reg[0][i][1];
      if (i > 0 && dataY_reg[0][i][1] != dataY_reg[0][i - 1][1])
      {
        has_different = true;
      }
    }
    if (total_points > 0)
    {
      mean = mean / (double)total_points;
      if (env->criterion == Environment::mae)
      {
        vector<double> vec(total_points);
        for (int i = 0; i < total_points; i++)
        {
          vec[i] = dataY_reg[0][i][1];
        }
        std::sort(vec.begin(), vec.end());
        int length = vec.size();
        if (length % 2 == 0)
        {
          median = (vec[length / 2 - 1] + vec[length / 2]) / 2;
        }
        else
        {
          median = vec[length / 2];
        }
      }
    }
    else
    {
      mean = 0;
      median = 0;
    }
    mse = 0;
    mae = 0;
    a_error = 0;
    sq_error = 0;
    for (int i = 0; i < total_points; i++)
    {
      double tmp;
      if (env->criterion == Environment::mse)
      {
        tmp = (mean - dataY_reg[0][i][1]);
      }
      else
      {
        tmp = (median - dataY_reg[0][i][1]);
      }
      mae += abs(tmp);
      a_error += abs(tmp);
      sq_error += tmp * tmp;
      mse += tmp * tmp;
    }
    if (total_points > 0)
    {
      mse = mse / (double)total_points;
      mae = mae / (double)total_points;
    }
    else
    {
      mse = 0;
      mae = 0;
    }
    if (!has_different)
    {
      shouldSplit = false;
    }
  }
  if (shouldSplit)
  {
    bool breakf = false;
    if (total_numeric_features == 0)
    {
      shouldSplit_numeric = false;
    }
    else
    {
      for (int j = 0; !breakf && j < total_numeric_features; j++)
      {
        for (int i = total_points - 1; !breakf && i < total_points; i++)
        {
          if (!regression)
          {
            if (X[j][dataY[j][i][0]] != X[j][dataY[j][0][0]])
            {
              breakf = true;
            }
          }
          else
          {
            if (X[j][(int)dataY_reg[j][i][0]] != X[j][(int)dataY_reg[j][0][0]])
            {
              breakf = true;
            }
          }
        }
      }
      if (!breakf)
      {
        shouldSplit_numeric = false;
      }
    }
    breakf = false;
    if (total_nominal_features == 0)
    {
      shouldSplit_nominal = false;
    }
    else
    {
      for (int j = 0; !breakf && j < total_nominal_features; j++)
      {
        for (int i = 0; !breakf && i < total_points; i++)
        {
          if (!regression)
          {
            if ((int) Xnominal[dataY[0][i][0]][j] != (int) Xnominal[dataY[0][0][0]][j])
            {
              breakf = true;
            }
          }
          else
          {
            if ((int) Xnominal[(int)dataY_reg[0][i][0]][j] != (int) Xnominal[(int)dataY_reg[0][0][0]][j])
            {
              breakf = true;
            }
          }
        }
      }
      if (!breakf)
      {
        shouldSplit_nominal = false;
      }
    }
    if (shouldSplit_nominal || shouldSplit_numeric)
    {
      shouldSplit = true;
    }
    else
    {
      shouldSplit = false;
    }
  }
}
void Node::getDataOC1()
{
  total_points = 0;
  maxVote = 0;
  int maxNi = 0;
  for (int m = 0; m < n_classes; m++)
  {
    total_points += Ni[m];
    if (Ni[m] > maxNi)
    {
      maxVote = m;
      maxNi = Ni[m];
    }
  }
  impurity = imp1(Ni, n_classes, env->criterion);
  int tmp = Ni[maxVote];
  if (total_points - tmp != 0 && env->criterion == Environment::gini_imp)
  {
    Ni[maxVote] = 0;
    LB = imp1(Ni, n_classes, env->criterion);
    Ni[maxVote] = tmp;
  }
  else
  {
    LB = 0;
  }
}
Node::Node() {}
Node::Node(double **X, int *Y, double *Y_reg, int total_points, int total_numeric_features, int n_classes, ctpl::thread_pool *p, Environment *env)
{
  vector<int> domain;
  vector<vector<bool>> is_nominal_val_in_dataset;
  initialize(X, NULL, Y, Y_reg, total_points, total_numeric_features, 0, domain, is_nominal_val_in_dataset, n_classes, p, env);
}
Node::Node(double **X, double **Xnominal, int *Y, double *Y_reg, int total_points, int total_numeric_features, int total_nominal_features, vector<int> nominal_domain, vector<vector<bool>> is_nominal_val_in_dataset, int n_classes, ctpl::thread_pool *p, Environment *env)
{
  initialize(X, Xnominal, Y, Y_reg, total_points, total_numeric_features, total_nominal_features, nominal_domain, is_nominal_val_in_dataset, n_classes, p, env);
}
void Node::initialize(double **X, double **Xnominal, int *Y, double *Y_reg, int total_points, int total_numeric_features, int total_nominal_features, vector<int> nominal_domain, vector<vector<bool>> is_nominal_val_in_dataset, int n_classes, ctpl::thread_pool *p, Environment *env)
{
  this->env = env;
  this->X = X;
  this->Y = Y;
  this->Y_reg = Y_reg;
  hasData = true;
  regression = env->regression;
  this->total_points = total_points;
  this->total_numeric_features = total_numeric_features;
  this->total_features = total_numeric_features + total_nominal_features;
  this->n_classes = n_classes;
  this->id = -1;
  this->p = p;
  this->total_nominal_features = total_nominal_features;
  this->Xnominal = Xnominal;
  this->nominal_domain = nominal_domain;
  this->is_nominal_val_in_dataset = is_nominal_val_in_dataset;
  shouldSplit = true;
  shouldSplit_numeric = true;
  shouldSplit_nominal = true;
  isLeaf = true;
  bestInd1 = 0;
  bestInd2 = 1;
  bestX1 = 0;
  bestX2 = 0;
  impS = 0;
  errorC = 0;
  if (total_numeric_features > 0)
  {
    double ***data = new double **[total_numeric_features];
    for (int j = 0; j < total_numeric_features; j++)
    {
      data[j] = new double *[total_points];
      for (int i = 0; i < total_points; i++)
      {
        data[j][i] = new double[3];
      }
    }
    for (int i = 0; i < total_points; i++)
    {
      for (int j = 0; j < total_numeric_features; j++)
      {
        data[j][i][0] = X[j][i];
        data[j][i][1] = i;
        if (!regression)
        {
          data[j][i][2] = Y[i];
        }
        else
        {
          data[j][i][2] = Y_reg[i];
        }
      }
    }
    for (int j = 0; j < total_numeric_features; j++)
    {
      qsort(data[j], total_points, sizeof data[j], compare_X);
    }
    if (!regression)
    {
      dataY = new int **[total_numeric_features];
      for (int j = 0; j < total_numeric_features; j++)
      {
        dataY[j] = new int *[total_points];
        for (int i = 0; i < total_points; i++)
        {
          dataY[j][i] = new int[2];
          dataY[j][i][0] = (int)data[j][i][1];
          dataY[j][i][1] = (int)data[j][i][2];
        }
      }
      for (int j = 0; j < total_numeric_features; j++)
      {
        for (int i = 0; i < total_points; i++)
        {
          delete[] data[j][i];
        }
        delete[] data[j];
      }
      delete[] data;
    }
    else
    {
      dataY_reg = new double **[total_numeric_features];
      for (int j = 0; j < total_numeric_features; j++)
      {
        dataY_reg[j] = new double *[total_points];
        for (int i = 0; i < total_points; i++)
        {
          dataY_reg[j][i] = new double[2];
          dataY_reg[j][i][0] = data[j][i][1];
          dataY_reg[j][i][1] = data[j][i][2];
        }
      }
      for (int j = 0; j < total_numeric_features; j++)
      {
        for (int i = 0; i < total_points; i++)
        {
          delete[] data[j][i];
        }
        delete[] data[j];
      }
      delete[] data;
    }
  }
  else
  {
    double ***data = new double **[1];
    for (int j = 0; j < 1; j++)
    {
      data[j] = new double *[total_points];
      for (int i = 0; i < total_points; i++)
      {
        data[j][i] = new double[3];
      }
    }
    for (int i = 0; i < total_points; i++)
    {
      for (int j = 0; j < 1; j++)
      {
        data[j][i][0] = i;
        data[j][i][1] = i;
        if (!regression)
        {
          data[j][i][2] = Y[i];
        }
        else
        {
          data[j][i][2] = Y_reg[i];
        }
      }
    }
    if (!regression)
    {
      dataY = new int **[1];
      for (int j = 0; j < 1; j++)
      {
        dataY[j] = new int *[total_points];
        for (int i = 0; i < total_points; i++)
        {
          dataY[j][i] = new int[2];
          dataY[j][i][0] = (int)data[j][i][1];
          dataY[j][i][1] = (int)data[j][i][2];
        }
      }
      for (int j = 0; j < total_numeric_features; j++)
      {
        for (int i = 0; i < total_points; i++)
        {
          delete[] data[j][i];
        }
        delete[] data[j];
      }
      delete[] data;
    }
    else
    {
      dataY_reg = new double **[1];
      for (int j = 0; j < 1; j++)
      {
        dataY_reg[j] = new double *[total_points];
        for (int i = 0; i < total_points; i++)
        {
          dataY_reg[j][i] = new double[2];
          dataY_reg[j][i][0] = data[j][i][1];
          dataY_reg[j][i][1] = data[j][i][2];
        }
      }
      for (int j = 0; j < total_numeric_features; j++)
      {
        for (int i = 0; i < total_points; i++)
        {
          delete[] data[j][i];
        }
        delete[] data[j];
      }
      delete[] data;
    }
  }
  isInTree = true;
  getData();
}
Node::Node(double **X, int *Y, double *Y_reg, int ***dataY, double ***dataY_reg, int total_points, int total_numeric_features, int n_classes, ctpl::thread_pool *p, Environment *env)
{
  vector<int> domain;
  vector<vector<bool>> is_nominal_val_in_dataset;
  initialize(X, NULL, Y, Y_reg, dataY, dataY_reg, total_points, total_numeric_features, 0, domain, is_nominal_val_in_dataset, n_classes, p, env);
}
Node::Node(double **X, double **Xnominal, int *Y, double *Y_reg, int ***dataY, double ***dataY_reg, int total_points, int total_numeric_features, int total_nominal_features, vector<int> nominal_domain, vector<vector<bool>> is_nominal_val_in_dataset, int n_classes, ctpl::thread_pool *p, Environment *env)
{
  initialize(X, Xnominal, Y, Y_reg, dataY, dataY_reg, total_points, total_numeric_features, total_nominal_features, nominal_domain, is_nominal_val_in_dataset, n_classes, p, env);
}
void Node::initialize(double **X, double **Xnominal, int *Y, double *Y_reg, int ***dataY, double ***dataY_reg, int total_points, int total_numeric_features, int total_nominal_features, vector<int> nominal_domain, vector<vector<bool>> is_nominal_val_in_dataset, int n_classes, ctpl::thread_pool *p, Environment *env)
{
  this->env = env;
  this->X = X;
  this->Y = Y;
  this->Y_reg = Y_reg;
  hasData = true;
  regression = env->regression;
  this->total_points = total_points;
  this->total_numeric_features = total_numeric_features;
  this->total_features = total_numeric_features + total_nominal_features;
  this->n_classes = n_classes;
  this->dataY = dataY;
  this->dataY_reg = dataY_reg;
  this->id = -1;
  this->p = p;
  this->total_nominal_features = total_nominal_features;
  this->Xnominal = Xnominal;
  this->nominal_domain = nominal_domain;
  this->is_nominal_val_in_dataset = is_nominal_val_in_dataset;
  a = vector<double>(0, 0);
  a.clear();
  shouldSplit = true;
  shouldSplit_numeric = true;
  shouldSplit_nominal = true;
  isLeaf = true;
  bestInd1 = 0;
  bestInd2 = 1;
  bestX1 = 0;
  bestX2 = 0;
  bestA = 0;
  bestB = 0;
  impS = 0;
  errorC = 0;
  isInTree = true;
  regression = env->regression;
  getData();
}
void Node::deleteDataOnly()
{
  if (hasData && possesses_data)
  {
    possesses_data = false;
    hasData = false;
    if (total_numeric_features > 0)
    {
      if (!regression)
      {
        for (int j = 0; j < total_numeric_features; j++)
        {
          for (int i = 0; i < total_points; i++)
          {
              delete[] dataY[j][i];
          }
          delete[] dataY[j];
        }
        delete[] dataY;
      }
      else
      {
        for (int j = 0; j < total_numeric_features; j++)
        {
          for (int i = 0; i < total_points; i++)
          {
            delete[] dataY_reg[j][i];
          }
          delete[] dataY_reg[j];
        }
        delete[] dataY_reg;
      }
    }
    else
    {
      if (!regression)
      {
        for (int j = 0; j < 1; j++)
        {
          for (int i = 0; i < total_points; i++)
          {
             delete[] dataY[j][i];
          }
          delete[] dataY[j];
        }
        delete[] dataY;
      }
      else
      {
        for (int j = 0; j < 1; j++)
        {
          for (int i = 0; i < total_points; i++)
          {
            delete[] dataY_reg[j][i];
          }
          delete[] dataY_reg[j];
        }
        delete[] dataY_reg;
      }
    }
  }
}
void Node::deleteData()
{
  if (!regression)
  {
    delete[] Ni;
  }
  if (hasData && possesses_data)
  {
    possesses_data = false;
    hasData = false;
    if (total_numeric_features > 0)
    {
      if (!regression)
      {
        for (int j = 0; j < total_numeric_features; j++)
        {
          for (int i = 0; i < total_points; i++)
          {
            delete[] dataY[j][i];
          }
          delete[] dataY[j];
        }
        delete[] dataY;
      }
      else
      {
        for (int j = 0; j < total_numeric_features; j++)
        {
          for (int i = 0; i < total_points; i++)
          {
            delete[] dataY_reg[j][i];
          }
          delete[] dataY_reg[j];
        }
        delete[] dataY_reg;
      }
    }
    else
    {
      if (!regression)
      {
        for (int j = 0; j < 1; j++)
        {
          for (int i = 0; i < total_points; i++)
          {
            delete[] dataY[j][i];
          }
          delete[] dataY[j];
        }
        delete[] dataY;
      }
      else
      {
        for (int j = 0; j < 1; j++)
        {
          for (int i = 0; i < total_points; i++)
          {
            delete[] dataY_reg[j][i];
          }
          delete[] dataY_reg[j];
        }
        delete[] dataY_reg;
      }
    }
  }
}
Node::~Node()
{
  deleteData();
}
vector<Node *> Node::distribute(vector<Node *> nodes)
{
  if (nodes.size() == 0)
  {
    return nodes;
  }
  vector<int> which_child;
  int *total_pointsS = new int[nodes.size()];
  for (int i = 0; i < nodes.size(); i++)
  {
    total_pointsS[i] = 0;
  }
  vector<double> x;
  vector<double> x_nominal;
  map<int, int> which_child2;
  for (int i = 0; i < total_points; i++)
  {
    x.clear();
    x_nominal.clear();
    int ind = 0;
    if (!regression)
    {
      for (int j = 0; j < total_numeric_features; j++)
      {
        x.push_back(X[j][dataY[0][i][0]]);
      }
      for (int j = 0; j < total_nominal_features; j++)
      {
        x_nominal.push_back(Xnominal[dataY[0][i][0]][j]);
      }
    }
    else
    {
      for (int j = 0; j < total_numeric_features; j++)
      {
        x.push_back(X[j][(int)dataY_reg[0][i][0]]);
      }
      for (int j = 0; j < total_nominal_features; j++)
      {
        x_nominal.push_back(Xnominal[(int)dataY_reg[0][i][0]][j]);
      }
    }
    for (int j = 0; j < nodes.size(); j++)
    {
      if (nodes[j]->evaluate(x, x_nominal))
      {
        ind = j;
        if (!regression)
        {
          which_child2[dataY[0][i][0]] = ind;
        }
        else
        {
          which_child2[(int)dataY_reg[0][i][0]] = ind;
        }
        break;
      }
    }
    which_child.push_back(ind);
    total_pointsS[ind] += 1;
  }
  bool empty_node = false;
  for (int i = 0; i < nodes.size(); i++)
  {
    if (total_pointsS[i] == 0)
    {
      empty_node = true;
    }
  }
  vector<int ***> data;
  vector<double ***> data_reg;
  if (!regression)
  {
    for (int l = 0; l < nodes.size(); l++)
    {
      int ***data0 = new int **[max(1, total_numeric_features)];
      for (int j = 0; j < max(1, total_numeric_features); j++)
      {
        data0[j] = new int *[total_pointsS[l]];
        for (int i = 0; i < total_pointsS[l]; i++)
        {
          data0[j][i] = new int[2];
        }
      }
      data.push_back(data0);
      data_reg.push_back(NULL);
    }
    int **cs = new int *[nodes.size()];
    for (int i = 0; i < nodes.size(); i++)
    {
      cs[i] = new int[max(1, total_numeric_features)];
      for (int j = 0; j < max(1, total_numeric_features); j++)
      {
        cs[i][j] = 0;
      }
    }
    for (int j = 0; j < max(1, total_numeric_features); j++)
    {
      for (int i = 0; i < total_points; i++)
      {
        int ind = which_child2[dataY[j][i][0]];
        data[ind][j][cs[ind][j]][0] = dataY[j][i][0];
        data[ind][j][cs[ind][j]][1] = dataY[j][i][1];
        cs[ind][j]++;
      }
    }
    for (int i = 0; i < nodes.size(); i++)
    {
      delete[] cs[i];
    }
    delete[] cs;
  }
  else
  {
    for (int l = 0; l < nodes.size(); l++)
    {
      double ***data0 = new double **[max(1, total_numeric_features)];
      for (int j = 0; j < max(1, total_numeric_features); j++)
      {
        data0[j] = new double *[total_pointsS[l]];
        for (int i = 0; i < total_pointsS[l]; i++)
        {
          data0[j][i] = new double[2];
        }
      }
      data_reg.push_back(data0);
      data.push_back(NULL);
    }
    int **cs = new int *[nodes.size()];
    for (int i = 0; i < nodes.size(); i++)
    {
      cs[i] = new int[max(1, total_numeric_features)];
      for (int j = 0; j < max(1, total_numeric_features); j++)
      {
        cs[i][j] = 0;
      }
    }
    for (int j = 0; j < max(1, total_numeric_features); j++)
    {
      for (int i = 0; i < total_points; i++)
      {
        int ind = 0;
        if (!regression)
        {
          ind = which_child2[dataY[j][i][0]];
        }
        else
        {
          ind = which_child2[(int)dataY_reg[j][i][0]];
        }
        data_reg[ind][j][cs[ind][j]][0] = dataY_reg[j][i][0];
        data_reg[ind][j][cs[ind][j]][1] = dataY_reg[j][i][1];
        cs[ind][j]++;
      }
    }
    for (int i = 0; i < nodes.size(); i++)
    {
      delete[] cs[i];
    }
    delete[] cs;
  }
  for (int l = 0; l < nodes.size(); l++)
  {
    nodes[l]->deleteData();
    nodes[l]->dataY = data[l];
    nodes[l]->dataY_reg = data_reg[l];
    nodes[l]->total_points = total_pointsS[l];
    nodes[l]->hasData = true;
    nodes[l]->possesses_data = true;
    nodes[l]->getData();
  }
  delete[] total_pointsS;
  return nodes;
}
vector<Node *> Node::divide(Rule rule)
{
  vector<Node *> nodes;
  if (!rule.split_found)
  {
    return nodes;
  }
  beautify(rule, env->max_margin);
  rule.normalize();
  vector<int> which_child;
  which_child.reserve(total_points);
  map<int, int> which_child2;
  int *total_pointsS = new int[rule.no_children];
  for (int i = 0; i < rule.no_children; i++)
  {
    total_pointsS[i] = 0;
  }
  for (int i = 0; i < total_points; i++)
  {
    int ind = 0;
    if (!regression)
    {
      ind = rule.evaluate(dataY[0][i][0], X, Xnominal);
      which_child2[dataY[0][i][0]] = ind;
    }
    else
    {
      ind = rule.evaluate((int)dataY_reg[0][i][0], X, Xnominal);
      which_child2[(int)dataY_reg[0][i][0]] = ind;
    }
    which_child.push_back(ind);
    total_pointsS[ind] += 1;
  }
  bool empty_node = false;
  int empty_nodes = 0;
  for (int i = 0; i < rule.no_children; i++)
  {
    if (total_pointsS[i] == 0)
    {
      empty_node = true;
      empty_nodes++;
    }
  }
  if (empty_nodes >= rule.no_children - 1)
  {
    delete[] total_pointsS;
    return nodes;
  }
  if (empty_node && rule.no_children == 2)
  {
    delete[] total_pointsS;
    return nodes;
  }
  if (empty_nodes >= rule.no_children - 2 && rule.is_cross)
  {
    delete[] total_pointsS;
    return divide(chooseBetter(rule));
  }
  vector<int ***> data;
  vector<double ***> data_reg;
  if (!regression)
  {
    for (int l = 0; l < rule.no_children; l++)
    {
      int ***data0 = new int **[max(1, total_numeric_features)];
      for (int j = 0; j < max(1, total_numeric_features); j++)
      {
        data0[j] = new int *[total_pointsS[l]];
        for (int i = 0; i < total_pointsS[l]; i++)
        {
          data0[j][i] = new int[2];
        }
      }
      data.push_back(data0);
      data_reg.push_back(NULL);
    }
    int **cs = new int *[rule.no_children];
    for (int i = 0; i < rule.no_children; i++)
    {
      cs[i] = new int[max(1, total_numeric_features)];
      for (int j = 0; j < max(1, total_numeric_features); j++)
      {
        cs[i][j] = 0;
      }
    }
    for (int j = 0; j < max(1, total_numeric_features); j++)
    {
      for (int i = 0; i < total_points; i++)
      {
        int ind = which_child2[dataY[j][i][0]];
        data[ind][j][cs[ind][j]][0] = dataY[j][i][0];
        data[ind][j][cs[ind][j]][1] = dataY[j][i][1];
        cs[ind][j]++;
      }
    }
    for (int i = 0; i < rule.no_children; i++)
    {
      delete[] cs[i];
    }
    delete[] cs;
  }
  else
  {
    for (int l = 0; l < rule.no_children; l++)
    {
      double ***data0 = new double **[max(1, total_numeric_features)];
      for (int j = 0; j < max(1, total_numeric_features); j++)
      {
        data0[j] = new double *[total_pointsS[l]];
        for (int i = 0; i < total_pointsS[l]; i++)
        {
          data0[j][i] = new double[2];
        }
      }
      data_reg.push_back(data0);
      data.push_back(NULL);
    }
    int **cs = new int *[rule.no_children];
    for (int i = 0; i < rule.no_children; i++)
    {
      cs[i] = new int[max(1, total_numeric_features)];
      for (int j = 0; j < max(1, total_numeric_features); j++)
      {
        cs[i][j] = 0;
      }
    }
    for (int j = 0; j < max(1, total_numeric_features); j++)
    {
      for (int i = 0; i < total_points; i++)
      {
        int ind = which_child2[(int)dataY_reg[j][i][0]];
        ind = rule.evaluate((int)dataY_reg[j][i][0], X, Xnominal);
        data_reg[ind][j][cs[ind][j]][0] = dataY_reg[j][i][0];
        data_reg[ind][j][cs[ind][j]][1] = dataY_reg[j][i][1];
        cs[ind][j]++;
      }
    }
    for (int i = 0; i < rule.no_children; i++)
    {
      delete[] cs[i];
    }
    delete[] cs;
  }
  if (!rule.is_nominal && !rule.is_cross)
  {
    if (rule.leq)
    {
      Node *node = new Node(X, Xnominal, Y, Y_reg, data[1], data_reg[1], total_pointsS[1], total_numeric_features, total_nominal_features, nominal_domain, is_nominal_val_in_dataset, n_classes, p, env);
      node->cond = Condition(rule.a, rule.b, Condition::leq);
      nodes.push_back(node);
      node->hasData = true;
      node->possesses_data = true;
      node = new Node(X, Xnominal, Y, Y_reg, data[0], data_reg[0], total_pointsS[0], total_numeric_features, total_nominal_features, nominal_domain, is_nominal_val_in_dataset, n_classes, p, env);
      node->cond = Condition(rule.a, rule.b, Condition::g);
      nodes.push_back(node);
      node->hasData = true;
      node->possesses_data = true;
    }
    else
    {
      Node *node = new Node(X, Xnominal, Y, Y_reg, data[0], data_reg[0], total_pointsS[0], total_numeric_features, total_nominal_features, nominal_domain, is_nominal_val_in_dataset, n_classes, p, env);
      node->cond = Condition(rule.a, rule.b, Condition::geq);
      nodes.push_back(node);
      node->hasData = true;
      node->possesses_data = true;
      node = new Node(X, Xnominal, Y, Y_reg, data[1], data_reg[1], total_pointsS[1], total_numeric_features, total_nominal_features, nominal_domain, is_nominal_val_in_dataset, n_classes, p, env);
      node->cond = Condition(rule.a, rule.b, Condition::l);
      nodes.push_back(node);
      node->hasData = true;
      node->possesses_data = true;
    }
  }
  else if (rule.is_nominal)
  {
    bool reduce_domain = true;
    for (int l = 0; l < rule.no_children; l++)
    {
      vector<vector<bool>> is_nominal_val_in_dataset2;
      for (int j = 0; j < total_nominal_features; j++)
      {
        vector<bool> d(is_nominal_val_in_dataset[j].size(), false);
        is_nominal_val_in_dataset2.push_back(d);
        if (j != rule.nominal_ind)
        {
          for (int k = 0; k < is_nominal_val_in_dataset[j].size(); k++)
          {
            is_nominal_val_in_dataset2[j][k] = is_nominal_val_in_dataset[j][k];
          }
        }
      }
      vector<bool> in_partition;
      for (int k = 0; k <= nominal_domain[rule.nominal_ind]; k++)
      {
        if (rule.which_partition[k] == l)
        {
          in_partition.push_back(true);
          if(is_nominal_val_in_dataset[rule.nominal_ind][k])
          {
            is_nominal_val_in_dataset2[rule.nominal_ind][k] = true;
          }
        }
        else
        {
          in_partition.push_back(false);
        }
      }
      if (!reduce_domain)
      {
        for (int k = 0; k < is_nominal_val_in_dataset[rule.nominal_ind].size(); k++)
        {
          is_nominal_val_in_dataset2[rule.nominal_ind][k] = is_nominal_val_in_dataset[rule.nominal_ind][k];
        }
      }
      Node *node = new Node(X, Xnominal, Y, Y_reg, data[l], data_reg[l], total_pointsS[l], total_numeric_features, total_nominal_features, nominal_domain, is_nominal_val_in_dataset2, n_classes, p, env);
      node->cond = Condition(rule.nominal_ind, in_partition, is_nominal_val_in_dataset2);
      node->hasData = true;
      node->possesses_data = true;
      nodes.push_back(node);
    }
  }
  else
  {
    Node *node = new Node(X, Xnominal, Y, Y_reg, data[0], data_reg[0], total_pointsS[0], total_numeric_features, total_nominal_features, nominal_domain, is_nominal_val_in_dataset, n_classes, p, env);
    node->cond = Condition(rule.ind1, rule.ind2, rule.x1, rule.x2, Condition::leq, Condition::leq);
    nodes.push_back(node);
    node->hasData = true;
    node->possesses_data = true;
    node = new Node(X, Xnominal, Y, Y_reg, data[1], data_reg[1], total_pointsS[1], total_numeric_features, total_nominal_features, nominal_domain, is_nominal_val_in_dataset, n_classes, p, env);
    node->cond = Condition(rule.ind1, rule.ind2, rule.x1, rule.x2, Condition::g, Condition::leq);
    nodes.push_back(node);
    node->hasData = true;
    node->possesses_data = true;
    node = new Node(X, Xnominal, Y, Y_reg, data[2], data_reg[2], total_pointsS[2], total_numeric_features, total_nominal_features, nominal_domain, is_nominal_val_in_dataset, n_classes, p, env);
    node->cond = Condition(rule.ind1, rule.ind2, rule.x1, rule.x2, Condition::leq, Condition::g);
    nodes.push_back(node);
    node->hasData = true;
    node->possesses_data = true;
    node = new Node(X, Xnominal, Y, Y_reg, data[3], data_reg[3], total_pointsS[3], total_numeric_features, total_nominal_features, nominal_domain, is_nominal_val_in_dataset, n_classes, p, env);
    node->cond = Condition(rule.ind1, rule.ind2, rule.x1, rule.x2, Condition::g, Condition::g);
    nodes.push_back(node);
    node->hasData = true;
    node->possesses_data = true;
  }
  delete[] total_pointsS;
  for (int i = 0; i < nodes.size(); i++)
  {
    if (nodes[i]->total_points == 0)
    {
      if (env->regression)
      {
        nodes[i]->median = median;
        nodes[i]->mean = mean;
        nodes[i]->mse = 0;
        nodes[i]->mae = 0;
      }
      else
      {
        nodes[i]->maxVote = maxVote;
        nodes[i]->impurity = 0;
      }
    }
  }
  return nodes;
}
vector<Node *> Node::branch_randomly()
{
  vector<Node *> nodes;
  vector<pair<int, int>> features;
  for (int j = 0; j < total_numeric_features; j++)
  {
    features.push_back(pair<int, int>(j, 0));
  }
  for (int j = 0; j < total_nominal_features; j++)
  {
    features.push_back(pair<int, int>(j, 1));
  }
  shuffle(features.begin(), features.end(), env->generator);
  for (int j2 = 0; j2 < features.size(); j2++)
  {
    bool split_found = false;
    if (features[j2].second == 0)
    {
      int j = features[j2].first;
      vector<double> split_vals;
      if (!regression)
      {
        if (X[j][dataY[j][total_points - 1][0]] != X[j][dataY[j][0][0]])
        {
          for (int i = 1; i < total_points; i++)
          {
            if (X[j][dataY[j][i][0]] != X[j][dataY[j][i - 1][0]])
            {
              split_vals.push_back((X[j][dataY[j][i][0]] + X[j][dataY[j][i - 1][0]]) / 2.0);
            }
          }
        }
      }
      else
      {
        if (X[j][(int)dataY_reg[j][total_points - 1][0]] != X[j][(int)dataY_reg[j][0][0]])
        {
          for (int i = 1; i < total_points; i++)
          {
            if (X[j][(int)dataY_reg[j][i][0]] != X[j][(int)dataY_reg[j][i - 1][0]])
            {
              split_vals.push_back((X[j][(int)dataY_reg[j][i][0]] + X[j][(int)dataY_reg[j][i - 1][0]]) / 2.0);
            }
          }
        }
      }
      if (split_vals.size() > 0)
      {
        shuffle(split_vals.begin(), split_vals.end(), env->generator);
        double b = split_vals[0];
        vector<double> a(total_numeric_features, 0);
        a[j] = 1;
        Rule rule(a, b);
        rule.leq = true;
        rule.split_found = true;
        rule.imp = 1;
        nodes = divide(rule);
        impS = 0;
        for (int l = 0; l < nodes.size(); l++)
        {
          if (!regression)
          {
            impS += nodes[l]->impurity;
          }
          else
          {
            impS += nodes[l]->total_points * nodes[l]->mse;
          }
        }
        return nodes;
      }
    }
    else
    {
      int j = features[j2].first;
      int N = nominal_domain[j];
      int L = env->nominal_partitions;
      vector<bool> found(N + 1, false);
      vector<int> values;
      for (int i = 0; i < total_points; i++)
      {
        if (!regression)
        {
          int v = (int)Xnominal[dataY[0][i][0]][j];
          if (!found[v])
          {
            values.push_back(v);
          }
          found[v] = true;
        }
        else
        {
          int v = (int)Xnominal[(int)dataY_reg[0][i][0]][j];
          if (!found[v])
          {
            values.push_back(v);
          }
          found[v] = true;
        }
      }
      if (values.size() >= 2)
      {
        vector<int> which_partition(N + 1, 0);
        bool splits = false;
        int next_partition = 0;
        while (!splits)
        {
          next_partition = 0;
          shuffle(values.begin(), values.end(), env->generator);
          for (int v = 0; v < values.size(); v++)
          {
            which_partition[v] = env->rand() % min(L, next_partition + 1);
            if (which_partition[v] == next_partition)
            {
              next_partition++;
            }
          }
          if (next_partition >= 2)
          {
            splits = true;
          }
        }
        for (int k = 0; k <= N; k++)
        {
          if (!found[k])
          {
            which_partition[k] = env->rand() % (min(L, next_partition));
          }
        }
        Rule rule(j, which_partition);
        rule.split_found = true;
        rule.imp = 1;
        nodes = divide(rule);
        impS = 0;
        for (int l = 0; l < nodes.size(); l++)
        {
          if (!regression)
          {
            impS += nodes[l]->impurity;
          }
          else
          {
            impS += nodes[l]->total_points * nodes[l]->mse;
          }
        }
        return nodes;
      }
    }
    if (split_found)
    {
      break;
    }
  }
  return nodes;
}
vector<Node *> Node::branch()
{
  if (env->branch_randomly)
  {
    return branch_randomly();
  }
  Rule bestRule;
  Rule bestSingle;
  Rule bestBi;
  Rule bestMult;
  Rule bestCat;
  bool has_bivariate = false;
  bool has_single = false;
  int max_nom_features = 0;
  int max_num_features = 0;
  if (env->max_features >= total_features)
  {
    max_nom_features = total_nominal_features;
    max_num_features = total_numeric_features;
  }
  else
  {
    vector<int> features(total_features, 0);
    for (int j = 0; j < total_nominal_features; j++)
    {
      features[j] = 1;
    }
    shuffle(features.begin(), features.end(), env->generator);
    for (int j = 0; j < env->max_features; j++)
    {
      if (features[j] == 1)
      {
        max_nom_features += 1;
      }
      else
      {
        max_num_features += 1;
      }
    }
  }
  if (total_points == 0 || (!regression && impurity == 0.0) || (regression && env->criterion == Environment::mse && mse == 0.0) || (regression && env->criterion == Environment::mae && mae == 0.0))
  {
    vector<Node *> nodes;
    return nodes;
  }
  if (total_numeric_features > 0)
  {
    vector<int> features;
    for (int j = 0; j < total_numeric_features; j++)
    {
      features.push_back(j);
    }
    shuffle(features.begin(), features.end(), env->generator);
    if (env->univariate)
    {
      bestSingle = branchUnivariate(features, max_num_features);
      if (bestSingle.split_found)
      {
        has_single = true;
      } 
    }
    else
    {
      bestSingle.attribute_order = features;
      for (int j = 0; j < total_features; j++)
      {
        bestSingle.a.push_back(0);
      }
    }
    bestRule = bestSingle;
    if (env->cross)
    {
      Rule bestCross = branchCross(bestSingle, max_num_features);
      if (bestCross.split_found)
      {
        if (true)
        {
          if (env->lookahead)
          {
            bestRule = chooseBetter(bestCross);
          }
          else
          {
            bestRule = bestCross;
          }
        }
      }
    }
    if ((bestRule.imp != 0 && env->bi && total_points >= env->oblique_split_threshold && max_num_features >= 2))
    {
      has_bivariate = true;
      bestBi = branchBivariate(bestSingle, max_num_features);
    }
    if (bestBi.imp < bestRule.imp)
    {
      bestRule = bestBi;
    }
    if (bestRule.imp != 0 && env->mult && total_points >= env->oblique_split_threshold && max_num_features >= 2)
    {
      if (has_bivariate && bestBi.imp < bestSingle.imp)
      {
        bestMult = branchMultivariate(bestBi, max_num_features);
      }
      else
      {
        bestMult = branchMultivariate(bestSingle, max_num_features);
      }
    }
    if (bestMult.imp < bestRule.imp)
    {
      bestRule = bestMult;
    }
  }
  if (total_nominal_features >= 1)
  {
    vector<int> features(total_nominal_features, 0);
    for (int j = 0; j < total_nominal_features; j++)
    {
      features[j] = j;
    }
    bestCat = branchNominal(features, max_nom_features);
  }
  if (bestCat.split_found && (!bestRule.split_found || bestCat.imp < bestRule.imp))
  {
    bestRule = bestCat;
  };
  vector<Node *> nodes;
  nodes = divide(bestRule);
  if (nodes.size() == 0 && has_single && !bestRule.is_cross)
  {
    cout<<"id: " << "No split found, using univariate instead.\n";
    nodes = divide(bestSingle);
  }
  if (nodes.size() == 0 && !bestRule.is_cross)
  {
    cout<<"id: " << "No split found. Trying univariate instead.\n";
    if (total_numeric_features > 0)
    {
      vector<int> features;
      for (int j = 0; j < total_numeric_features; j++)
      {
        features.push_back(j);
      }
      bestSingle = branchUnivariate(features, max_num_features);
      if (bestSingle.split_found)
      {
        has_single = true;
      }
    }
    nodes = divide(bestSingle);
    if (nodes.size() == 0)
    {
      cout << id << ": "
           << "... but no split found. Giving up on that path!\n";
    }
  }
  impS = bestRule.imp;
  return nodes;
}
Rule Node::branchNominal(vector<int> features, int max_features)
{
  shuffle(features.begin(), features.end(), env->generator);
  Rule *rule3 = new Rule();
  if (!env->regression)
  {
    for (int j2 = 0; j2 < features.size(); j2++)
    {
      if (rule3->split_found && j2 >= max_features)
      {
        continue;
      }
      int j = features[j2];
      NBranchAndBound_class(Xnominal, Y, dataY, total_points, n_classes, nominal_domain[j], j, is_nominal_val_in_dataset, rule3, env, p);
    }
  }
  else
  {
    for (int j2 = 0; j2 < features.size(); j2++)
    {
      if (rule3->split_found && j2 >= max_features)
      {
        continue;
      }
      int j = features[j2];
      NBranchAndBound_reg(Xnominal, Y_reg, dataY_reg, total_points, nominal_domain[j], j, is_nominal_val_in_dataset, rule3, env);
    }
  }
  opt_gap = rule3->opt_gap;
  time = rule3->time;
  if (!rule3->split_found)
  {
    delete rule3;
    Rule rule = Rule();
    rule.split_found = false;
    return rule;
  }
  int ind = rule3->nominal_ind;
  vector<int> which_partition;
  for (int k = 0; k <= nominal_domain[ind]; k++)
  {
    which_partition.push_back(0);
  }
  for (int l = 0; l < rule3->partition.size(); l++)
  {
    for (int i2 = 0; i2 < rule3->partition[l].size(); i2++)
    {
      which_partition[rule3->partition[l][i2]] = l;
    }
  }
  Rule rule(ind, which_partition);
  rule.imp = rule3->imp;
  rule.split_found = true;
  delete rule3;
  return rule;
}
Rule Node::branchUnivariate(vector<int> features2, int max_features)
{
  int *N1;
  int *N2;
  int *N;
  if (!env->regression)
  {
    N1 = new int[n_classes];
    N2 = new int[n_classes];
    N = new int[n_classes];
  }
  double bestImp;
  double bestImp2;
  if (!env->regression)
  {
    bestImp = impurity + 1;
    bestImp2 = impurity + 1;
  }
  else
  {
    if (env->criterion == Environment::mse)
    {
      bestImp = total_points * mse + 1;
    }
    else
    {
      bestImp = a_error + 1;
    }
    bestImp2 = std::numeric_limits<double>::max();
  }
  int bestInd = 0;
  double bestX = 0;
  bool split_found = false;
  bool split_possible = false;
  vector<pair<int, double>> best_splits;
  vector<double> bestXs;
  for (int j = 0; j < total_numeric_features; j++)
  {
    bestXs.push_back(0);
  }
  if (!env->regression)
  {
    for (int m = 0; m < n_classes; m++)
    {
      N1[m] = 0;
      N[m] = 0;
      N2[m] = 0;
    }
    for (int i = 0; i < total_points; i++)
    {
      N[dataY[0][i][1]] += 1;
      N2[dataY[0][i][1]] += 1;
    }
  }
  vector<pair<int, double>> features;
  vector<bool> single_considered(total_numeric_features);
  for (int j = 0; j < single_considered.size(); j++)
  {
    single_considered[j] = false;
  }
  vector<Vector_reg> reg;
  for (int j2 = 0; j2 < features2.size(); j2++)
  {
    double bound;
    if (!env->regression)
    {
      bestImp2 = impurity + 1;
      bound = impurity;
    }
    else
    {
      if (env->criterion == Environment::mse)
      {
        bestImp2 = total_points * mse + 1;
        bound = total_points * mse;
      }
      else
      {
        bestImp2 = a_error + 1;
        bound = a_error;
      }
    }
    int j = features2[j2];
    double val = total_points;
    if (j2 >= max_features && bestImp < bound)
    {
      features.push_back(std::pair<int, double>(j, bestImp2 + j2 * numeric_limits<double>::min()));
      single_considered[j] = false;
      continue;
    }
    single_considered[j] = true;
    if (!env->regression)
    {
      for (int m = 0; m < n_classes; m++)
      {
        N1[m] = 0;
        N2[m] = N[m];
      }
    }
    else
    {
      reg.clear();
      reg.push_back(Vector_reg());
      reg.push_back(Vector_reg());
      reg[0].reserve(total_points);
      reg[1].reserve(total_points);
      if (env->criterion == Environment::mae)
      {
        reg[0].use_mae();
        reg[1].use_mae();
      }
      for (int i = total_points - 1; i >= 0; i--)
      {
        reg[1].push_back(dataY_reg[j][i][1]);
      }
    }
    int right = total_points;
    bool split_on_feat_possible = false;
    double x0;
    if (!env->regression)
    {
      x0 = X[j][dataY[j][0][0]];
    }
    else
    {
      x0 = X[j][(int)dataY_reg[j][0][0]];
    }
    for (int i = 0; i < total_points;)
    {
      double x;
      if (!env->regression)
      {
        x = X[j][dataY[j][i][0]];
      }
      else
      {
        x = X[j][(int)dataY_reg[j][i][0]];
      }
      if (x != x0)
      {
        split_possible = true;
        split_on_feat_possible = true;
      }
      if (!env->regression)
      {
        while (i < total_points && X[j][dataY[j][i][0]] == x)
        {
          right--;
          N1[dataY[j][i][1]] += 1;
          N2[dataY[j][i][1]] -= 1;
          i++;
        }
      }
      else
      {
        while (i < total_points && X[j][(int)dataY_reg[j][i][0]] == x)
        {
          right--;
          reg[0].push_back(dataY_reg[j][i][1]);
          reg[1].pop_back();
          i++;
        }
      }
      if (right == 0)
      {
        break;
      }
      if (!env->regression)
      {
        val = imp(N1, N2, n_classes, env->criterion);
      }
      else
      {
        val = reg[0].error(env->criterion) + reg[1].error(env->criterion);
      }
      if (val <= bestImp2 && i < total_points)
      {
        bestImp2 = val;
        if (!env->regression)
        {
          if (i < total_points)
          {
            bestXs[j] = (X[j][dataY[j][i - 1][0]] + X[j][dataY[j][i][0]]) / 2;
          }
          else
          {
            bestXs[j] = (X[j][dataY[j][i - 1][0]]);
          }
        }
        else
        {
          if (i < total_points)
          {
            bestXs[j] = (X[j][(int)dataY_reg[j][i - 1][0]] + X[j][(int)dataY_reg[j][i][0]]) / 2;
          }
          else
          {
            bestXs[j] = (X[j][dataY[j][i - 1][0]]);
          }
        }
      }
      if (val <= bestImp && i < total_points)
      {
        split_found = true;
        bestInd = j;
        bestX = x;
        if (!env->regression)
        {
          if (i < total_points)
          {
            bestX = (X[j][dataY[j][i - 1][0]] + X[j][dataY[j][i][0]]) / 2;
          }
          else
          {
            bestX = (X[j][dataY[j][i - 1][0]]);
          }
        }
        else
        {
          if (i < total_points)
          {
            bestX = (X[j][(int)dataY_reg[j][i - 1][0]] + X[j][(int)dataY_reg[j][i][0]]) / 2;
          }
          else
          {
            bestX = (X[j][(int)dataY_reg[j][i - 1][0]]);
          }
        }
        if (val < bestImp)
        {
          best_splits.clear();
        }
        best_splits.push_back(std::pair<int, double>(j, bestX));
        bestImp = val;
      }
    }
    features.push_back(std::pair<int, double>(j, bestImp2 + j * numeric_limits<double>::min()));
  }
  Rule rule;
  if (best_splits.size() > 0 && split_possible)
  {
    std::sort(features.begin(), features.end(), compareP);
    shuffle(best_splits.begin(), best_splits.end(), env->generator);
    bestInd = best_splits[0].first;
    bestX = best_splits[0].second;
    vector<double> a;
    for (int j = 0; j < total_numeric_features; j++)
    {
      a.push_back(0);
    }
    a[bestInd] = 1;
    double b = bestX;
    Rule rule_tmp(a, b);
    rule = rule_tmp;
    rule.leq = true;
    rule.imp = bestImp;
    rule.split_found = true;
    vector<bool> flags(features.size());
    for (int j = 0; j < features.size(); j++)
    {
      flags[j] = false;
    }
    for (int j = 0; j < features.size(); j++)
    {
      if (!flags[features[j].first])
      {
        rule.attribute_order.push_back(features[j].first);
      }
    }
    rule.single_considered = single_considered;
    if (features2.size() > 1)
    {
      rule.ind1 = features[0].first;
      rule.x1 = bestXs[rule.ind1];
      rule.ind2 = features[1].first;
      rule.x2 = bestXs[rule.ind2];
    }
    else
    {
      rule.ind1 = features[0].first;
      rule.x1 = bestXs[rule.ind1];
      rule.ind2 = features[0].first;
      rule.x2 = bestXs[rule.ind2];
    }
  }
  else
  {
    rule.imp = std::numeric_limits<double>::max();
    rule.split_found = false;
  }
  if (!env->regression)
  {
    delete[] N1;
    delete[] N2;
    delete[] N;
  }
  return rule;
}
Rule Node::branchBivariate(vector<int> features, double bound)
{
  bool split_found = false;
  Rule rule;
  rule.imp = std::numeric_limits<double>::max();
  rule.split_found = split_found;
  Result *res = new Result(split_found, bound, 0.0, rule.x1, rule.ind2, rule.ind1, total_numeric_features);
  vector<future<void>> fs;
  vector<bool> all_equal(features.size());
  int usable_features = 0;
  for (int j2 = 0; j2 < features.size(); j2++)
  {
    int j = features[j2];
    all_equal[j2] = true;
    if (!regression)
    {
      double x0 = X[j][dataY[j][0][0]];
      for (int i = total_points - 1; i < total_points; i++)
      {
        double x = X[j][dataY[j][i][0]];
        if (x != x0)
        {
          all_equal[j2] = false;
          usable_features++;
          break;
        }
      }
    }
    else
    {
      double x0 = X[j][(int)dataY_reg[j][0][0]];
      for (int i = total_points - 1; i < total_points; i++)
      {
        double x = X[j][(int)dataY_reg[j][i][0]];
        if (x != x0)
        {
          all_equal[j2] = false;
          usable_features++;
          break;
        }
      }
    }
  }
  if (usable_features == 0)
  {
    return rule;
  }
  if (usable_features == 1)
  {
    return branchUnivariate(features, features.size());
  }
  vector<bool> considered(features.size());
  for (int j = 0; j < features.size(); j++)
  {
    considered[j] = false;
  }
  int bisplits = 0;
  int bisplits2 = 0;
  int not_bisplits = 0;
  for (int k = 0; k < features.size(); k++)
  {
    int j1 = features[k];
    for (int l = k + 1; l < features.size(); l++)
    {
      bisplits2++;
      int j2 = features[l];
      if ((all_equal[k] && all_equal[l]))
      {
        not_bisplits++;
        continue;
      }
      if (all_equal[k] && considered[l] == true)
      {
        not_bisplits++;
        continue;
      }
      if (all_equal[l] && considered[k] == true)
      {
        not_bisplits++;
        continue;
      }
      if (!all_equal[k] && !all_equal[l])
      {
        bisplits++;
      }
      considered[k] = true;
      considered[l] = true;
      if (!regression)
      {
        if (env->n_threads > 1)
        {
          fs.push_back(p->push(branchAndBound, X[j1], X[j2], j1, j2, n_classes, total_points, dataY[j1], dataY[j2], (double **)NULL, (double **)NULL, res, env));
          fs.push_back(p->push(branchAndBound, X[j2], X[j1], j2, j1, n_classes, total_points, dataY[j2], dataY[j1], (double **)NULL, (double **)NULL, res, env));
        }
        else
        {
          branchAndBound(0, X[j1], X[j2], j1, j2, n_classes, total_points, dataY[j1], dataY[j2], (double **)NULL, (double **)NULL, res, env);
          branchAndBound(0, X[j2], X[j1], j2, j1, n_classes, total_points, dataY[j2], dataY[j1], (double **)NULL, (double **)NULL, res, env);
        }
      }
      else
      {
        if (env->n_threads > 1)
        {
          fs.push_back(p->push(branchAndBound, X[j1], X[j2], j1, j2, n_classes, total_points, (int **)NULL, (int **)NULL, dataY_reg[j1], dataY_reg[j2], res, env));
          fs.push_back(p->push(branchAndBound, X[j2], X[j1], j2, j1, n_classes, total_points, (int **)NULL, (int **)NULL, dataY_reg[j2], dataY_reg[j1], res, env));
        }
        else
        {
          branchAndBound(0, X[j1], X[j2], j1, j2, n_classes, total_points, (int **)NULL, (int **)NULL, dataY_reg[j1], dataY_reg[j2], res, env);
          branchAndBound(0, X[j2], X[j1], j2, j1, n_classes, total_points, (int **)NULL, (int **)NULL, dataY_reg[j2], dataY_reg[j1], res, env);
        }
      }
    }
  }
  for (int k = 0; k < fs.size(); k++)
  {
    fs[k].wait();
  }
  if (!res->split_found)
  {
    cout << "No split found\n";
    delete res;
    rule.split_found = false;
    return rule;
  }
  vector<double> a;
  double b = res->bestB;
  for (int j = 0; j < total_numeric_features; j++)
  {
    if (j == res->bestInd1)
    {
      if (j == res->bestInd2)
      {
        a.push_back(1);
        b = res->bestB / (1 - res->bestA);
      }
      else
      {
        a.push_back(-res->bestA);
      }
    }
    else if (j == res->bestInd2)
    {
      a.push_back(1);
    }
    else
    {
      a.push_back(0);
    }
  }
  bool to_leq = true;
  int j1 = res->bestInd1;
  int j2 = res->bestInd2;
  bool has_eq = false;
  double min_b = std::numeric_limits<double>::lowest();
  double max_b = std::numeric_limits<double>::max();
  bool divides = false;
  bool upper = false;
  bool lower = false;
  for (int i = 0; i < total_points; i++)
  {
    double tmp;
    if (!regression)
    {
      tmp = a[j1] * X[j1][dataY[0][i][0]] + a[j2] * X[j2][dataY[0][i][0]];
    }
    else
    {
      tmp = a[j1] * X[j1][(int)dataY_reg[0][i][0]] + a[j2] * X[j2][(int)dataY_reg[0][i][0]];
    }
    if (tmp == b)
    {
      to_leq = false;
      has_eq = true;
    }
    if (tmp >= b)
    {
      upper = true;
      if (tmp < max_b)
      {
        max_b = tmp;
      }
    }
    else
    {
      lower = true;
      if (tmp > min_b)
      {
        min_b = tmp;
      }
    }
  }
  if (upper && lower)
  {
    divides = true;
  }
  double b_orig = b;
  if (min_b > std::numeric_limits<double>::lowest() && max_b < std::numeric_limits<double>::max() && max_b > min_b)
  {
    to_leq = true;
    b = (min_b + max_b) / 2.;
  }
  if (divides)
  {
  }
  Rule rule2(a, b);
  rule.leq = to_leq;
  rule2.split_found = res->split_found;
  rule2.imp = res->bestImp;
  rule2.divides = res->divides;
  delete res;
  return rule2;
}
Rule Node::branchBivariate(Rule rule, int max_features)
{
  vector<int> features;
  if (rule.attribute_order.size() > 0)
  {
    features = rule.attribute_order;
  }
  else
  {
    for (int j = 0; j < total_numeric_features; j++)
    {
      features.push_back(j);
    }
    shuffle(features.begin(), features.end(), env->generator);
  }
  vector<int> unused_features;
  vector<int> used_features;
  for (int j = 0; j < features.size(); j++)
  {
    if (j >= max_features)
    {
      unused_features.push_back(j);
    }
    else
    {
      used_features.push_back(j);
    }
  }
  bool split_found = rule.split_found;
  double bound = rule.imp;
  Result *res = new Result(split_found, bound, 0.0, rule.x1, rule.ind2, rule.ind1, total_numeric_features);
  vector<future<void>> fs;
  vector<bool> all_equal(total_numeric_features);
  int usable_features = 0;
  for (int j = 0; j < total_numeric_features; j++)
  {
    all_equal[j] = true;
    if (!regression)
    {
      double x0 = X[j][dataY[j][0][0]];
      for (int i = total_points - 1; i < total_points; i++)
      {
        double x = X[j][dataY[j][i][0]];
        if (x != x0)
        {
          all_equal[j] = false;
          usable_features++;
          break;
        }
      }
    }
    else
    {
      double x0 = X[j][(int)dataY_reg[j][0][0]];
      for (int i = total_points - 1; i < total_points; i++)
      {
        double x = X[j][(int)dataY_reg[j][i][0]];
        if (x != x0)
        {
          all_equal[j] = false;
          usable_features++;
          break;
        }
      }
    }
  }
  vector<bool> considered(total_numeric_features);
  for (int j = 0; j < total_numeric_features; j++)
  {
    considered[j] = false;
  }
  if (rule.single_considered.size() == total_numeric_features)
  {
    for (int j = 0; j < rule.single_considered.size(); j++)
    {
      if (rule.single_considered[j])
      {
        considered[j] = true;
      }
    }
  }
  int bisplits = 0;
  int bisplits2 = 0;
  int not_bisplits = 0;
  for (int k = 0; k < used_features.size(); k++)
  {
    int j1 = used_features[k];
    for (int l = k + 1; l < used_features.size(); l++)
    {
      bisplits2++;
      if (res->bestImp > LB)
      {
        int j2 = used_features[l];
        if ((all_equal[j1] && all_equal[j2]))
        {
          not_bisplits++;
          continue;
        }
        if (all_equal[j1] && considered[j2] == true)
        {
          not_bisplits++;
          continue;
        }
        if (all_equal[j2] && considered[j1] == true)
        {
          not_bisplits++;
          continue;
        }
        if (!all_equal[j1] && !all_equal[j2])
        {
          bisplits++;
        }
        considered[j1] = true;
        considered[j2] = true;
        if (!regression)
        {
          if (env->n_threads > 1)
          {
            fs.push_back(p->push(branchAndBound, X[j1], X[j2], j1, j2, n_classes, total_points, dataY[j1], dataY[j2], (double **)NULL, (double **)NULL, res, env));
            fs.push_back(p->push(branchAndBound, X[j2], X[j1], j2, j1, n_classes, total_points, dataY[j2], dataY[j1], (double **)NULL, (double **)NULL, res, env));
          }
          else
          {
            branchAndBound(0, X[j1], X[j2], j1, j2, n_classes, total_points, dataY[j1], dataY[j2], (double **)NULL, (double **)NULL, res, env);
            branchAndBound(0, X[j2], X[j1], j2, j1, n_classes, total_points, dataY[j2], dataY[j1], (double **)NULL, (double **)NULL, res, env);
          }
        }
        else
        {
          if (env->n_threads > 1)
          {
            fs.push_back(p->push(branchAndBound, X[j1], X[j2], j1, j2, n_classes, total_points, (int **)NULL, (int **)NULL, dataY_reg[j1], dataY_reg[j2], res, env));
            fs.push_back(p->push(branchAndBound, X[j2], X[j1], j2, j1, n_classes, total_points, (int **)NULL, (int **)NULL, dataY_reg[j2], dataY_reg[j1], res, env));
          }
          else
          {
            branchAndBound(0, X[j1], X[j2], j1, j2, n_classes, total_points, (int **)NULL, (int **)NULL, dataY_reg[j1], dataY_reg[j2], res, env);
            branchAndBound(0, X[j2], X[j1], j2, j1, n_classes, total_points, (int **)NULL, (int **)NULL, dataY_reg[j2], dataY_reg[j1], res, env);
          }
        }
      }
    }
  }
  for (int k = 0; k < fs.size(); k++)
  {
    fs[k].wait();
  }
  if (unused_features.size() > 0 && ((!regression && res->bestImp >= bound) || (regression && res->bestImp >= bound)))
  {
    for (int k = 0; k < unused_features.size(); k++)
    {
      int j1 = used_features[k];
      for (int l = 0; l < used_features.size(); l++)
      {
        int j2 = used_features[l];
        if ((all_equal[j1] && all_equal[j2]))
        {
          continue;
        }
        if (all_equal[j1] && considered[j2] == true)
        {
          continue;
        }
        if (all_equal[j2] && considered[j1] == true)
        {
          continue;
        }
        fs.clear();
        considered[j1] = true;
        considered[j2] = true;
        if (!regression)
        {
          fs.push_back(p->push(branchAndBound, X[j1], X[j2], j1, j2, n_classes, total_points, dataY[j1], dataY[j2], (double **)NULL, (double **)NULL, res, env));
          fs.push_back(p->push(branchAndBound, X[j2], X[j1], j2, j1, n_classes, total_points, dataY[j2], dataY[j1], (double **)NULL, (double **)NULL, res, env));
        }
        else
        {
          fs.push_back(p->push(branchAndBound, X[j1], X[j2], j1, j2, n_classes, total_points, (int **)NULL, (int **)NULL, dataY_reg[j1], dataY_reg[j2], res, env));
          fs.push_back(p->push(branchAndBound, X[j2], X[j1], j2, j1, n_classes, total_points, (int **)NULL, (int **)NULL, dataY_reg[j2], dataY_reg[j1], res, env));
        }
      }
      for (int k = 0; k < fs.size(); k++)
      {
        fs[k].wait();
      }
      if ((!regression && res->bestImp < impurity) || (regression && res->bestImp < total_points * mse))
      {
        break;
      }
      if (!all_equal[j1])
      {
        used_features.push_back(j1);
      }
    }
  }
  vector<double> a;
  double b = res->bestB;
  for (int j = 0; j < total_numeric_features; j++)
  {
    if (j == res->bestInd1)
    {
      if (j == res->bestInd2)
      {
        a.push_back(1);
        b = res->bestB / (1 - res->bestA);
      }
      else
      {
        a.push_back(-res->bestA);
      }
    }
    else if (j == res->bestInd2)
    {
      a.push_back(1);
    }
    else
    {
      a.push_back(0);
    }
  }
  bool to_leq = true;
  int j1 = res->bestInd1;
  int j2 = res->bestInd2;
  bool has_eq = false;
  double min_b = std::numeric_limits<double>::lowest();
  double max_b = std::numeric_limits<double>::max();
  bool divides = false;
  bool upper = false;
  bool lower = false;
  for (int i = 0; i < total_points; i++)
  {
    double tmp;
    if (!regression)
    {
      tmp = a[j1] * X[j1][dataY[0][i][0]] + a[j2] * X[j2][dataY[0][i][0]];
    }
    else
    {
      tmp = a[j1] * X[j1][(int)dataY_reg[0][i][0]] + a[j2] * X[j2][(int)dataY_reg[0][i][0]];
    }
    if (tmp == b)
    {
      to_leq = false;
      has_eq = true;
    }
    if (tmp >= b)
    {
      upper = true;
      if (tmp < max_b)
      {
        max_b = tmp;
      }
    }
    else
    {
      lower = true;
      if (tmp > min_b)
      {
        min_b = tmp;
      }
    }
  }
  if (upper && lower)
  {
    divides = true;
  }
  double b_orig = b;
  if (min_b > std::numeric_limits<double>::lowest() && max_b < std::numeric_limits<double>::max() && max_b > min_b)
  {
    to_leq = true;
    b = (min_b + max_b) / 2.;
  }
  if (divides)
  {
  }
  Rule rule2(a, b);
  rule.leq = to_leq;
  rule2.split_found = res->split_found;
  rule2.imp = res->bestImp;
  rule2.divides = res->divides;
  delete res;
  return rule2;
}
Rule Node::branchMultivariate(vector<int> features)
{
  bool split_found = false;
  double bound = std::numeric_limits<double>::max();
  Rule rule;
  bool remove_linear_dependent = true;
  if (remove_linear_dependent)
  {
    vector<int> features2 = getLinearIndependentFeatures(rule, 1e-8, features, true);
  }
  vector<int> order;
  order.push_back(features.size());
  for (int k = 0; k < features.size(); k++)
  {
    order.push_back(features[k]);
  }
  double **X2 = new double *[total_points];
  int *Y2;
  double *Y2_reg;
  if (!regression)
  {
    Y2 = new int[total_points];
  }
  else
  {
    Y2_reg = new double[total_points];
  }
  vector<int> used_features;
  vector<int> unused_features;
  vector<bool> all_equal(total_numeric_features);
  vector<bool> is_used(total_numeric_features);
  int usable_features = 0;
  for (int j2 = 0; j2 < features.size(); j2++)
  {
    int j = features[j2];
    all_equal[j] = true;
    if (!regression)
    {
      double x0 = X[j][dataY[j][0][0]];
      for (int i = total_points - 1; i < total_points; i++)
      {
        double x = X[j][dataY[j][i][0]];
        if (x != x0)
        {
          all_equal[j] = false;
          if (usable_features < features.size())
          {
            used_features.push_back(j);
            usable_features++;
            is_used[j] = true;
          }
          else
          {
            unused_features.push_back(j);
            is_used[j] = false;
          }
          break;
        }
        else
        {
          unused_features.push_back(j);
          is_used[j] = false;
        }
      }
    }
    else
    {
      double x0 = X[j][(int)dataY_reg[j][0][0]];
      for (int i = total_points - 1; i < total_points; i++)
      {
        double x = X[j][(int)dataY_reg[j][i][0]];
        if (x != x0)
        {
          all_equal[j] = false;
          if (usable_features < features.size())
          {
            used_features.push_back(j);
            usable_features++;
            is_used[j] = true;
          }
          else
          {
            unused_features.push_back(j);
            is_used[j] = false;
          }
          break;
        }
        else
        {
          unused_features.push_back(j);
          is_used[j] = false;
        }
      }
    }
  }
  if (usable_features == 0)
  {
    Rule rule;
    return rule;
  }
  if (usable_features == 1)
  {
    return branchUnivariate(features, features.size());
  }
  double imp_t = std::numeric_limits<double>::max();
  bool start_solution = false;
  double b_t = 0;
  vector<int> prefer;
  vector<double> a_t(used_features.size());
  for (int j = 0; j < used_features.size(); j++)
  {
    a_t[j] = 0;
  }
  Result *res = new Result(start_solution, imp_t, a_t, b_t);
  double **X_t = new double *[used_features.size()];
  int ***data_t = new int **[used_features.size()];
  double ***data_reg_t = new double **[used_features.size()];
  for (int j2 = 0; j2 < used_features.size(); j2++)
  {
    int j = used_features[j2];
    X_t[j2] = X[j];
    if (regression)
    {
      data_reg_t[j2] = dataY_reg[j];
    }
    else
    {
      data_t[j2] = dataY[j];
    }
  }
  for (int i = 0; i < total_points; i++)
  {
    X2[i] = new double[used_features.size() + 1];
    if (!regression)
    {
      Y2[i] = data_t[0][i][1];
      for (int j2 = 0; j2 < used_features.size(); j2++)
      {
        int j = used_features[j2];
        X2[i][j2] = X[j][data_t[0][i][0]];
      }
    }
    else
    {
      Y2_reg[i] = data_reg_t[0][i][1];
      for (int j2 = 0; j2 < used_features.size(); j2++)
      {
        int j = used_features[j2];
        X2[i][j2] = X[j][(int)data_reg_t[0][i][0]];
      }
    }
    X2[i][used_features.size()] = 1;
  }
  double *distr_S_G;
  if (!regression)
  {
    distr_S_G = new double[2 * n_classes];
    for (int i = 0; i < 2 * n_classes; i++)
    {
      distr_S_G[i] = 0;
    }
  }
  if (env->cross_entropy)
  {
    SimpleCrossEntropy(p, data_t, data_reg_t, X_t, order, X2, Y2, Y2_reg, total_points, used_features.size() + 1, features.size() + 1, n_classes, distr_S_G, res, env, 0, 0);
  }
  if (env->simplex)
  {
    BrightSide(p, dataY, dataY_reg, X, order, X2, Y2, Y2_reg, total_points, used_features.size() + 1, n_classes, distr_S_G, res, env, 0, 0);
  }
  if (!regression)
  {
    delete[] distr_S_G;
  }
  for (int i = 0; i < total_points; i++)
  {
    delete[] X2[i];
  }
  delete[] X2;
  if (!regression)
  {
    delete[] Y2;
  }
  else
  {
    delete[] Y2_reg;
  }
  vector<double> a(total_numeric_features);
  for (int j = 0; j < total_numeric_features; j++)
  {
    a[j] = 0;
  }
  for (int j2 = 0; j2 < used_features.size(); j2++)
  {
    int j = used_features[j2];
    a[j] = res->a[j2];
  }
  double b = res->bestB;
  Rule rule2(a, b);
  rule2.split_found = res->split_found;
  rule2.imp = res->bestImp;
  delete res;
  return rule2;
}
Rule Node::branchCross(Rule rule, int max_features)
{
  int *count = new int[total_numeric_features];
  int *is = new int[total_numeric_features];
  int **cX = new int *[total_numeric_features];
  for (int j = 0; j < total_numeric_features; j++)
  {
    cX[j] = new int[total_points + 1];
  }
  bool flag = false;
  for (int j = 0; j < total_numeric_features; j++)
  {
    count[j] = 2;
    is[j] = 1;
    cX[j][0] = 0;
    flag = true;
  }
  if (!env->regression)
  {
    for (int i = 0; i < total_points - 1; i++)
    {
      for (int j = 0; j < total_numeric_features; j++)
      {
        if (X[j][dataY[j][i][0]] < X[j][dataY[j][i + 1][0]])
        {
          cX[j][is[j]] = i;
          is[j] += 1;
          count[j] += 1;
        }
      }
    }
  }
  else
  {
    for (int i = 0; i < total_points - 1; i++)
    {
      for (int j = 0; j < total_numeric_features; j++)
      {
        if (X[j][(int)dataY_reg[j][i][0]] < X[j][(int)dataY_reg[j][i + 1][0]])
        {
          cX[j][is[j]] = i;
          is[j] += 1;
          count[j] += 1;
        }
      }
    }
  }
  for (int j = 0; j < total_numeric_features; j++)
  {
    cX[j][is[j]] = total_points - 1;
  }
  vector<int> features;
  if (rule.attribute_order.size() > 0 && false)
  {
    features = rule.attribute_order;
  }
  else
  {
    for (int j = 0; j < total_numeric_features; j++)
    {
      features.push_back(j);
    }
    shuffle(features.begin(), features.end(), env->generator);
  }
  vector<int> unused_features;
  vector<int> used_features;
  for (int j = 0; j < features.size(); j++)
  {
    if (j >= max_features)
    {
      unused_features.push_back(features[j]);
    }
    else
    {
      used_features.push_back(features[j]);
    }
  }
  bool split_found = rule.split_found;
  double bound = rule.imp;
  Result *res = new Result(split_found, bound, 0.0, rule.x1, rule.ind2, rule.ind1, total_numeric_features);
  vector<future<void>> fs;
  vector<bool> all_equal(total_numeric_features);
  int usable_features = 0;
  for (int j = 0; j < total_numeric_features; j++)
  {
    all_equal[j] = true;
    if (!regression)
    {
      double x0 = X[j][dataY[j][0][0]];
      for (int i = total_points - 1; i < total_points; i++)
      {
        double x = X[j][dataY[j][i][0]];
        if (x != x0)
        {
          all_equal[j] = false;
          usable_features++;
          break;
        }
      }
    }
    else
    {
      double x0 = X[j][(int)dataY_reg[j][0][0]];
      for (int i = total_points - 1; i < total_points; i++)
      {
        double x = X[j][(int)dataY_reg[j][i][0]];
        if (x != x0)
        {
          all_equal[j] = false;
          usable_features++;
          break;
        }
      }
    }
  }
  vector<bool> considered(total_numeric_features);
  for (int j = 0; j < total_numeric_features; j++)
  {
    considered[j] = false;
  }
  bool preprocess = true;
  if (preprocess)
  {
    bestInd1 = 0;
    bestInd2 = 1;
    bestX1 = 0;
    bestX2 = 0;
    bool split_found2 = false;
    vector<pair<int, double>> features2;
    double bestImp = std::numeric_limits<double>::max();
    double bestImp2 = std::numeric_limits<double>::max();
    int *N1;
    int *N2;
    int *N3;
    int *N4;
    int *N;
    if (!env->regression)
    {
      N1 = new int[n_classes];
      N2 = new int[n_classes];
      N3 = new int[n_classes];
      N4 = new int[n_classes];
      N = new int[n_classes];
      for (int m = 0; m < n_classes; m++)
      {
        N[m] = 0;
      }
      for (int i = 0; i < total_points; i++)
      {
        N[dataY[0][i][1]] += 1;
      }
    }
    else
    {
    }
    for (int j2 = 0; j2 < used_features.size(); j2++)
    {
      int j = used_features[j2];
      double val = std::numeric_limits<double>::max();
      double bestVal = std::numeric_limits<double>::max();
      for (int m = 0; m < n_classes; m++)
      {
        N1[m] = 0;
        N2[m] = N[m];
      }
      int right = total_points;
      for (int i = 0; i < total_points;)
      {
        double x;
        if (!env->regression)
        {
          x = X[j][dataY[j][i][0]];
        }
        else
        {
        }
        if (!env->regression)
        {
          while (i < total_points && X[j][dataY[j][i][0]] == x)
          {
            right--;
            N1[dataY[j][i][1]] += 1;
            N2[dataY[j][i][1]] -= 1;
            i++;
          }
          val = imp(N1, N2, n_classes, env->criterion);
        }
        else
        {
        }
        if (val <= bestImp && i < total_points)
        {
          bestImp2 = bestImp;
          bestInd2 = bestInd1;
          bestX2 = bestX1;
          split_found2 = true;
          bestImp = val;
          bestInd1 = j;
          bestX1 = x;
          if (i < total_points)
          {
            bestX1 = (X[j][dataY[j][i - 1][0]] + X[j][dataY[j][i][0]]) / 2;
          }
          else
          {
            bestX1 = (X[j][dataY[j][i - 1][0]]);
          }
        }
        else if (val <= bestImp2 && i < total_points)
        {
          bestInd2 = j;
          bestImp2 = val;
          bestX2 = x;
          if (i < total_points)
          {
            bestX2 = (X[j][dataY[j][i - 1][0]] + X[j][dataY[j][i][0]]) / 2;
          }
          else
          {
            bestX2 = (X[j][dataY[j][i - 1][0]]);
          }
        }
        if (val < bestVal)
        {
          bestVal = val;
        }
      }
      features2.push_back(std::pair<int, double>(j, bestVal + j * numeric_limits<double>::min()));
    }
    std::sort(features2.begin(), features2.end(), compareP);
    for (int m = 0; m < features2.size(); m++)
    {
      used_features[m] = features2[m].first;
    }
    if (split_found2)
    {
      Rule rule2(bestInd1, bestInd2, bestX1, bestX2);
      if (!env->regression)
      {
        for (int m = 0; m < n_classes; m++)
        {
          N1[m] = 0;
          N2[m] = 0;
          N3[m] = 0;
          N4[m] = 0;
        }
      }
      for (int i = 0; i < total_points; i++)
      {
        int ind = 0;
        if (!env->regression)
        {
          ind = rule2.evaluate(dataY[0][i][0], X, Xnominal);
        }
        else
        {
          ind = rule2.evaluate((int)dataY_reg[0][i][0], X, Xnominal);
        }
        if (ind == 0)
        {
          N1[dataY[0][i][1]]++;
        }
        else if (ind == 1)
        {
          N2[dataY[0][i][1]]++;
        }
        else if (ind == 2)
        {
          N3[dataY[0][i][1]]++;
        }
        else if (ind == 3)
        {
          N4[dataY[0][i][1]]++;
        }
      }
      double val = 0;
      if (!env->regression)
      {
        val = imp(N1, N2, n_classes, env->criterion) + imp(N3, N4, n_classes, env->criterion);
      }
      else
      {
      }
      res->split_found = true;
      res->bestImp = val;
      res->bestX1 = bestX1;
      res->bestX2 = bestX2;
      res->bestInd1 = bestInd1;
      res->bestInd2 = bestInd2;
      res->isCross = true;
    }
    delete[] N1;
    delete[] N2;
    delete[] N3;
    delete[] N4;
    delete[] N;
  }
  int bisplits = 0;
  int bisplits2 = 0;
  int not_bisplits = 0;
  if (res->bestImp != 0)
  {
    for (int k = 0; k < used_features.size(); k++)
    {
      int j1 = used_features[k];
      for (int l = k; l < used_features.size(); l++)
      {
        bisplits2++;
        if (true)
        {
          int j2 = used_features[l];
          considered[j1] = true;
          considered[j2] = true;
          if (!env->regression)
          {
            if (env->n_threads > 1)
            {
              fs.push_back(p->push(CrossBranchAndBound, X[j1], X[j2], j1, j2, n_classes, total_points, dataY[j1], dataY[j2], (double **)NULL, (double **)NULL, cX[j1], cX[j2], count[j1], count[j2], res, env));
            }
            else
            {
              CrossBranchAndBound(0, X[j1], X[j2], j1, j2, n_classes, total_points, dataY[j1], dataY[j2], (double **)NULL, (double **)NULL, cX[j1], cX[j2], count[j1], count[j2], res, env);
            }
          }
          else
          {
            if (env->n_threads > 1)
            {
              fs.push_back(p->push(CrossBranchAndBound, X[j1], X[j2], j1, j2, n_classes, total_points, (int **)NULL, (int **)NULL, dataY_reg[j1], dataY_reg[j2], cX[j1], cX[j2], count[j1], count[j2], res, env));
            }
            else
            {
              CrossBranchAndBound(0, X[j1], X[j2], j1, j2, n_classes, total_points, (int **)NULL, (int **)NULL, dataY_reg[j1], dataY_reg[j2], cX[j1], cX[j2], count[j1], count[j2], res, env);
            }
          }
        }
      }
    }
  }
  for (int k = 0; k < fs.size(); k++)
  {
    fs[k].wait();
  }
  delete[] count;
  delete[] is;
  for (int j = 0; j < total_numeric_features; j++)
  {
    delete[] cX[j];
  }
  delete[] cX;
  if (!res->isCross)
  {
    delete res;
    Rule rule2;
    return rule2;
  }
  Rule rule2((int)res->bestInd1, (int)res->bestInd2, res->bestX1, res->bestX2);
  rule2.imp = res->bestImp;
  rule2.split_found = res->split_found;
  delete res;
  return rule2;
}
vector<int> Node::getLinearIndependentFeatures(Rule &rule, double precision, vector<int> features, bool use_order)
{
  int count = 0;
  for (int j = 0; j < rule.single_considered.size(); j++)
  {
    if (rule.single_considered[j])
    {
      count++;
    }
  }
  vector<int> permutation(features.size() + 1);
  permutation[0] = 0;
  for (int j = 0; j < features.size(); j++)
  {
    permutation[j + 1] = features[j] + 1;
  }
  if (use_order && count == total_numeric_features && rule.attribute_order.size() == total_numeric_features)
  {
    permutation[0] = 0;
    for (int k = 0; k < rule.attribute_order.size(); k++)
    {
      permutation[k + 1] = rule.attribute_order[k] + 1;
    }
  }
  else
  {
    shuffle(permutation.begin() + 1, permutation.end(), env->generator);
  }
  Eigen::MatrixXd X_tmp(total_points, features.size() + 1);
  Eigen::VectorXd a_tmp(features.size() + 1);
  if (rule.split_found)
  {
    a_tmp(0) = -rule.b;
    for (int k = 1; k < permutation.size(); k++)
    {
      a_tmp(k) = rule.a[permutation[k] - 1];
    }
  }
  else
  {
    a_tmp(0) = 0;
    for (int k = 1; k < permutation.size(); k++)
    {
      a_tmp(k) = 0;
    }
  }
  for (int i = 0; i < total_points; i++)
  {
    X_tmp(i, 0) = 1;
  }
  for (int k = 1; k < permutation.size(); k++)
  {
    int j = permutation[k];
    for (int i = 0; i < total_points; i++)
    {
      if (!regression)
      {
        X_tmp(i, k) = X[j - 1][dataY[0][i][0]];
      }
      else
      {
        X_tmp(i, k) = X[j - 1][(int)dataY_reg[0][i][0]];
      }
    }
  }
  Eigen::VectorXd b = X_tmp * a_tmp;
  Eigen::VectorXd b2 = b;
  vector<int> features2;
  int n = X_tmp.rows();
  int d = X_tmp.cols();
  bool stop = false;
  vector<pair<int, int>> pivots;
  int k = 0;
  for (int i = 0; i < n; i++)
  {
    if (k >= d)
    {
      break;
    }
    int i_max = i;
    double max = abs(X_tmp(i, k));
    for (int i2 = i + 1; i2 < n; i2++)
    {
      if (max > precision)
      {
        break;
      }
      if (abs(X_tmp(i2, k)) > max)
      {
        i_max = i2;
        max = abs(X_tmp(i2, k));
      }
    }
    if (abs(X_tmp(i_max, k)) <= precision)
    {
      k = k + 1;
      i = i - 1;
    }
    else
    {
      if (permutation[k] >= 1)
      {
        features2.push_back(permutation[k] - 1);
      }
      Eigen::VectorXd vec_tmp = X_tmp.row(i);
      X_tmp.row(i) = X_tmp.row(i_max);
      X_tmp.row(i_max) = vec_tmp;
      double tmp = b(i);
      b(i) = b(i_max);
      b(i_max) = tmp;
      pivots.push_back(pair<int, int>(i, k));
      for (int i2 = i + 1; i2 < n; i2++)
      {
        double f = X_tmp(i2, k) / X_tmp(i, k);
        X_tmp(i2, k) = 0;
        for (int j = k + 1; j < d; j++)
        {
          X_tmp(i2, j) -= f * X_tmp(i, j);
        }
        b(i2) -= f * b(i);
      }
      k++;
    }
  }
  if (rule.split_found)
  {
    a_tmp = Eigen::VectorXd(pivots.size());
    for (int j = 0; j < total_numeric_features; j++)
    {
      rule.a[j] = 0;
    }
    rule.b = 0;
    for (int j = pivots.size() - 1; j >= 0; j--)
    {
      int i = pivots[j].first;
      int k = pivots[j].second;
      a_tmp(j) = b(i) / X_tmp(i, k);
      for (int i2 = i; i2 >= 0; i2--)
      {
        b(i2) -= a_tmp(j) * X_tmp(i2, k);
      }
      if (permutation[k] == 0)
      {
        rule.b = -a_tmp(j);
      }
      else
      {
        rule.a[permutation[k] - 1] = a_tmp(j);
      }
    }
  }
  return features2;
}
Rule Node::branchMultivariate(Rule rule, int max_features)
{
  vector<int> features;
  for (int j = 0; j < total_numeric_features; j++)
  {
    features.push_back(j);
  }
  bool remove_linear_dependent = true;
  vector<bool> is_linear_dependent(total_numeric_features, false);
  if (remove_linear_dependent)
  {
    for (int j = 0; j < total_numeric_features; j++)
    {
      is_linear_dependent[j] = true;
    }
    vector<int> features2 = getLinearIndependentFeatures(rule, 1e-8, features, true);
    sort(features2.begin(), features2.end());
    for (int k = 0; k < features2.size(); k++)
    {
      is_linear_dependent[features2[k]] = false;
    }
  }
  if (max_features < total_numeric_features)
  {
    shuffle(features.begin(), features.end(), env->generator);
  }
  bool split_found = rule.split_found;
  double bound = rule.imp;
  vector<int> order;
  order.push_back(features.size());
  for (int k = 0; k < features.size(); k++)
  {
    order.push_back(features[k]);
  }
  double **X2 = new double *[total_points];
  int *Y2;
  double *Y2_reg;
  if (!regression)
  {
    Y2 = new int[total_points];
  }
  else
  {
    Y2_reg = new double[total_points];
  }
  vector<int> used_features;
  vector<int> unused_features;
  vector<bool> all_equal(total_numeric_features);
  vector<bool> is_used(total_numeric_features);
  int usable_features = 0;
  for (int j = 0; j < total_numeric_features; j++)
  {
    all_equal[j] = true;
    if (!regression)
    {
      double x0 = X[j][dataY[j][0][0]];
      if (x0 < X[j][dataY[j][total_points - 1][0]])
      {
        all_equal[j] = false;
      }
      else
      {
        all_equal[j] = true;
      }
    }
    else
    {
      double x0 = X[j][(int)dataY_reg[j][0][0]];
      if (x0 < X[j][(int)dataY_reg[j][total_points - 1][0]])
      {
        all_equal[j] = false;
      }
      else
      {
        all_equal[j] = true;
      }
    }
  }
  int n_equal = 0;
  int n_linear_dependent = 0;
  for (int k = 0; k < features.size(); k++)
  {
    if (is_linear_dependent[features[k]])
    {
      n_linear_dependent++;
    }
    if (all_equal[features[k]])
    {
      n_equal++;
    }
  }
  for (int k = 0; k < features.size(); k++)
  {
    if (!is_linear_dependent[features[k]] && !all_equal[features[k]] && usable_features < max_features)
    {
      usable_features++;
      is_used[features[k]] = true;
      used_features.push_back(features[k]);
    }
    else
    {
      if (!all_equal[features[k]] && !is_linear_dependent[features[k]])
      {
        unused_features.push_back(features[k]);
      }
    }
  }
  if(used_features.size()<2){
    Rule rule3;
    rule3.split_found=false;
    return rule3;
  }
  double imp_t = std::numeric_limits<double>::max();
  bool start_solution = false;
  double b_t = 0;
  vector<int> prefer;
  if (rule.split_found)
  {
    for (int j = 0; j < total_numeric_features; j++)
    {
      is_used[j] = false;
    }
    for (int j = 0; j < total_numeric_features; j++)
    {
      if (abs(rule.a[j]) > 1e-14 && all_equal[j] == false)
      {
        prefer.push_back(j);
        is_used[j] = true;
      }
      else if (abs(rule.a[j]) > 1e-14)
      {
        if (!regression)
        {
          b_t -= rule.a[j] * X[j][dataY[j][0][0]];
        }
        else
        {
          b_t -= rule.a[j] * X[j][(int)dataY_reg[j][0][0]];
        }
      }
    }
    if (prefer.size() <= max_features)
    {
      for (int j2 = 0; j2 < used_features.size() && prefer.size() < max_features; j2++)
      {
        if (!is_used[used_features[j2]])
        {
          prefer.push_back(used_features[j2]);
        }
      }
      used_features = prefer;
      start_solution = true;
    }
  }
  sort(used_features.begin(), used_features.end());
  vector<double> a_t(used_features.size());
  for (int j = 0; j < used_features.size(); j++)
  {
    a_t[j] = 0;
  }
  if (start_solution)
  {
    for (int j = 0; j < used_features.size(); j++)
    {
      a_t[j] = rule.a[used_features[j]];
    }
    b_t += rule.b;
    imp_t = rule.imp;
  }
  Result *res = new Result(start_solution, imp_t, a_t, b_t);
  double **X_t = new double *[used_features.size()];
  int ***data_t = new int **[used_features.size()];
  double ***data_reg_t = new double **[used_features.size()];
  for (int j2 = 0; j2 < used_features.size(); j2++)
  {
    int j = used_features[j2];
    X_t[j2] = X[j];
    if (regression)
    {
      data_reg_t[j2] = dataY_reg[j];
    }
    else
    {
      data_t[j2] = dataY[j];
    }
  }
  for (int i = 0; i < total_points; i++)
  {
    X2[i] = new double[used_features.size() + 1];
    if (!regression)
    {
      Y2[i] = data_t[0][i][1];
      for (int j2 = 0; j2 < used_features.size(); j2++)
      {
        int j = used_features[j2];
        X2[i][j2] = X[j][data_t[0][i][0]];
      }
    }
    else
    {
      Y2_reg[i] = data_reg_t[0][i][1];
      for (int j2 = 0; j2 < used_features.size(); j2++)
      {
        int j = used_features[j2];
        X2[i][j2] = X[j][(int)data_reg_t[0][i][0]];
      }
    }
    X2[i][used_features.size()] = 1;
  }
  double *distr_S_G;
  if (!regression)
  {
    distr_S_G = new double[2 * n_classes];
    for (int i = 0; i < 2 * n_classes; i++)
    {
      distr_S_G[i] = 0;
    }
  }
  double impX = 0;
  double impSA = 0;
  if (env->cross_entropy || env->simplex)
  {
    if (env->cross_entropy)
    {
      SimpleCrossEntropy(p, data_t, data_reg_t, X_t, order, X2, Y2, Y2_reg, total_points, used_features.size() + 1, features.size() + 1, n_classes, distr_S_G, res, env, 0, 0);
    }
    if (env->simplex)
    {
      double imp_start = res->bestImp;
      BrightSide(p, dataY, dataY_reg, X, order, X2, Y2, Y2_reg, total_points, used_features.size() + 1, n_classes, distr_S_G, res, env, 0, 0);
    }
  }
  delete[] X_t;
  delete[] data_t;
  delete[] data_reg_t;
  for (int i = 0; i < total_points; i++)
  {
    delete[] X2[i];
  }
  delete[] X2;
  if (!regression)
  {
    delete[] Y2;
  }
  else
  {
    delete[] Y2_reg;
  }
  vector<double> a(total_numeric_features);
  for (int j = 0; j < total_numeric_features; j++)
  {
    a[j] = 0;
  }
  for (int j2 = 0; j2 < used_features.size(); j2++)
  {
    int j = used_features[j2];
    a[j] = res->a[j2];
  }
  double b = res->bestB;
  Rule rule2(a, b);
  rule2.split_found = res->split_found;
  rule2.imp = res->bestImp;
  if (!regression)
  {
    delete[] distr_S_G;
  }
  delete res;
  return rule2;
}
bool Node::evaluate(vector<double> &x, vector<double> &x_nominal)
{
  return cond.evaluate(x, x_nominal);
}
int Node::classify()
{
  return maxVote;
}
Node *Node::copy(int id, double **X, double **Xnominal, int *Y, double *Y_reg, Environment *env, bool take_foreign, bool data)
{
  Node *node = new Node();
  node->regression = regression;
  node->total_numeric_features = total_numeric_features;
  node->total_nominal_features = total_nominal_features;
  node->env = env;
  node->total_features = total_features;
  node->total_points = total_points;
  node->n_classes = n_classes;
  node->isLeaf = isLeaf;
  node->hasData = false;
  node->possesses_data = false;
  node->id = id;
  node->str_id = "" + str_id;
  node->maxVote = maxVote;
  node->maxVote_str = "" + maxVote_str;
  node->impurity = impurity;
  node->shouldSplit = shouldSplit;
  node->shouldSplit_numeric = shouldSplit_numeric;
  node->shouldSplit_nominal = shouldSplit_nominal;
  node->random_state = random_state;
  node->mean = mean;
  node->mse = mse;
  node->mae = mae;
  node->sq_error = sq_error;
  node->a_error = a_error;
  node->classification_error = classification_error;
  node->errorC = errorC;
  node->LB = LB;
  node->bestInd1 = bestInd1;
  node->bestInd2 = bestInd2;
  node->bestA = bestA;
  node->bestB = bestB;
  node->bestX1 = bestX1;
  node->bestX2 = bestX2;
  node->type = type;
  node->split = split;
  node->errorC = errorC;
  node->depth = depth;
  node->impS1 = impS1;
  node->impS2 = impS2;
  node->impS = impS;
  node->isInTree = isInTree;
  node->p = p;
  for (int i = 0; i < nominal_domain.size(); i++)
  {
    node->nominal_domain.push_back(nominal_domain[i]);
  }
  for (int i = 0; i < is_nominal_val_in_dataset.size(); i++)
  {
    node->is_nominal_val_in_dataset.push_back(vector<bool>());
    for (int j = 0; j < is_nominal_val_in_dataset[i].size(); j++)
    {
      node->is_nominal_val_in_dataset[i].push_back(is_nominal_val_in_dataset[i][j]);
    }
  }
  for (int i = 0; i < a.size(); i++)
  {
    node->a.push_back(a[i]);
  }
  node->cond = cond.copy();
  if (!regression)
  {
    node->Ni = new double[n_classes];
    for (int i = 0; i < n_classes; i++)
    {
      node->Ni[i] = Ni[i];
    }
  }
  if (data)
  {
    if (regression)
    {
      node->Y_reg = new double[total_points];
      for (int i = 0; i < total_points; i++)
      {
        node->Y_reg[i] = this->Y_reg[i];
      }
    }
    else
    {
      node->Y = new int[total_points];
      for (int i = 0; i < total_points; i++)
      {
        node->Y[i] = this->Y[i];
      }
    }
    if (total_numeric_features > 0)
    {
      node->X = new double *[total_numeric_features];
      for (int j = 0; j < total_numeric_features; j++)
      {
        node->X[j] = new double[total_points];
      }
      for (int j = 0; j < total_numeric_features; j++)
      {
        for (int i = 0; i < total_points; i++)
        {
          node->X[j][i] = this->X[j][i];
        }
      }
    }
    if (total_nominal_features > 0)
    {
      node->Xnominal = new double *[total_points];
      for (int i = 0; i < total_points; i++)
      {
        node->Xnominal[i] = new double[total_nominal_features];
        for (int j = 0; j < total_nominal_features; j++)
        {
          node->Xnominal[i][j] = this->Xnominal[i][j];
        }
      }
    }
    if (!regression)
    {
      node->dataY = new int **[max(1, total_numeric_features)];
      for (int j = 0; j < max(1, total_numeric_features); j++)
      {
        node->dataY[j] = new int *[total_points];
        for (int i = 0; i < total_points; i++)
        {
          node->dataY[j][i] = new int[2];
          node->dataY[j][i][0] = dataY[j][i][0];
          node->dataY[j][i][1] = dataY[j][i][1];
        }
      }
    }
    else
    {
      node->dataY_reg = new double **[max(1, total_numeric_features)];
      for (int j = 0; j < max(1, total_numeric_features); j++)
      {
        node->dataY_reg[j] = new double *[total_points];
        for (int i = 0; i < total_points; i++)
        {
          node->dataY_reg[j][i] = new double[2];
          node->dataY_reg[j][i][0] = dataY_reg[j][i][0];
          node->dataY_reg[j][i][1] = dataY_reg[j][i][1];
        }
      }
    }
    node->hasData = true;
    node->possesses_data = true;
  }
  else
  {
    if (take_foreign)
    {
      node->X = this->X;
      node->Xnominal = this->Xnominal;
      node->Y = this->Y;
      node->Y_reg = this->Y_reg;
      node->dataY = this->dataY;
      node->dataY_reg = this->dataY_reg;
      node->possesses_data = false;
      node->hasData = true;
    }
    else
    {
      node->X = X;
      node->Xnominal = Xnominal;
      node->Y = Y;
      node->Y_reg = Y_reg;
    }
  }
  return node;
}
#define Malloc(type, n) (type *)malloc((n) * sizeof(type))
void print_null(const char *s)
{
}
Rule Node::beautify(Rule &rule, bool max_margin)
{
  if (rule.is_cross || rule.is_nominal || !max_margin)
  {
    return rule;
  }
  bool reduce = true;
  int count = 0;
  for (int j = 0; j < rule.a.size(); j++)
  {
    if (abs(rule.a[j]) > 0)
    {
      count++;
    }
  }
  if (count == 0)
  {
    reduce = false;
  }
  if (count == 1)
  {
    reduce = false;
  }
  if (reduce)
  {
    bool reduced = false;
    count = 0;
    for (int j = 0; j < rule.single_considered.size(); j++)
    {
      if (rule.single_considered[j])
      {
        count++;
      }
    }
    vector<int> order;
    for (int j = 0; j < total_numeric_features + 1; j++)
    {
      order.push_back(j);
    }
    if (count == total_numeric_features && rule.attribute_order.size() == total_numeric_features)
    {
      for (int j = total_numeric_features - 1; j >= 0; j--)
      {
        order[j] = rule.attribute_order[j];
      }
    }
    count = 0;
    for (int k = 0; k < total_numeric_features; k++)
    {
      int j = order[k];
      if (rule.a[j] == 0)
      {
        continue;
      }
      double min = std::numeric_limits<double>::lowest();
      double max = std::numeric_limits<double>::max();
      double sum;
      for (int i = 0; i < total_points; i++)
      {
        if (!env->regression && X[j][dataY[0][i][0]] == 0)
        {
          continue;
        }
        if (env->regression && X[j][(int)dataY_reg[0][i][0]] == 0)
        {
          continue;
        }
        sum = -rule.b;
        for (int j2 = 0; j2 < total_numeric_features; j2++)
        {
          if (!env->regression)
          {
            sum += rule.a[j2] * X[j2][dataY[0][i][0]];
          }
          else
          {
            sum += rule.a[j2] * X[j2][(int)dataY_reg[0][i][0]];
          }
        }
        bool case1 = true;
        if (sum < 0)
        {
          case1 = false;
        }
        double tmp;
        double x;
        if (!env->regression)
        {
          sum -= rule.a[j] * X[j][dataY[0][i][0]];
          tmp = -sum / X[j][dataY[0][i][0]];
          x = X[j][dataY[0][i][0]];
        }
        else
        {
          sum -= rule.a[j] * X[j][(int)dataY_reg[0][i][0]];
          tmp = -sum / X[j][(int)dataY_reg[0][i][0]];
          x = X[j][(int)dataY_reg[0][i][0]];
        }
        if (case1)
        {
          if (x > 0)
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
          if (x < 0)
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
        rule.a[j] = 0;
        count++;
      }
      else
      {
        if (min > std::numeric_limits<double>::lowest() && max < std::numeric_limits<double>::max() && (max - min) > 1e-6)
        {
          rule.a[j] = min + (max - min) / 2.;
        }
      }
    }
  }
  double min = std::numeric_limits<double>::lowest();
  double max = std::numeric_limits<double>::max();
  double sum;
  vector<double> pm(total_points);
  vector<int> used_features;
  for (int j = 0; j < total_numeric_features; j++)
  {
    if (rule.a[j] != 0)
    {
      used_features.push_back(j);
    }
  }
  double mindiff = std::numeric_limits<double>::max();
  if (!regression)
  {
    for (int i = 0; i < total_points; i++)
    {
      sum = -rule.b;
      for (int j = 0; j < total_numeric_features; j++)
      {
        sum += rule.a[j] * X[j][dataY[0][i][0]];
      }
      if (abs(sum) < mindiff)
      {
        mindiff = abs(sum);
      }
      if (sum < 0)
      {
        pm[i] = -1;
        if (sum > min)
        {
          min = sum;
        }
      }
      if (sum >= 0)
      {
        pm[i] = 1;
        if (sum < max)
        {
          max = sum;
        }
      }
    }
  }
  else
  {
    for (int i = 0; i < total_points; i++)
    {
      sum = -rule.b;
      for (int j = 0; j < total_numeric_features; j++)
      {
        sum += rule.a[j] * X[j][(int)dataY_reg[0][i][0]];
      }
      if (abs(sum) < mindiff)
      {
        mindiff = abs(sum);
      }
      if (sum < 0)
      {
        pm[i] = -1;
        if (sum > min)
        {
          min = sum;
        }
      }
      if (sum >= 0)
      {
        pm[i] = 1;
        if (sum < max)
        {
          max = sum;
        }
      }
    }
  }
  rule.b += (max + min) / 2;
  rule.leq = true;
  if (used_features.size() == 1 || !max_margin)
  {
    return rule;
  }
  bool liblinear = false;
  bool clp = false;
  bool gurobi = true;
  #ifndef use_gurobi
    gurobi=false;
  #endif
  if(env->margin_norm==2 && gurobi==false){
    liblinear=true;
    clp=false;
  }
  else{
    clp=true;
  }
  if (liblinear)
  {
    struct parameter param;
    struct problem prob;
    struct model *model;
    struct feature_node *x_space;
    int elements = (used_features.size() + 1) * total_points;
    prob.l = total_points;
    prob.bias = 1;
    prob.y = Malloc(double, prob.l);
    prob.x = Malloc(struct feature_node *, prob.l);
    x_space = Malloc(struct feature_node, elements + prob.l);
    int j2 = 0;
    int max_index = 0;
    for (int i = 0; i < total_points; i++)
    {
      prob.x[i] = &x_space[j2];
      prob.y[i] = pm[i];
      for (int j = 0; j < used_features.size(); j++)
      {
        if (!regression)
        {
          x_space[j2].index = j + 1;
          x_space[j2].value = X[used_features[j]][dataY[0][i][0]];
        }
        else
        {
          x_space[j2].index = j + 1;
          x_space[j2].value = X[used_features[j]][(int)dataY_reg[0][i][0]];
        }
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
    param.solver_type = L2R_L2LOSS_SVC;
    param.C = 1e10;
    param.p = 0.1;
    param.nu = 0.;
    param.eps = 1e-3;
    param.nr_weight = 0;
    param.regularize_bias = 0;
    param.weight_label = NULL;
    param.weight = NULL;
    param.init_sol = NULL;
    double tmp = 0;
    param.init_sol = Malloc(double, used_features.size() + 1);
    if (mindiff > 0)
    {
      mindiff = 0.9 * mindiff;
    }
    else
    {
      mindiff = 1;
    }
    for (int j = 0; j < used_features.size(); j++)
    {
      param.init_sol[j] = rule.a[used_features[j]] / mindiff;
      tmp += (rule.a[used_features[j]] * rule.a[used_features[j]]) / 2.;
    }
    tmp += (rule.b * rule.b) / 2.;
    param.init_sol[used_features.size()] = -rule.b / mindiff;
    set_print_string_function(print_null);
    model = train(&prob, &param);
    vector<double> a_t(total_numeric_features, 0);
    double b_t = 0;
    for (int j = 0; j < used_features.size(); j++)
    {
      a_t[used_features[j]] = model->w[j];
    }
    b_t = -model->w[used_features.size()];
    int errors = 0;
    for (int i = 0; i < total_points; i++)
    {
      double sum1 = -rule.b;
      double sum2 = -b_t;
      int sig1 = 1;
      int sig2 = 1;
      for (int j = 0; j < total_numeric_features; j++)
      {
        if (!regression)
        {
          sum1 += rule.a[j] * X[j][dataY[0][i][0]];
          sum2 += a_t[j] * X[j][dataY[0][i][0]];
        }
        else
        {
          sum1 += rule.a[j] * X[j][(int)dataY_reg[0][i][0]];
          sum2 += a_t[j] * X[j][(int)dataY_reg[0][i][0]];
        }
      }
      if (sum1 <= 0)
      {
        sig1 = -1;
      }
      if (sum2 <= 0)
      {
        sig2 = -1;
      }
      if (sum1 == 0)
      {
        sig1 = 0;
      }
      if (sum2 == 0)
      {
        sig2 = 0;
      }
      if (sig1 != sig2)
      {
        errors++;
      }
    }
    if (errors == 0 || errors == total_points)
    {
      rule.a = a_t;
      rule.b = b_t;
    }
    free_and_destroy_model(&model);
    destroy_param(&param);
    free(prob.y);
    free(prob.x);
    free(x_space);
  }
  #ifdef use_gurobi
  if (gurobi)
  {
    grb_env.start();
    GRBModel model = GRBModel(grb_env);
    model.set(GRB_IntParam_NumericFocus, 3);
    model.set(GRB_IntParam_OutputFlag, 0);
    model.set(GRB_IntParam_LogToConsole, 0);
    int d=total_numeric_features+1;
    vector<GRBVar> x(d);
    for (int j = 0; j < d; j++)
    {
        if (rule.a[j] == 0 && j < d - 1)
        {
            x[j] = model.addVar(-0, 0, 0.0, GRB_CONTINUOUS, "x" + j);
        }
        else
        {
            x[j] = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "x" + j);
        }
    }
    if(env->margin_norm==1){
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
          obj += abs[j];
      }
      model.setObjective(obj, GRB_MINIMIZE);
    }
    else{
      GRBQuadExpr obj = GRBQuadExpr();
      for (int j = 0; j < d - 1; j++)
      {
          obj += x[j] * x[j];
      }
    }
    for (int i = 0; i < total_points; i++)
    {
        double sum = -rule.b;
        GRBLinExpr e = GRBLinExpr();
        e+=x[d-1];
        if(!regression){
            for (int j = 0; j < d-1; j++)
            {
                sum += X[j][dataY[0][i][0]] * rule.a[j];
                e+=x[j]*X[j][dataY[0][i][0]];
            }
        }
        else{
          for (int j = 0; j < d-1; j++)
            {
                sum += X[j][(int) dataY_reg[0][i][0]] * rule.a[j];
                e+=X[j][(int) dataY_reg[0][i][0]]*x[j];
            }
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
    int status = model.get(GRB_IntAttr_Status);
    if (status == GRB_OPTIMAL || status == GRB_SUBOPTIMAL)
    {
      double b_t = -x[d-1].get(GRB_DoubleAttr_X);
      vector<double> a_t(total_numeric_features, 0);
      for(int j=0;j<d-1;j++){
        a_t[j]=x[j].get(GRB_DoubleAttr_X);
      }
      int errors = 0;
      for (int i = 0; i < total_points; i++)
      {
        double sum1 = -rule.b;
        double sum2 = -b_t;
        int sig1 = 1;
        int sig2 = 1;
        for (int j = 0; j < total_numeric_features; j++)
        {
          if (!regression)
          {
            sum1 += rule.a[j] * X[j][dataY[0][i][0]];
            sum2 += a_t[j] * X[j][dataY[0][i][0]];
          }
          else
          {
            sum1 += rule.a[j] * X[j][(int)dataY_reg[0][i][0]];
            sum2 += a_t[j] * X[j][(int)dataY_reg[0][i][0]];
          }
        }
        if (sum1 <= 0)
        {
          sig1 = -1;
        }
        if (sum2 <= 0)
        {
          sig2 = -1;
        }
        if (sum1 == 0)
        {
          sig1 = 0;
        }
        if (sum2 == 0)
        {
          sig2 = 0;
        }
        if (sig1 != sig2)
        {
          errors++;
        }
      }
      if (errors == 0 || errors == total_points)
      {
          rule.a = a_t;
          rule.b = b_t;
      }
    }
  }
  #endif
  #ifndef use_gurobi
  if(clp){
    int d=total_numeric_features+1;
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
        if (rule.a[j] == 0 && j < d - 1)
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
        clpmodel.setObjectiveCoefficient(abs_ind[j],1);
    }
    CoinBuild buildObject;
    int** rowInds=new int*[total_points];
    double** rowVals=new double*[total_points];
    for(int i=0;i<total_points;i++){
        rowInds[i]=new int[d];
        rowVals[i]=new double[d];
        rowInds[i][d-1]=a_ind[d-1];
        rowVals[i][d-1]=1;
        double sum = -rule.b;
        if(!regression){
          for (int j = 0; j < d-1; j++)
          {
              sum += X[j][dataY[0][i][0]] * rule.a[j];
              rowInds[i][j]=a_ind[j];
              rowVals[i][j]=X[j][dataY[0][i][0]];
          }
        }
        else{
          for (int j = 0; j < d-1; j++)
          {
              sum += X[j][(int) dataY_reg[0][i][0]] * rule.a[j];
              rowInds[i][j]=a_ind[j];
              rowVals[i][j]=X[j][(int) dataY_reg[0][i][0]];
          }
        }
        if(sum>=0){
            buildObject.addRow(d, rowInds[i], rowVals[i],1.0, COIN_DBL_MAX);
        }
        else{
            buildObject.addRow(d, rowInds[i], rowVals[i],-COIN_DBL_MAX, -1.0);
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
    if(has_solution){
        double * columnPrimal = clpmodel.primalColumnSolution();
        double b_t = -columnPrimal[d-1];
        vector<double> a_t(total_numeric_features, 0);
        for(int j=0;j<d-1;j++){
          a_t[j]= columnPrimal[a_ind[j]];
        }
        int errors = 0;
        for (int i = 0; i < total_points; i++)
        {
          double sum1 = -rule.b;
          double sum2 = -b_t;
          int sig1 = 1;
          int sig2 = 1;
          for (int j = 0; j < total_numeric_features; j++)
          {
            if (!regression)
            {
              sum1 += rule.a[j] * X[j][dataY[0][i][0]];
              sum2 += a_t[j] * X[j][dataY[0][i][0]];
            }
            else
            {
              sum1 += rule.a[j] * X[j][(int)dataY_reg[0][i][0]];
              sum2 += a_t[j] * X[j][(int)dataY_reg[0][i][0]];
            }
          }
          if (sum1 <= 0)
          {
            sig1 = -1;
          }
          if (sum2 <= 0)
          {
            sig2 = -1;
          }
          if (sum1 == 0)
          {
            sig1 = 0;
          }
          if (sum2 == 0)
          {
            sig2 = 0;
          }
          if (sig1 != sig2)
          {
            errors++;
          }
        }
        if (errors == 0 || errors == total_points)
        {
          rule.a = a_t;
          rule.b = b_t;
        }
    }
  }
  #endif
  return rule;
}
Rule Node::chooseBetter(Rule &rule)
{
  if (!rule.is_cross)
  {
    return rule;
  }
  vector<double> imps(2, 0);
  vector<double> splitfound(2, false);
  if (!env->regression)
  {
    for (int k = 0; k < 2; k++)
    {
      int ind = rule.ind1;
      double x = rule.x1;
      if (k == 1)
      {
        ind = rule.ind2;
        x = rule.x2;
      }
      vector<double> N(2 * n_classes, 0);
      bool left = false;
      bool right = false;
      bool is_equ = false;
      for (int i = 0; i < total_points; i++)
      {
        if (X[ind][dataY[0][i][0]] == x)
        {
          is_equ = true;
        }
        if (X[ind][dataY[0][i][0]] <= x)
        {
          N[dataY[0][i][1]] += 1;
          left = true;
        }
        else
        {
          N[n_classes + dataY[0][i][1]] += 1;
          right = true;
        }
      }
      if (left && right)
      {
        splitfound[k] = true;
        imps[k] = imp(N, n_classes, env->criterion);
      }
      else
      {
        imps[k] = std::numeric_limits<double>::max();
      }
    }
    if (!splitfound[0] && !splitfound[1])
    {
      rule.split_found = false;
      return rule;
    }
    int ind = rule.ind1;
    double x = rule.x1;
    int ind2 = 0;
    if (imps[1] <= imps[0] && splitfound[1])
    {
      ind = rule.ind2;
      x = rule.x2;
      ind2 = 1;
    }
    vector<double> a(total_numeric_features, 0);
    a[ind] = 1;
    Rule rule2(a, x);
    rule2.leq = true;
    rule2.imp = imps[ind2];
    rule2.split_found = true;
    return rule2;
  }
  else
  {
    for (int k = 0; k < 2; k++)
    {
      int ind = rule.ind1;
      double x = rule.x1;
      if (k == 1)
      {
        ind = rule.ind2;
        x = rule.x2;
      }
      vector<Vector_reg> N_reg(2);
      if (env->criterion == Environment::mae)
      {
        N_reg[0].use_mae();
        N_reg[1].use_mae();
      }
      bool left = false;
      bool right = false;
      for (int i = 0; i < total_points; i++)
      {
        if (X[ind][(int)dataY_reg[0][i][0]] <= x)
        {
          N_reg[0].push_back(dataY_reg[0][i][1]);
          left = true;
        }
        else
        {
          N_reg[1].push_back(dataY_reg[0][i][1]);
          right = true;
        }
      }
      if (left && right)
      {
        splitfound[k] = true;
        imps[k] = N_reg[0].error(env->criterion) + N_reg[1].error(env->criterion);
      }
      else
      {
        imps[k] = std::numeric_limits<double>::max();
      }
    }
    if (!splitfound[0] && !splitfound[1])
    {
      rule.split_found = false;
      return rule;
    }
    int ind = rule.ind1;
    double x = rule.x1;
    int ind2 = 0;
    if (imps[1] <= imps[0] && splitfound[1])
    {
      ind = rule.ind2;
      x = rule.x2;
      ind2 = 1;
    }
    vector<double> a(total_numeric_features, 0);
    a[ind] = 1;
    Rule rule2(a, x);
    rule2.leq = true;
    rule2.imp = imps[ind2];
    rule2.split_found = true;
    return rule2;
  }
}
double Node::intrinsicInformation(Rule &rule)
{
  vector<int> total_points_tmp(rule.no_children, 0);
  for (int i = 0; i < total_points; i++)
  {
    int ind = 0;
    if (!regression)
    {
      ind = rule.evaluate(dataY[0][i][0], X, Xnominal);
    }
    else
    {
      ind = rule.evaluate((int)dataY_reg[0][i][0], X, Xnominal);
    }
    total_points_tmp[ind] += 1;
  }
  double iI = 0;
  double tmp;
  for (int l = 0; l < total_points_tmp.size(); l++)
  {
    if (total_points_tmp[l] != 0)
    {
      tmp = ((double)total_points_tmp[l] / (double)total_points);
      iI -= tmp * log2(tmp);
    }
  }
  return iI;
}
