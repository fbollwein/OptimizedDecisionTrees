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
#include <iostream>
#include <sys/time.h>
#include <string>
#include <regex>
#include <fstream>
#include <limits>
#include <math.h>
#include "General.h"
#include "DecisionTree.h"
#include "vMF.h"
#include <filesystem>
using namespace std;
bool is_number(const std::string &s)
{
	std::string::const_iterator it = s.begin();
	while (it != s.end() && std::isdigit(*it))
	{
		it++;
	}
	return !s.empty() && it == s.end();
}
void DecisionTree::setParam(vector<string> params, vector<string> values)
{
	for (int i = 0; i < params.size(); i++)
	{
		setParam(params[i], values[i]);
	}
}
void DecisionTree::setParam(vector<string> params, vector<double> values)
{
	for (int i = 0; i < params.size(); i++)
	{
		setParam(params[i], values[i]);
	}
}
void DecisionTree::setParam(string param, string value)
{
	if(is_number(value)){
		setParam(param,atof(value.c_str()));
	}
	if (param.compare("criterion") == 0)
	{
		if (value.compare("gini") == 0)
		{
			env.criterion = env.gini_imp;
		}
		if (value.compare("entropy") == 0)
		{
			env.criterion = env.entropy_imp;
		}
		if (value.compare("twoing") == 0)
		{
			env.criterion = env.twoing;
		}
		if (value.compare("mse") == 0)
		{
			env.criterion_reg = env.mse;
		}
		if (value.compare("mae") == 0)
		{
			env.criterion_reg = env.mae;
		}
	}
	else if (param.compare("cross_entropy_samples") == 0 || param.compare("ce_samples") == 0)
	{
		env.cross_entropy_samples = -1;
		env.cross_entropy_samples_form = string(value.c_str());
		env.cross_entropy_samples_expression = FunctionParser();
		env.cross_entropy_samples_expression.parse(string(value.c_str()));
	}
	else if (param.compare("cross_entropy_no_improvement") == 0 || param.compare("ce_no_improvement") == 0)
	{
		env.cross_entropy_no_improvement = -1;
		env.cross_entropy_no_improvement_form = string(value.c_str());
		env.cross_entropy_no_improvement_expression = FunctionParser();
		env.cross_entropy_no_improvement_expression.parse(string(value.c_str()));
	}
	else if (param.compare("sa_start_iterations") == 0)
	{
		env.spx_start_iterations = -1;
		env.spx_start_iterations_form = string(value.c_str());
		env.spx_start_iterations_expression = FunctionParser();
		env.spx_start_iterations_expression.parse(string(value.c_str()));
	}
	else if (param.compare("sa_no_improvement") == 0)
	{
		env.spx_no_improvement = -1;
		env.spx_no_improvement_form = string(value.c_str());
		env.spx_no_improvement_expression = FunctionParser();
		env.spx_no_improvement_expression.parse(string(value.c_str()));
	}
	else if (param.compare("sa_scale") == 0)
	{
		if (value.compare("no") == 0)
		{
			env.spx_scale = -1;
		}
		else if (value.compare("geom") == 0)
		{
			env.spx_scale = 0;
		}
		else if (value.compare("equ") == 0)
		{
			env.spx_scale = 1;
		}
	}
	else
	{
		cout << "String parameter \"" << param << "\" unknown!\n";
	}
}
void DecisionTree::setParam(string param, double value)
{
	if (param.compare("criterion") == 0)
	{
		if (value==Environment::gini_imp)
		{
			env.criterion = env.gini_imp;
		}
		if (value==Environment::entropy_imp)
		{
			env.criterion = env.entropy_imp;
		}
		if (value==Environment::twoing)
		{
			env.criterion = env.twoing;
		}
		if (value==Environment::mse)
		{
			env.criterion_reg = env.mse;
		}
		if (value==Environment::mae)
		{
			env.criterion_reg = env.mae;
		}
	}
	else if (param.compare("max_depth") == 0)
	{
		env.max_depth = value;
	}
	else if (param.compare("min_samples_split") == 0)
	{
		env.min_samples_split = value;
	}
	else if (param.compare("min_cases") == 0)
	{
		env.min_cases = value;
	}
	else if (param.compare("min_weight_fraction_leaf") == 0)
	{
		env.min_weight_fraction_leaf = value;
	}
	else if (param.compare("max_features") == 0)
	{
		env.max_features = value;
	}
	else if (param.compare("max_leaf_nodes") == 0)
	{
		env.max_leaf_nodes = value;
	}
	else if (param.compare("min_impurity_decrease") == 0)
	{
		env.min_impurity_decrease = value;
	}
	else if (param.compare("min_impurity_split") == 0)
	{
		env.min_impurity_split = value;
	}
	else if (param.compare("maximize_margin") == 0 || param.compare("maximise_margin") == 0 || param.compare("max_margin") == 0)
	{
		if (value==1)
		{
			env.max_margin = true;
			env.margin_norm = 1;
		}
		else if (value == 2)
		{
			env.max_margin = true;
			env.margin_norm = 2;
		}
		else
		{
			env.max_margin = false;
		}
	}
	else if (param.compare("univariate") == 0 || param.compare("1d") == 0)
	{
		if (value==0)
		{
			if (!env.simplex)
			{
				env.univariate = false;
			}
		}
		else
		{
			env.univariate = true;
		}
	}
	else if (param.compare("multi") == 0 || param.compare("nd") == 0)
	{
		if (value==0)
		{
			env.mult = false;
		}
		else
		{
			env.mult = true;
			env.cross_entropy = true;
			env.cross = false;
		}
	}
	else if (param.compare("bi") == 0 || param.compare("2d") == 0)
	{
		if (value==0)
		{
			env.bi = false;
		}
		else
		{
			env.bi = true;
			env.cross = false;
		}
	}
	else if (param.compare("simplex") == 0 || param.compare("simulated_annealing") == 0 || param.compare("sa") == 0)
	{
		if (value==0)
		{
			env.simplex = true;
			env.mult = true;
			env.cross = false;
			env.univariate = true;
		}
		else
		{
			env.simplex = false;
		}
	}
	else if (param.compare("xsplit") == 0 || param.compare("cross-split") == 0)
	{
		if (value==0)
		{
			env.cross = false;
		}
		else
		{
			env.cross = true;
			env.bi = false;
			env.mult = false;
		}
	}
	else if (param.compare("lookahead") == 0)
	{
		if (value==0)
		{
			env.lookahead = false;
		}
		else
		{
			env.lookahead = true;
		}
	}
	else if (param.compare("random_state") == 0)
	{
		env.seed(value);
	}
	else if (param.compare("nominal_timelimit") == 0)
	{
		env.nom_timelimit = value;
	}
	else if (param.compare("n_threads") == 0)
	{
		env.n_threads = value;
	}
	else if (param.compare("oblique_split_threshold") == 0)
	{
		env.oblique_split_threshold = value;
	}
	else if (param.compare("cross_entropy") == 0 || param.compare("ce") == 0)
	{
		if(value==0){
			env.cross_entropy = false;
		}
		else
		{
			env.cross_entropy = true;
			env.mult = true;
		}
	}
	else if (param.compare("random_splits") == 0)
	{
		if (value==0)
		{
			env.branch_randomly = false;
		}
		else
		{
			env.branch_randomly = true;
		}
	}
	else if (param.compare("normalize") == 0)
	{
		if (value==0)
		{
			env.normalize = false;
		}
		else
		{
			env.normalize = true;
		}
	}
	else if (param.compare("nominal_partitions") == 0 && value>0)
	{
		env.nominal_partitions = value;
	}
	else if (param.compare("cross_entropy_samples") == 0 || param.compare("ce_samples") == 0)
	{
		if(value>0)
		{
			env.cross_entropy_samples = value;
			env.cross_entropy_samples_form = "";
			env.cross_entropy_samples_expression.destroy();
		}
	}
	else if (param.compare("cross_entropy_no_improvement") == 0)
	{
		if(value>0)
		{
			env.cross_entropy_no_improvement = value;
			env.cross_entropy_no_improvement_form = "";
			env.cross_entropy_no_improvement_expression.destroy();
		}
	}
	else if (param.compare("cross_entropy_alpha") == 0)
	{
		if(value>0 && value <=1){
			env.cross_entropy_alpha = value;
		}
	}
	else if (param.compare("cross_entropy_rho") == 0)
	{
		if(value>0 && value <1){
			env.cross_entropy_rho = value;
		}
	}
	else if (param.compare("check_residual") == 0)
	{
		if (value==0)
		{
			env.check_residual = false;
		}
		else
		{
			env.check_residual = true;
		}
	}
	else if (param.compare("sa_start_iterations") == 0){
		if(value>1)
		{	
			env.spx_start_iterations=value;
			env.spx_start_iterations_form="";
			env.spx_start_iterations_expression.destroy();
		}
	}
	else if (param.compare("sa_no_improvement") == 0){
		if(value>1)
		{
			env.spx_no_improvement=value;
		}
	}
	else if (param.compare("sa_beta") == 0)
	{
		if (value < 1 && value > 0)
		{
			env.spx_beta_T = value;
		}
	}
	else if (param.compare("sa_max_eta") == 0)
	{
		if(value>0)
		{
			env.spx_max_eta == value;
		}
	}
	else
	{
		cout << "Numerical parameter \"" << param << "\" unknown!\n";
	}
}
DecisionTree::DecisionTree()
{
	env = Environment();
}
DecisionTree::~DecisionTree()
{
	for (std::map<int, Node *>::iterator it = nodeList.begin(); it != nodeList.end(); it++)
	{
		int i = it->first;
		delete nodeList[i];
	}
	if (fitted)
	{
		if (root_possesses_data)
		{
			if (!regression)
			{
				delete[] Y;
			}
			else
			{
				delete[] Y_reg;
			}
			if (total_numeric_features > 0)
			{
				for (int j = 0; j < total_numeric_features; j++)
				{
					delete[] X[j];
				}
				delete[] X;
			}
			if (total_nominal_features > 0)
			{
				for (int i = 0; i < total_points; i++)
				{
					delete[] Xnominal[i];
				}
				delete[] Xnominal;
			}
		}
		p->stop(true);
		delete p;
	}
	nodeList.clear();
	adjList.clear();
}
void DecisionTree::removeSubTree(int id)
{
	for (int i = adjList[id].size() - 1; i >= 0; i--)
	{
		int id2 = adjList[id][i]->id;
		removeSubTree(id2);
		adjList[id].erase(adjList[id].begin() + i);
		adjList.erase(id2);
		delete nodeList[id2];
		nodeList.erase(id2);
	}
	std::vector<Node *>().swap(adjList[id]);
	nodeList[id]->isLeaf = true;
}
void DecisionTree::prune(vector<vector<double>> X, vector<int> Y)
{
	bool flag = true;
	int error = evaluate(X, Y);
	while (flag)
	{
		flag = false;
		map<int, Node *>::reverse_iterator rit = nodeList.rbegin();
		int id = 1;
		double min = 0;
		while (rit != nodeList.rend() && !flag)
		{
			int i = rit->first;
			Node *node = nodeList[i];
			if (node->isLeaf)
			{
				rit++;
				continue;
			}
			double reducedError = node->errorC - getSubTreeError(nodeList[i]->id, false, true);
			if (reducedError <= min)
			{
				flag = true;
				id = i;
				min = reducedError;
			}
			rit++;
		}
		if (flag)
		{
			nodeList[id]->isLeaf = true;
			removeSubTree(id);
		}
	}
}
void DecisionTree::prune_reg(vector<vector<double>> X, vector<double> Y)
{
	bool flag = true;
	double error = evaluate_reg(X, Y);
	while (flag)
	{
		flag = false;
		map<int, Node *>::reverse_iterator rit = nodeList.rbegin();
		int id = 1;
		double min = 0;
		while (rit != nodeList.rend() && !flag)
		{
			int i = rit->first;
			Node *node = nodeList[i];
			if (node->isLeaf)
			{
				rit++;
				continue;
			}
			double reducedError = node->errorC - getSubTreeError(nodeList[i]->id, false, true);
			if (reducedError <= min)
			{
				flag = true;
				id = i;
				min = reducedError;
			}
			rit++;
		}
		if (flag)
		{
			nodeList[id]->isLeaf = true;
			removeSubTree(id);
		}
	}
}
int DecisionTree::MCCevaluate(vector<vector<double>> &X, vector<int> &Y2, std::map<int, bool> &T)
{
	vector<vector<double>> X_nominal;
	vector<vector<double>> X2;
	vector<int> isNominal;
	for (int i = 0; i < total_features; i++)
	{
		isNominal.push_back(0);
	}
	for (int i = 0; i < nominal_indices.size(); i++)
	{
		isNominal[nominal_indices[i]] = 1;
	}
	for (int i = 0; i < X.size(); i++)
	{
		vector<double> x_numeric;
		vector<double> x_nominal;
		for (int j = 0; j < total_features; j++)
		{
			if (isNominal[j] == 0)
			{
				x_numeric.push_back(X[i][j]);
			}
			else
			{
				x_nominal.push_back(X[i][j]);
			}
		}
		X_nominal.push_back(x_nominal);
		X2.push_back(x_numeric);
	}
	X = X2;
	int count = 0;
	vector<int> Y;
	for (int i = 0; i < X.size(); i++)
	{
		Node *node = root;
		while (node->isLeaf == false && adjList[node->id][0]->isInTree && T[adjList[node->id][0]->id])
		{
			for (int j = 0; j < adjList[node->id].size(); j++)
			{
				if (adjList[node->id][j]->evaluate(X[i], X_nominal[i]))
				{
					node = adjList[node->id][j];
					break;
				}
			}
		}
		Y.push_back(node->maxVote);
		if (Y[i] != Y2[i])
		{
			count++;
		}
	}
	return count;
}
int DecisionTree::MCCevaluate_reg(vector<vector<double>> &X, vector<double> &Y2, std::map<int, bool> &T)
{
	vector<vector<double>> X_nominal;
	vector<vector<double>> X2;
	vector<int> isNominal;
	for (int i = 0; i < total_features; i++)
	{
		isNominal.push_back(0);
	}
	for (int i = 0; i < nominal_indices.size(); i++)
	{
		isNominal[nominal_indices[i]] = 1;
	}
	for (int i = 0; i < X.size(); i++)
	{
		vector<double> x_numeric;
		vector<double> x_nominal;
		for (int j = 0; j < total_features; j++)
		{
			if (isNominal[j] == 0)
			{
				x_numeric.push_back(X[i][j]);
			}
			else
			{
				x_nominal.push_back(X[i][j]);
			}
		}
		X_nominal.push_back(x_nominal);
		X2.push_back(x_numeric);
	}
	X = X2;
	double count = 0;
	vector<double> Y;
	for (int i = 0; i < X.size(); i++)
	{
		Node *node = root;
		while (node->isLeaf == false && adjList[node->id][0]->isInTree && T[adjList[node->id][0]->id])
		{
			for (int j = 0; j < adjList[node->id].size(); j++)
			{
				if (adjList[node->id][j]->evaluate(X[i], X_nominal[i]))
				{
					node = adjList[node->id][j];
					break;
				}
			}
		}
		Y.push_back(node->mean);
		count += (Y2[i] - Y[i]) * (Y2[i] - Y[i]);
	}
	return count;
}
vector<int> DecisionTree::MCCgetLeaves(int id, std::map<int, bool> &T)
{
	vector<int> qu;
	qu.push_back(id);
	vector<int> l;
	while (qu.size() >= 1)
	{
		int id2 = qu[qu.size() - 1];
		Node *node = nodeList[id2];
		qu.pop_back();
		if (node->isLeaf || !adjList[node->id][0]->isInTree || !T[adjList[node->id][0]->id])
		{
			l.push_back(node->id);
		}
		else
		{
			for (int i = 0; i < adjList[node->id].size(); i++)
			{
				qu.push_back(adjList[node->id][i]->id);
			}
		}
	}
	return l;
}
vector<int> DecisionTree::MCCgetSubTree(int id, std::map<int, bool> &T, bool inc)
{
	vector<int> qu;
	qu.push_back(id);
	vector<int> l;
	if (inc)
	{
		l.push_back(id);
	}
	while (qu.size() >= 1)
	{
		int id2 = qu[qu.size() - 1];
		Node *node = nodeList[id2];
		qu.pop_back();
		if (node->id != id)
		{
			l.push_back(node->id);
		}
		if (node->isLeaf || !adjList[node->id][0]->isInTree || !T[adjList[node->id][0]->id])
		{
		}
		else
		{
			for (int i = 0; i < adjList[node->id].size(); i++)
			{
				qu.push_back(adjList[node->id][i]->id);
			}
		}
	}
	return l;
}
vector<int> DecisionTree::MCCgetParents(int id, std::map<int, int> &parent)
{
	vector<int> p;
	Node *node = nodeList[id];
	while (parent[node->id] != -1)
	{
		p.push_back(parent[node->id]);
		node = nodeList[parent[node->id]];
	}
	return p;
}
int DecisionTree::MCCgetSubTreeError(vector<int> ids, std::map<int, bool> &T, std::map<int, double> &error)
{
	int err = 0;
	for (int i = 0; i < ids.size(); i++)
	{
		int id = ids[i];
		if (nodeList[id]->isLeaf || adjList[id][0]->isInTree == false || T[adjList[id][0]->id] == false)
		{
			err += error[id];
		}
	}
	return err;
}
void DecisionTree::MCCprune(vector<vector<double>> X, vector<int> Y)
{
	vector<std::map<int, bool>> T;
	std::map<int, bool> T_tmp;
	std::map<int, double> error;
	std::map<int, double> r;
	std::map<int, double> p;
	std::map<int, double> g;
	std::map<int, int> parent;
	map<int, bool> recompute;
	parent[root->id] = -1;
	for (std::map<int, Node *>::iterator it = nodeList.begin(); it != nodeList.end(); it++)
	{
		int i = it->first;
		recompute[i] = true;
		if (!nodeList[i]->isLeaf)
		{
			for (int j = 0; j < adjList[i].size(); j++)
			{
				parent[adjList[i][j]->id] = i;
			}
		}
	}
	int N = 0;
	for (std::map<int, Node *>::iterator it = nodeList.begin(); it != nodeList.end(); it++)
	{
		int i = it->first;
		g[i] = 0;
		if (nodeList[i]->isInTree)
		{
			T_tmp[i] = true;
			double N_tmp = 0;
			for (int k = 0; k < n_classes; k++)
			{
				if (nodeList[i]->isLeaf)
				{
					N += nodeList[i]->Ni[k];
				}
				N_tmp += nodeList[i]->Ni[k];
			}
			if (N_tmp != 0)
			{
				double max_tmp = 0;
				int max_k = 0;
				for (int k = 0; k < n_classes; k++)
				{
					if (nodeList[i]->Ni[k] > max_tmp)
					{
						max_tmp = nodeList[i]->Ni[k];
						max_k = k;
					}
				}
				if (N_tmp != 0)
				{
					nodeList[i]->maxVote = max_k;
					error[i] = N_tmp - max_tmp;
					r[i] = error[i] / N_tmp;
					p[i] = N_tmp;
				}
				else
				{
					error[i] = 0;
					r[i] = 0;
					p[i] = 0;
				}
			}
		}
		else
		{
			T_tmp[i] = false;
		}
	}
	for (std::map<int, Node *>::iterator it = nodeList.begin(); it != nodeList.end(); it++)
	{
		int i = it->first;
		p[i] = p[i] / N;
	}
	T.push_back(T_tmp);
	vector<double> alphas;
	vector<double> prune_ids;
	vector<double> miss_error;
	double minerror = std::numeric_limits<double>::max();
	vector<int> cut_ids;
	vector<int> ids;
	vector<int> l;
	vector<int> subtree;
	vector<int> parents_tmp;
	double alpha = std::numeric_limits<double>::max();
	int prune_id = 0;
	int nodes = 0;
	std::clock_t start;
	int counter = 0;
	while (true)
	{
		counter += 1;
		start = clock();
		ids.clear();
		int best_k = 0;
		miss_error.push_back(MCCevaluate(X, Y, T_tmp));
		if (miss_error[miss_error.size() - 1] <= minerror)
		{
			minerror = miss_error[miss_error.size() - 1];
			cut_ids.clear();
			for (int i = 0; i < prune_ids.size(); i++)
			{
				cut_ids.push_back(prune_ids[i]);
			}
		}
		alpha = std::numeric_limits<double>::max();
		prune_id = 0;
		nodes = 0;
		ids = MCCgetSubTree(root->id, T_tmp, true);
		if (ids.size() <= adjList[root->id].size() + 1)
		{
			break;
		}
		double duration = 0;
		bool prune_all = true;
		for (int h = 0; h < ids.size(); h++)
		{
			int i = ids[h];
			if (T_tmp[i] && !nodeList[i]->isLeaf && adjList[i][0]->isInTree && T_tmp[adjList[i][0]->id])
			{
				double R = 0;
				l = MCCgetLeaves(i, T_tmp);
				for (int j = 0; j < l.size(); j++)
				{
					R += r[l[j]] * p[l[j]];
				}
				g[i] = (r[i] * p[i] - R) / (l.size() - 1);
				recompute[i] = false;
				if (g[i] < alpha && i != root->id)
				{
					alpha = g[i];
					prune_id = i;
				}
				else if (g[i] == alpha && !prune_all && i != root->id)
				{
					if (MCCgetSubTree(i, T_tmp, false).size() < MCCgetSubTree(prune_id, T_tmp, false).size())
					{
						alpha = g[i];
						prune_id = i;
					}
				}
			}
		}
		if (prune_all)
		{
			for (int h = 0; h < ids.size(); h++)
			{
				int i = ids[h];
				if (T_tmp[i] && !nodeList[i]->isLeaf && adjList[i][0]->isInTree && T_tmp[adjList[i][0]->id])
				{
					if (g[i] == alpha && i != root->id)
					{
						parents_tmp = MCCgetParents(i, parent);
						for (int g = 0; g < parents_tmp.size(); g++)
						{
							recompute[parents_tmp[g]] = true;
						}
						prune_ids.push_back(i);
						subtree = MCCgetSubTree(i, T_tmp, false);
						for (int j = 0; j < subtree.size(); j++)
						{
							T_tmp[subtree[j]] = false;
						}
					}
				}
			}
		}
		else
		{
			parents_tmp = MCCgetParents(prune_id, parent);
			for (int g = 0; g < parents_tmp.size(); g++)
			{
				recompute[parents_tmp[g]] = true;
			}
			prune_ids.push_back(prune_id);
			subtree = MCCgetSubTree(prune_id, T_tmp, false);
			for (int j = 0; j < subtree.size(); j++)
			{
				T_tmp[subtree[j]] = false;
			}
		}
		if (prune_id != 0)
		{
			alphas.push_back(alpha);
		}
	}
	for (int k = 0; k < cut_ids.size(); k++)
	{
		removeSubTree(cut_ids[k]);
	}
}
void DecisionTree::MCCprune_reg(vector<vector<double>> &X, vector<double> &Y)
{
	vector<std::map<int, bool>> T;
	std::map<int, bool> T_tmp;
	std::map<int, double> error;
	std::map<int, double> r;
	std::map<int, double> p;
	std::map<int, double> g;
	std::map<int, int> parent;
	map<int, bool> recompute;
	parent[root->id] = -1;
	for (std::map<int, Node *>::iterator it = nodeList.begin(); it != nodeList.end(); it++)
	{
		int i = it->first;
		recompute[i] = true;
		if (!nodeList[i]->isLeaf)
		{
			for (int j = 0; j < adjList[i].size(); j++)
			{
				parent[adjList[i][j]->id] = i;
			}
		}
	}
	int N = 0;
	for (std::map<int, Node *>::iterator it = nodeList.begin(); it != nodeList.end(); it++)
	{
		int i = it->first;
		g[i] = 0;
		if (nodeList[i]->isInTree)
		{
			T_tmp[i] = true;
			double N_tmp = 0;
			if (nodeList[i]->isLeaf)
			{
				N += nodeList[i]->total_points;
			}
			N_tmp += nodeList[i]->total_points;
			if (N_tmp != 0)
			{
				error[i] = N_tmp * nodeList[i]->mse;
				r[i] = error[i] / N_tmp;
				p[i] = N_tmp;
			}
			else
			{
				error[i] = 0;
				r[i] = 0;
				p[i] = 0;
			}
		}
		else
		{
			T_tmp[i] = false;
		}
	}
	for (std::map<int, Node *>::iterator it = nodeList.begin(); it != nodeList.end(); it++)
	{
		int i = it->first;
		p[i] = p[i] / N;
	}
	T.push_back(T_tmp);
	vector<double> alphas;
	vector<double> prune_ids;
	vector<double> miss_error;
	double minerror = std::numeric_limits<double>::max();
	vector<int> cut_ids;
	vector<int> ids;
	vector<int> l;
	vector<int> subtree;
	vector<int> parents_tmp;
	double alpha = std::numeric_limits<double>::max();
	int prune_id = 0;
	int nodes = 0;
	std::clock_t start;
	while (true)
	{
		start = clock();
		ids.clear();
		int best_k = 0;
		miss_error.push_back(MCCevaluate_reg(X, Y, T_tmp));
		if (miss_error[miss_error.size() - 1] <= minerror)
		{
			minerror = miss_error[miss_error.size() - 1];
			cut_ids.clear();
			for (int i = 0; i < prune_ids.size(); i++)
			{
				cut_ids.push_back(prune_ids[i]);
			}
		}
		alpha = std::numeric_limits<double>::max();
		prune_id = 0;
		nodes = 0;
		ids = MCCgetSubTree(root->id, T_tmp, true);
		if (ids.size() <= adjList[root->id].size() + 1)
		{
			break;
		}
		double duration = 0;
		bool prune_all = true;
		for (int h = 0; h < ids.size(); h++)
		{
			int i = ids[h];
			if (T_tmp[i] && !nodeList[i]->isLeaf && adjList[i][0]->isInTree && T_tmp[adjList[i][0]->id])
			{
				if (recompute[i])
				{
					double R = 0;
					std::clock_t start2 = clock();
					l = MCCgetLeaves(i, T_tmp);
					duration += (std::clock() - start2);
					for (int j = 0; j < l.size(); j++)
					{
						R += r[l[j]] * p[l[j]];
					}
					g[i] = (r[i] * p[i] - R) / (l.size() - 1);
					recompute[i] = false;
				}
				if (g[i] < alpha && i != root->id)
				{
					alpha = g[i];
					prune_id = i;
				}
				else if (g[i] == alpha && !prune_all && i != root->id)
				{
					if (MCCgetSubTree(i, T_tmp, false).size() < MCCgetSubTree(prune_id, T_tmp, false).size())
					{
						alpha = g[i];
						prune_id = i;
					}
				}
			}
		}
		if (prune_all)
		{
			for (int h = 0; h < ids.size(); h++)
			{
				int i = ids[h];
				if (T_tmp[i] && !nodeList[i]->isLeaf && adjList[i][0]->isInTree && T_tmp[adjList[i][0]->id])
				{
					if (g[i] == alpha && i != root->id)
					{
						parents_tmp = MCCgetParents(i, parent);
						for (int g = 0; g < parents_tmp.size(); g++)
						{
							recompute[parents_tmp[g]] = true;
						}
						prune_ids.push_back(i);
						subtree = MCCgetSubTree(i, T_tmp, false);
						for (int j = 0; j < subtree.size(); j++)
						{
							T_tmp[subtree[j]] = false;
						}
					}
				}
			}
		}
		else
		{
			parents_tmp = MCCgetParents(prune_id, parent);
			for (int g = 0; g < parents_tmp.size(); g++)
			{
				recompute[parents_tmp[g]] = true;
			}
			prune_ids.push_back(prune_id);
			subtree = MCCgetSubTree(prune_id, T_tmp, false);
			for (int j = 0; j < subtree.size(); j++)
			{
				T_tmp[subtree[j]] = false;
			}
		}
		if (prune_id != 0)
		{
			alphas.push_back(alpha);
		}
	}
	for (int k = 0; k < cut_ids.size(); k++)
	{
		removeSubTree(cut_ids[k]);
	}
}
int DecisionTree::evaluate(vector<vector<double>> X, vector<int> Y2)
{
	vector<vector<double>> X_nominal;
	vector<vector<double>> X2;
	vector<int> isNominal;
	for (int i = 0; i < total_features; i++)
	{
		isNominal.push_back(0);
	}
	for (int i = 0; i < nominal_indices.size(); i++)
	{
		isNominal[nominal_indices[i]] = 1;
	}
	for (int i = 0; i < X.size(); i++)
	{
		vector<double> x_numeric;
		vector<double> x_nominal;
		for (int j = 0; j < total_features; j++)
		{
			if (isNominal[j] == 0)
			{
				x_numeric.push_back(X[i][j]);
			}
			else
			{
				x_nominal.push_back(X[i][j]);
			}
		}
		X_nominal.push_back(x_nominal);
		X2.push_back(x_numeric);
	}
	X = X2;
	int count = 0;
	vector<int> Y;
	for (std::map<int, Node *>::iterator it = nodeList.begin(); it != nodeList.end(); it++)
	{
		int i = it->first;
		nodeList[i]->errorC = 0;
	}
	for (int i = 0; i < X.size(); i++)
	{
		Node *node = root;
		while (node->isLeaf == false)
		{
			if (Y2[i] != node->maxVote)
			{
				node->errorC++;
			}
			for (int j = 0; j < adjList[node->id].size(); j++)
			{
				adjList[node->id][j];
				if (adjList[node->id][j]->evaluate(X[i], X_nominal[i]))
				{
					node = adjList[node->id][j];
					break;
				}
			}
		}
		if (Y2[i] != node->maxVote)
		{
			node->errorC++;
		}
		Y.push_back(node->maxVote);
		if (Y[i] != Y2[i])
		{
			count++;
		}
	}
	return count;
}
double DecisionTree::evaluate_reg(vector<vector<double>> X, vector<double> Y2)
{
	vector<vector<double>> X_nominal;
	vector<vector<double>> X2;
	vector<int> isNominal;
	for (int i = 0; i < total_features; i++)
	{
		isNominal.push_back(0);
	}
	for (int i = 0; i < nominal_indices.size(); i++)
	{
		isNominal[nominal_indices[i]] = 1;
	}
	for (int i = 0; i < X.size(); i++)
	{
		vector<double> x_numeric;
		vector<double> x_nominal;
		for (int j = 0; j < total_features; j++)
		{
			if (isNominal[j] == 0)
			{
				x_numeric.push_back(X[i][j]);
			}
			else
			{
				x_nominal.push_back(X[i][j]);
			}
		}
		X_nominal.push_back(x_nominal);
		X2.push_back(x_numeric);
	}
	X = X2;
	double err = 0;
	vector<int> Y;
	for (std::map<int, Node *>::iterator it = nodeList.begin(); it != nodeList.end(); it++)
	{
		int i = it->first;
		nodeList[i]->errorC = 0;
	}
	for (int i = 0; i < X.size(); i++)
	{
		Node *node = root;
		while (node->isLeaf == false)
		{
			node->errorC += (Y2[i] - node->mean) * (Y2[i] - node->mean);
			for (int j = 0; j < adjList[node->id].size(); j++)
			{
				adjList[node->id][j];
				if (adjList[node->id][j]->evaluate(X[i], X_nominal[i]))
				{
					node = adjList[node->id][j];
					break;
				}
			}
		}
		node->errorC += (Y2[i] - node->mean) * (Y2[i] - node->mean);
		Y.push_back(node->mean);
		err += (Y2[i] - Y[i]) * (Y2[i] - Y[i]);
	}
	return err;
}
vector<int> DecisionTree::predict(vector<vector<double>> X)
{
	vector<vector<double>> X_nominal;
	vector<vector<double>> X2;
	vector<int> isNominal;
	for (int i = 0; i < total_features; i++)
	{
		isNominal.push_back(0);
	}
	for (int i = 0; i < nominal_indices.size(); i++)
	{
		isNominal[nominal_indices[i]] = 1;
	}
	for (int i = 0; i < X.size(); i++)
	{
		vector<double> x_numeric;
		vector<double> x_nominal;
		for (int j = 0; j < total_features; j++)
		{
			if (isNominal[j] == 0)
			{
				x_numeric.push_back(X[i][j]);
			}
			else
			{
				x_nominal.push_back(X[i][j]);
			}
		}
		X_nominal.push_back(x_nominal);
		X2.push_back(x_numeric);
	}
	X = X2;
	int count = 0;
	vector<int> Y;
	for (int i = 0; i < X.size(); i++)
	{
		Node *node = root;
		while (node->isLeaf == false)
		{
			for (int j = 0; j < adjList[node->id].size(); j++)
			{
				if (adjList[node->id][j]->evaluate(X[i], X_nominal[i]))
				{
					node = adjList[node->id][j];
					break;
				}
			}
		}
		Y.push_back(node->maxVote);
	}
	return Y;
}
vector<double> DecisionTree::predict_reg(vector<vector<double>> X)
{
	vector<vector<double>> X_nominal;
	vector<vector<double>> X2;
	vector<int> isNominal;
	for (int i = 0; i < total_features; i++)
	{
		isNominal.push_back(0);
	}
	for (int i = 0; i < nominal_indices.size(); i++)
	{
		isNominal[nominal_indices[i]] = 1;
	}
	for (int i = 0; i < X.size(); i++)
	{
		vector<double> x_numeric;
		vector<double> x_nominal;
		for (int j = 0; j < total_features; j++)
		{
			if (isNominal[j] == 0)
			{
				x_numeric.push_back(X[i][j]);
			}
			else
			{
				x_nominal.push_back(X[i][j]);
			}
		}
		X_nominal.push_back(x_nominal);
		X2.push_back(x_numeric);
	}
	X = X2;
	int count = 0;
	vector<double> Y;
	for (int i = 0; i < X.size(); i++)
	{
		Node *node = root;
		while (node->isLeaf == false)
		{
			for (int j = 0; j < adjList[node->id].size(); j++)
			{
				if (adjList[node->id][j]->evaluate(X[i], X_nominal[i]))
				{
					node = adjList[node->id][j];
					break;
				}
			}
		}
		if (env.criterion == Environment::mse)
		{
			Y.push_back(node->mean);
		}
		else
		{
			Y.push_back(node->median);
		}
	}
	return Y;
}
void DecisionTree::setAttributeNames(vector<string> names)
{
	this->names = names;
}
void DecisionTree::setNominalNames(int index, vector<int> values, vector<string> categories)
{
	map<int, string> tmp;
	if (values.size() == categories.size())
	{
		for (int j = 0; j < values.size(); j++)
		{
			tmp[values[j]] = categories[j];
		}
	}
	cat_to_name[index] = tmp;
}
void DecisionTree::setYNames(vector<int> values, vector<string> categories)
{
	if (values.size() == categories.size())
	{
		for (int j = 0; j < values.size(); j++)
		{
			y_to_name.insert(pair<int, string>(values[j], categories[j]));
		}
	}
}
void DecisionTree::fit(vector<vector<double>> x, vector<int> y, vector<int> nominal)
{
	if (fitted)
	{
		cout << "Decision Tree is already fitted\n";
		return;
	}
	fitted = true;
	regression = false;
	env.regression = false;
	this->nominal_indices = nominal;
	p = new ctpl::thread_pool(env.n_threads);
	depth = 0;
	vector<int> isNominal;
	total_features = x[0].size();
	total_numeric_features = x[0].size() - nominal.size();
	total_nominal_features = nominal.size();
	for (int i = 0; i < total_features; i++)
	{
		isNominal.push_back(0);
	}
	for (int i = 0; i < nominal.size(); i++)
	{
		isNominal[nominal[i]] = 1;
	}
	if (names.size() == 0 || names.size() != total_features)
	{
		for (int j = 0; j < total_features; j++)
		{
			if (isNominal[j])
			{
				nominalNames.push_back("x" + std::string("<sub>") + itos(j + 1) + "</sub>");
			}
			else
			{
				numericNames.push_back("x" + std::string("<sub>") + itos(j + 1) + "</sub>");
			}
		}
	}
	else
	{
		for (int j = 0; j < total_features; j++)
		{
			if (isNominal[j])
			{
				nominalNames.push_back(names[j]);
			}
			else
			{
				numericNames.push_back(names[j]);
			}
		}
	}
	total_points = x.size();
	if (env.min_samples_split < 1.0)
	{
		env.min_samples_split = round(env.min_samples_split * total_points);
	}
	if (env.min_cases < 1.0)
	{
		env.min_cases = round(env.min_cases * total_points);
	}
	n_classes = 0;
	if (env.max_features > total_features)
	{
		env.max_features = total_features;
	}
	Y = new int[total_points];
	if (nominal.size() < total_features)
	{
		X = new double *[total_numeric_features];
		for (int j = 0; j < total_numeric_features; j++)
		{
			X[j] = new double[total_points];
		}
		int j2 = 0;
		for (int j = 0; j < total_features; j++)
		{
			if (isNominal[j] == 0)
			{
				for (int i = 0; i < total_points; i++)
				{
					X[j2][i] = x[i][j];
				}
				j2++;
			}
		}
	}
	if (nominal.size() > 0)
	{
		Xnominal = new double *[total_points];
		for (int i = 0; i < total_points; i++)
		{
			Xnominal[i] = new double[nominal.size()];
			for (int j2 = 0; j2 < nominal.size(); j2++)
			{
				Xnominal[i][j2] = x[i][nominal[j2]];
			}
		}
		for (int j = 0; j < nominal.size(); j++)
		{
			int max = 0;
			for (int i = 0; i < total_points; i++)
			{
				if (Xnominal[i][j] > max)
				{
					max = Xnominal[i][j];
				}
			}
			nominal_domain.push_back(max);
		}
	}
	for (int j = 0; j < total_nominal_features; j++)
	{
		vector<bool> d(nominal_domain[j] + 1, false);
		is_nominal_val_in_dataset.push_back(d);
	}
	for (int j = 0; j < total_nominal_features; j++)
	{
		for (int i = 0; i < total_points; i++)
		{
			int val = (int)Xnominal[i][j];
			is_nominal_val_in_dataset[j][val] = true;
		}
	}
	map<int, map<int, string>> cat_to_name_tmp;
	for (int j = 0; j < nominal.size(); j++)
	{
		for (int v = 0; v <= nominal_domain[j]; v++)
		{
			cat_to_name_tmp[j][v] = itos(v);
		}
		for (map<int, string>::iterator it = cat_to_name[nominal[j]].begin(); it != cat_to_name[nominal[j]].end(); it++)
		{
			cat_to_name_tmp[j][it->first] = it->second;
		}
	}
	cat_to_name = cat_to_name_tmp;
	map<int, string> y_to_name_tmp;
	for (int i = 0; i < total_points; i++)
	{
		if (y[i] + 1 >= n_classes)
		{
			n_classes = y[i] + 1;
		}
		Y[i] = y[i];
		y_to_name_tmp[y[i]] = itos(y[i]);
	}
	for (map<int, string>::iterator it = y_to_name.begin(); it != y_to_name.end(); it++)
	{
		if (y_to_name.find(it->first) == y_to_name.end())
		{
		}
		else
		{
			y_to_name_tmp[it->first] = y_to_name[it->first];
		}
	}
	y_to_name = y_to_name_tmp;
	root = new Node(X, Xnominal, Y, Y_reg, total_points, total_numeric_features, total_nominal_features, nominal_domain, is_nominal_val_in_dataset, n_classes, p, &env);
	vector<double> a;
	for (int j = 0; j < total_numeric_features; j++)
	{
		a.push_back(0);
	}
	root->cond = Condition(a, 0, Condition::t);
	root->depth = 0;
	int count = 0;
	root->id = count;
	adjList[count] = vector<Node *>();
	nodeList[count] = root;
	leaves.push_back(root);
	count++;
	n_leaf_nodes = 1;
	if (total_points < env.min_samples_split)
	{
		return;
	}
	bool flag = false;
	for (int t = 0; t < env.max_depth; t++)
	{
		bool inc_depth = true;
		int n_leaves = leaves.size();
		if (n_leaves == 0 || flag)
		{
			break;
		}
		vector<Node *> leaves2;
		for (int i = 0; i < n_leaves; i++)
		{
			int n_leaves2 = leaves.size();
			Node *node = leaves.back();
			leaves.pop_back();
			if (node->total_points < env.min_samples_split)
			{
				node->isLeaf = true;
				continue;
			}
			vector<Node *> nodes = node->branch();
			if (n_leaf_nodes + nodes.size() - 1 > env.max_leaf_nodes)
			{
				for (int l = 0; l < nodes.size(); l++)
				{
					delete nodes[l];
				}
				node->isLeaf = true;
				continue;
			}
			if (node->impurity < env.min_impurity_split)
			{
				for (int l = 0; l < nodes.size(); l++)
				{
					delete nodes[l];
				}
				node->isLeaf = true;
				continue;
			}
			if (nodes.size() == 0)
			{
				node->isLeaf = true;
				continue;
			}
			double impurity_decrease = node->impurity / node->total_points;
			for (int j = 0; j < nodes.size(); j++)
			{
				impurity_decrease -= (nodes[j]->total_points / (double)node->total_points) * (nodes[j]->impurity / nodes[j]->total_points);
			}
			impurity_decrease = (node->total_points / (double)total_points) * impurity_decrease;
			if (impurity_decrease < env.min_impurity_decrease)
			{
				for (int l = 0; l < nodes.size(); l++)
				{
					delete nodes[l];
				}
				node->isLeaf = true;
				continue;
			}
			int c2 = 0;
			for (int j = 0; j < nodes.size(); j++)
			{
				if (nodes[j]->total_points >= env.min_cases && nodes[j]->total_points > 0)
				{
					c2++;
				}
			}
			if (c2 < 2)
			{
				for (int l = 0; l < nodes.size(); l++)
				{
					delete nodes[l];
				}
				node->isLeaf = true;
				continue;
			}
			if (inc_depth)
			{
				depth++;
				inc_depth = false;
			}
			n_leaf_nodes--;
			node->isLeaf = false;
			for (int j = 0; j < nodes.size(); j++)
			{
				nodes[j]->id = count;
				nodes[j]->depth = depth;
				adjList[count] = vector<Node *>();
				nodeList[count] = nodes[j];
				count++;
				adjList[node->id].push_back(nodes[j]);
				n_leaf_nodes++;
				if (nodes[j]->shouldSplit)
				{
					leaves2.push_back(nodes[j]);
				}
			}
		}
		reverse(leaves2.begin(), leaves2.end());
		if (env.max_leaf_nodes != std::numeric_limits<int>::max())
		{
			shuffle(leaves2.begin(), leaves2.end(), env.generator);
		}
		leaves = leaves2;
	}
	leaves.clear();
	map<int, Node *>::iterator it = nodeList.begin();
	while (it != nodeList.end())
	{
		int i = it->first;
		nodeList[i]->cond.names = numericNames;
		nodeList[i]->cond.nominalNames = nominalNames;
		nodeList[i]->cond.cat_to_name = cat_to_name;
		nodeList[i]->maxVote_str = y_to_name[nodeList[i]->maxVote];
		if (nodeList[i]->isLeaf)
		{
			leaves.push_back(nodeList[i]);
		}
		it++;
	}
}
void DecisionTree::fit(vector<vector<double>> x, vector<int> y)
{
	vector<int> nominal_tmp;
	fit(x, y, nominal_tmp);
}
void DecisionTree::fit_reg(vector<vector<double>> x, vector<double> y, vector<int> nominal)
{
	if (fitted)
	{
		cout << "Decision Tree is already fitted\n";
		return;
	}
	env.criterion = env.criterion_reg;
	fitted = true;
	regression = true;
	env.regression = true;
	this->nominal_indices = nominal;
	p = new ctpl::thread_pool(env.n_threads);
	depth = 0;
	vector<int> isNominal;
	total_features = x[0].size();
	total_numeric_features = x[0].size() - nominal.size();
	total_nominal_features = nominal.size();
	for (int i = 0; i < total_features; i++)
	{
		isNominal.push_back(0);
	}
	for (int i = 0; i < nominal.size(); i++)
	{
		isNominal[nominal[i]] = 1;
	}
	if (names.size() == 0 || names.size() != total_features)
	{
		for (int j = 0; j < total_features; j++)
		{
			if (isNominal[j])
			{
				nominalNames.push_back("x" + std::string("<sub>") + itos(j + 1) + "</sub>");
			}
			else
			{
				numericNames.push_back("x" + std::string("<sub>") + itos(j + 1) + "</sub>");
			}
		}
	}
	else
	{
		for (int j = 0; j < total_features; j++)
		{
			if (isNominal[j])
			{
				nominalNames.push_back(names[j]);
			}
			else
			{
				numericNames.push_back(names[j]);
			}
		}
	}
	total_points = x.size();
	if (env.min_samples_split < 1.0)
	{
		env.min_samples_split = round(env.min_samples_split * total_points);
	}
	if (env.min_cases < 1.0)
	{
		env.min_cases = round(env.min_cases * total_points);
	}
	n_classes = 0;
	if (env.max_features > total_features)
	{
		env.max_features = total_features;
	}
	Y_reg = new double[total_points];
	if (nominal.size() < total_features)
	{
		X = new double *[total_numeric_features];
		for (int j = 0; j < total_numeric_features; j++)
		{
			X[j] = new double[total_points];
		}
		for (int i = 0; i < total_points; i++)
		{
			int j2 = 0;
			for (int j = 0; j < total_features; j++)
			{
				if (isNominal[j] == 0)
				{
					X[j2][i] = x[i][j];
					j2++;
				}
			}
		}
	}
	if (nominal.size() > 0)
	{
		Xnominal = new double *[total_points];
		for (int i = 0; i < total_points; i++)
		{
			Xnominal[i] = new double[nominal.size()];
			for (int j2 = 0; j2 < nominal.size(); j2++)
			{
				Xnominal[i][j2] = x[i][nominal[j2]];
			}
		}
		for (int j = 0; j < nominal.size(); j++)
		{
			int max = 0;
			for (int i = 0; i < total_points; i++)
			{
				if (Xnominal[i][j] > max)
				{
					max = Xnominal[i][j];
				}
			}
			nominal_domain.push_back(max);
		}
	}
	for (int j = 0; j < total_nominal_features; j++)
	{
		vector<bool> d(nominal_domain[j] + 1, false);
		is_nominal_val_in_dataset.push_back(d);
	}
	for (int j = 0; j < total_nominal_features; j++)
	{
		for (int i = 0; i < total_points; i++)
		{
			int val = (int)Xnominal[i][j];
			is_nominal_val_in_dataset[j][val] = true;
		}
	}
	map<int, map<int, string>> cat_to_name_tmp;
	for (int j = 0; j < nominal.size(); j++)
	{
		for (int v = 0; v <= nominal_domain[j]; v++)
		{
			cat_to_name_tmp[j][v] = itos(v);
		}
		for (map<int, string>::iterator it = cat_to_name[nominal[j]].begin(); it != cat_to_name[nominal[j]].end(); it++)
		{
			cat_to_name_tmp[j][it->first] = it->second;
		}
	}
	cat_to_name = cat_to_name_tmp;
	for (int i = 0; i < total_points; i++)
	{
		Y_reg[i] = y[i];
	}
	root = new Node(X, Xnominal, Y, Y_reg, total_points, total_numeric_features, total_nominal_features, nominal_domain, is_nominal_val_in_dataset, n_classes, p, &env);
	vector<double> a;
	for (int j = 0; j < total_numeric_features; j++)
	{
		a.push_back(0);
	}
	root->cond = Condition(a, 0, Condition::t);
	root->depth = 0;
	int count = 0;
	root->id = count;
	adjList[count] = vector<Node *>();
	nodeList[count] = root;
	leaves.push_back(root);
	count++;
	n_leaf_nodes = 1;
	if (total_points < env.min_samples_split)
	{
		return;
	}
	bool flag = false;
	for (int t = 0; t < env.max_depth; t++)
	{
		bool inc_depth = true;
		int n_leaves = leaves.size();
		if (n_leaves == 0 || flag)
		{
			break;
		}
		vector<Node *> leaves2;
		for (int i = 0; i < n_leaves; i++)
		{
			int n_leaves2 = leaves.size();
			Node *node = leaves.back();
			leaves.pop_back();
			if (node->total_points < env.min_samples_split)
			{
				node->isLeaf = true;
				continue;
			}
			vector<Node *> nodes = node->branch();
			if (n_leaf_nodes + nodes.size() - 1 > env.max_leaf_nodes)
			{
				node->isLeaf = true;
				continue;
			}
			if (node->mse < env.min_impurity_split)
			{
				node->isLeaf = true;
				continue;
			}
			if (nodes.size() == 0)
			{
				node->isLeaf = true;
				continue;
			}
			double impurity_decrease = node->mse;
			for (int j = 0; j < nodes.size(); j++)
			{
				impurity_decrease -= (nodes[j]->total_points / (double)node->total_points) * nodes[j]->mse;
			}
			impurity_decrease = (node->total_points / (double)total_points) * impurity_decrease;
			if (impurity_decrease < env.min_impurity_decrease)
			{
				node->isLeaf = true;
				continue;
			}
			int c2 = 0;
			for (int j = 0; j < nodes.size(); j++)
			{
				if (nodes[j]->total_points >= env.min_cases && nodes[j]->total_points > 0)
				{
					c2++;
				}
			}
			if (c2 < 2)
			{
				node->isLeaf = true;
				continue;
			}
			if (inc_depth)
			{
				depth++;
				inc_depth = false;
			}
			n_leaf_nodes--;
			node->isLeaf = false;
			for (int j = 0; j < nodes.size(); j++)
			{
				nodes[j]->id = count;
				nodes[j]->depth = depth;
				adjList[count] = vector<Node *>();
				nodeList[count] = nodes[j];
				count++;
				adjList[node->id].push_back(nodes[j]);
				n_leaf_nodes++;
				if (nodes[j]->shouldSplit)
				{
					leaves2.push_back(nodes[j]);
				}
			}
		}
		reverse(leaves2.begin(), leaves2.end());
		if (env.max_leaf_nodes != std::numeric_limits<int>::max())
		{
			shuffle(leaves2.begin(), leaves2.end(), env.generator);
		}
		leaves = leaves2;
	}
	leaves.clear();
	map<int, Node *>::iterator it = nodeList.begin();
	while (it != nodeList.end())
	{
		int i = it->first;
		nodeList[i]->cond.names = numericNames;
		nodeList[i]->cond.nominalNames = nominalNames;
		nodeList[i]->cond.cat_to_name = cat_to_name;
		nodeList[i]->maxVote_str = y_to_name[nodeList[i]->maxVote];
		if (nodeList[i]->isLeaf)
		{
			leaves.push_back(nodeList[i]);
		}
		it++;
	}
}
void DecisionTree::fit_reg(vector<vector<double>> x, vector<double> y)
{
	vector<int> nominal_tmp;
	fit_reg(x, y, nominal_tmp);
}
double DecisionTree::getError(int id)
{
	double e = 0;
	double f = 0;
	double zs[10] = {4.0, 3.09, 2.58, 2.33, 1.65, 1.28, 0.84, 0.69, 0.25, 0.0};
	double cs[10] = {0, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.25, 0.4, 1.0};
	int ind = 0;
	for (int i = 0; i < 9; i++)
	{
		if (env.cf >= cs[i] && env.cf <= cs[i + 1])
		{
			if (env.cf < (cs[i] + cs[i + 1]) / 2)
			{
				ind = i;
			}
			else
			{
				ind = i + 1;
			}
		}
	}
	double z = zs[ind];
	Node *node = nodeList[id];
	double N = node->total_points;
	if (N == 0)
	{
		return 0;
	}
	if (!regression)
	{
		for (int m = 0; m < n_classes; m++)
		{
			if (m != node->maxVote)
			{
				f += node->Ni[m];
			}
		}
	}
	else
	{
		if (env.criterion == Environment::mae)
		{
			f += node->mae * node->total_points;
		}
		else
		{
			f += node->mse * node->total_points;
		}
	}
	f = f / N;
	e = (f + (z * z) / (2 * N) + z * sqrt(f / N - (f * f) / N + (z * z) / (4 * N * N))) / (1 + (z * z) / N);
	return e;
}
double DecisionTree::getSubTreeError(int id, bool collect, bool eval)
{
	Node *node = nodeList[id];
	double e = 0;
	if (node->isLeaf)
	{
		if (collect == false)
		{
			return 0;
		}
		else
		{
			if (!eval)
			{
				return (node->total_points * getError(node->id));
			}
			else
			{
				return (node->errorC);
			}
		}
	}
	for (int i = 0; i < adjList[id].size(); i++)
	{
		e += getSubTreeError(adjList[id][i]->id, true, eval);
	}
	if (collect)
	{
		return e;
	}
	else
	{
		if (!eval)
		{
			return (e / node->total_points);
		}
		else
		{
			return e;
		}
	}
}
string DecisionTree::export_graphviz(bool on_edge)
{
	string out = "digraph Tree {\n node [shape=box] ; \n";
	map<int, Node *>::iterator it2 = nodeList.begin();
	for (it2 = nodeList.begin(); it2 != nodeList.end(); it2++)
	{
		if (nodeList[it2->first]->cond.isCross)
		{
			on_edge = true;
		}
		if (nodeList[it2->first]->cond.isNominal)
		{
			on_edge = true;
		}
	}
	it2 = nodeList.begin();
	while (it2 != nodeList.end())
	{
		int i = it2->first;
		if (!nodeList[i]->isInTree)
		{
			it2++;
			continue;
		}
		out += itos(i);
		out += " [label=<";
		out += string("id = ") + itos(nodeList[i]->id);
		if (!regression)
		{
			if (!(env.criterion == env.twoing))
			{
				out += "<br/>";
				if (env.criterion == env.gini_imp)
				{
					out += string("gini = ");
				}
				else if (env.criterion == env.entropy_imp)
				{
					out += string("entropy = ");
				}
				else
				{
					out += string("impurity = ");
				}
				if (nodeList[i]->total_points != 0)
				{
					out += dtos(std::round((nodeList[i]->impurity / nodeList[i]->total_points) * 1000.0) / 1000.0);
				}
				else
				{
					out += dtos(0);
				}
			}
			out += "<br/>" + string("values = [");
			for (int j = 0; j < n_classes; j++)
			{
				out += itos(nodeList[i]->Ni[j]);
				if (j != n_classes - 1)
				{
					out += ", ";
				}
			}
			out += "]";
		}
		else
		{
			if (env.criterion == Environment::mse)
			{
				out += "<br/>" + string("n = ") + itos(nodeList[i]->total_points);
				out += "<br/>" + string("mse = ") + dtos(std::round((nodeList[i]->mse) * 1000.0) / 1000.0);
			}
			else
			{
				out += "<br/>" + string("n = ") + itos(nodeList[i]->total_points);
				out += "<br/>" + string("mae = ") + dtos(std::round((nodeList[i]->mae) * 1000.0) / 1000.0);
			}
		}
		if (adjList[i].size() == 0)
		{
			if (!regression)
			{
				out += "<br/>" + string("class: ") + y_to_name[nodeList[i]->maxVote];
			}
			else
			{
				if (env.criterion == Environment::mse)
				{
					out += "<br/>" + string("value: ") + dtos(std::round(nodeList[i]->mean * 1000.) / 1000.);
				}
				else
				{
					out += "<br/>" + string("value: ") + dtos(std::round(nodeList[i]->median * 1000.) / 1000.);
				}
			}
		}
		if (!on_edge)
		{
			if (adjList[i].size() > 0 && adjList[i][0]->isInTree)
			{
				if (adjList[i][0]->cond.op == Condition::geq)
				{
					adjList[i][0]->cond.invert();
					adjList[i][1]->cond.invert();
				}
				out += "<br/>" + adjList[i][0]->cond.toString() + "<br/>";
			}
		}
		out += ">] ;\n";
		it2++;
	}
	if (!on_edge)
	{
		map<int, vector<Node *>>::iterator it = adjList.begin();
		while (it != adjList.end())
		{
			int i = it->first;
			bool flag = true;
			for (int j = 0; j < adjList[i].size(); j++)
			{
				if (nodeList[i]->isInTree && adjList[i][j]->isInTree)
				{
					out += itos(i) + " -> " + itos(adjList[i][j]->id);
					if (i == root->id)
					{
						if (adjList[i][j]->cond.op == adjList[i][0]->cond.op)
						{
							out += " [labeldistance=2.2, labelangle=45, headlabel=\"True  \"] ;\n";
						}
						else
						{
							out += " [labeldistance=2.2, labelangle=-45, headlabel=\"   False\"] ;\n";
						}
					}
					else
					{
						out += " ;\n";
					}
				}
			}
			it++;
		}
	}
	else
	{
		map<int, vector<Node *>>::iterator it = adjList.begin();
		while (it != adjList.end())
		{
			int i = it->first;
			if (adjList[i].size() > 0 && (adjList[i][0]->cond.op == Condition::geq))
			{
				for (int j = 0; j < adjList[i].size(); j++)
				{
					adjList[i][j]->cond.invert();
				}
			}
			for (int j = 0; j < adjList[i].size(); j++)
			{
				if (nodeList[i]->isInTree && adjList[i][j]->isInTree)
				{
					out += itos(i) + " -> " + itos(adjList[i][j]->id);
					out += " [label=<" + adjList[i][j]->cond.toString();
					out += ">] ;\n";
				}
			}
			it++;
		}
	}
	out += "}";
	return out;
}
int DecisionTree::nodeCount()
{
	int count = 0;
	map<int, Node *>::iterator it = nodeList.begin();
	while (it != nodeList.end())
	{
		int i = it->first;
		if (nodeList[i]->isInTree)
		{
			count++;
		}
		it++;
	}
	return count;
}
int DecisionTree::leafCount()
{
	int count = 0;
	map<int, Node *>::iterator it = nodeList.begin();
	while (it != nodeList.end())
	{
		int i = it->first;
		if (nodeList[i]->isInTree && nodeList[i]->isLeaf)
		{
			count++;
		}
		it++;
	}
	return count;
}
int DecisionTree::getDepth()
{
	int depth = 0;
	map<int, Node *>::iterator it = nodeList.begin();
	while (it != nodeList.end())
	{
		int i = it->first;
		if (nodeList[i]->isInTree && nodeList[i]->depth > depth)
		{
			depth = nodeList[i]->depth;
		}
		it++;
	}
	return depth;
}
string replaceAll(string str, string from, string to)
{
	string res;
	size_t start = 0;
	size_t to_len = to.length();
	size_t from_len = from.length();
	size_t pos = str.find(from, start);
	while (pos != string::npos)
	{
		res += str.substr(start, pos - start);
		res += to;
		start = pos + from_len;
		pos = str.find(from, start);
	}
	res += str.substr(start, pos);
	return res;
}
string DecisionTree::save()
{
	std::ostringstream f;
	std::map<int, int> parent;
	parent[root->id] = -1;
	for (std::map<int, Node *>::iterator it = nodeList.begin(); it != nodeList.end(); it++)
	{
		int i = it->first;
		if (!nodeList[i]->isLeaf)
		{
			for (int j = 0; j < adjList[i].size(); j++)
			{
				parent[adjList[i][j]->id] = i;
			}
		}
	}
	f << "Identifier: " << id << "\n";
	if (regression)
	{
		f << "Regression\n";
	}
	else
	{
		f << "Classification: " << n_classes << "\n";
	}
	f << "Crit: " << env.criterion << "\n";
	f << "Attributes: " << total_features << "\n";
	f << "AttributeNames:";
	for (int i = 0; i < names.size(); i++)
	{
		string name_tmp = replaceAll(names[i], " ", "[blank]");
		f << " '" << name_tmp << "'";
	}
	f << "\n";
	f << "Nominal:";
	for (int i = 0; i < nominal_indices.size(); i++)
	{
		f << " " << nominal_indices[i];
	}
	f << "\n";
	f << "NominalDomain:";
	for (int i = 0; i < nominal_domain.size(); i++)
	{
		f << " " << nominal_domain[i];
	}
	f << "\n";
	for (int j = 0; j < is_nominal_val_in_dataset.size(); j++)
	{
		f << "NomInDataset:";
		for (int v = 0; v < is_nominal_val_in_dataset[j].size(); v++)
		{
			f << " " << is_nominal_val_in_dataset[j][v];
		}
		f << "\n";
	}
	/*
		Namen der Kategorien der nominalen Features
	*/
	for (map<int, map<int, string>>::iterator it1 = cat_to_name.begin(); it1 != cat_to_name.end(); it1++)
	{
		std::map<int, string> tmp = it1->second;
		f << "CatNames: " << it1->first;
		for (map<int, string>::iterator it2 = tmp.begin(); it2 != tmp.end(); it2++)
		{
			f << " " << it2->first << " " << it2->second << " ";
		}
		f << "\n";
	}
	f << "ClassNames:";
	for (std::map<int, string>::iterator it = y_to_name.begin(); it != y_to_name.end(); it++)
	{
		f << " " << it->first << " " << replaceAll(it->second, " ", "[blank]");
	}
	f << "\n";
	f << "RootId: " << root->id << "\n";
	for (std::map<int, Node *>::iterator it = nodeList.begin(); it != nodeList.end(); it++)
	{
		Node *node = nodeList[it->first];
		f << "Node: \n";
		f << "id: " << node->id << "\n";
		f << "depth: " << node->depth << "\n";
		f << "parent: " << parent[node->id] << "\n";
		f << "n: " << node->total_points << "\n";
		f << "isLeaf: " << node->isLeaf << "\n";
		f << "isInTree: " << node->isInTree << "\n";
		if (regression)
		{
			f << "median: " << node->median << "\n";
			f << "mean: " << node->mean << "\n";
			f << "mae: " << node->mae << "\n";
			f << "mse: " << node->mse << "\n";
		}
		else
		{
			f << "pred: " << node->maxVote << "\n";
			f << "distr:";
			for (int j = 0; j < n_classes; j++)
			{
				f << " " << node->Ni[j];
			}
			f << "\n";
			f << "imp: " << std::hexfloat << node->impurity << std::defaultfloat << "\n";
		}
		if (node->id != root->id)
		{
			f << "rule: ";
			if (node->cond.isCross)
			{
				f << "cross: ";
				f << node->cond.ind1 << " ";
				f << node->cond.ind2 << " ";
				f << node->cond.x1 << " ";
				f << node->cond.x2 << " ";
				f << node->cond.op << " ";
				f << node->cond.op2 << " ";
				f << "\n";
			}
			else if (node->cond.isNominal)
			{
				f << "nom: ";
				f << node->cond.nom_feature;
				for (int v = 0; v < node->cond.in_partition.size(); v++)
				{
					f << " " << node->cond.in_partition[v];
				}
				f << "\n";
				f << node->cond.is_nominal_val_in_dataset.size() << "\n";
				for (int j = 0; j < node->cond.is_nominal_val_in_dataset.size(); j++)
				{
					for (int k = 0; k < node->cond.is_nominal_val_in_dataset[j].size(); k++)
					{
						f << " " << node->cond.is_nominal_val_in_dataset[j][k];
					}
					f << "\n";
				}
			}
			else
			{
				f << "obl:";
				for (int i = 0; i < node->cond.a.size(); i++)
				{
					f << " " << std::hexfloat << node->cond.a[i] << std::defaultfloat;
				}
				f << " " << node->cond.op;
				f << " " << std::hexfloat << node->cond.b << std::defaultfloat;
			}
			f << "\n";
		}
		f << "End Node\n";
	}
	f << "Leaves:";
	for (int i = 0; i < leaves.size(); i++)
	{
		f << " " << leaves[i]->id;
	}
	f << "\n";
	for (std::map<int, vector<Node *>>::iterator it = adjList.begin(); it != adjList.end(); it++)
	{
		f << "AdjList: " << it->first;
		for (int i = 0; i < adjList[it->first].size(); i++)
		{
			f << " " << adjList[it->first][i]->id;
		}
		f << "\n";
	}
	f<<"End Tree\n";
	return f.str();
}
void DecisionTree::save(string filename, int mode)
{
	bool exists = std::filesystem::exists(filename);
	bool is_direct = std::filesystem::is_directory(filename);
	if (exists && mode==0)
	{
		cout << "File already exists.\n";
		return;
	}
	if (exists && is_direct)
	{
		cout <<filename<<" is a directory.\n";
		return;
	}
	ofstream f;
	if(mode==2){
		f.open(filename, ios::out | ios::app);
	}
	else{
		f.open(filename, ios::out | ios::trunc);
	}
	if (f.is_open())
	{
		f << save();
		f.close();
	}
	else
	{
		cout << "Could not open file.\n";
	}
}
void DecisionTree::save(string filename,bool overwrite)
{
	if(overwrite){
		save(filename,1);	
	}
	else{
		save(filename,0);
	}
}
vector<string> split(string str, string delim)
{
	vector<string> res;
	size_t start = 0;
	size_t delim_len = delim.length();
	size_t end = str.find(delim, start);
	while (end != string::npos)
	{
		res.push_back(str.substr(start, end - start));
		start = end + delim_len;
		end = str.find(delim, start);
	}
	res.push_back(str.substr(start, str.length() - start));
	return res;
}
void DecisionTree::load(string filename)
{
	ifstream myfile(filename);
	if (myfile.is_open())
	{
		streamload(myfile);
		myfile.close();
	}
}
void DecisionTree::streamload(istream &in)
{
	if (fitted)
	{
		cout << "Decision Tree is already fitted\n";
		return;
	}
	string line;
	p = new ctpl::thread_pool(env.n_threads);
	std::map<int, int> parent;
	int rootId = 0;
	while (getline(in, line) && line.compare("End Tree")!=0)
	{
		vector<string> tokens = split(line, " ");
		if (tokens[0].compare("Identifier:") == 0)
		{
			id = atoi(tokens[1].c_str());
		}
		else if (tokens[0].compare("Classification:") == 0)
		{
			env.regression = false;
			regression = false;
			n_classes = atoi(tokens[1].c_str());
		}
		else if (tokens[0].compare("Regression") == 0)
		{
			env.regression = true;
		}
		else if (tokens[0].compare("Crit:") == 0)
		{
			env.criterion = atoi(tokens[1].c_str());
		}
		else if (tokens[0].compare("Attributes:") == 0)
		{
			total_features = atoi(tokens[1].c_str());
		}
		else if (tokens[0].compare("AttributeNames:") == 0)
		{
			names.clear();
			for (int i = 1; i < tokens.size(); i++)
			{
				names.push_back(replaceAll(tokens[i], "[blank]", " "));
			}
		}
		else if (tokens[0].compare("Nominal:") == 0)
		{
			nominal_indices.clear();
			for (int i = 1; i < tokens.size(); i++)
			{
				nominal_indices.push_back(atoi(tokens[i].c_str()));
			}
		}
		else if (tokens[0].compare("NominalDomain:") == 0)
		{
			nominal_domain.clear();
			for (int i = 1; i < tokens.size(); i++)
			{
				nominal_domain.push_back(atoi(tokens[i].c_str()));
			}
		}
		else if (tokens[0].compare("NomInDataset:") == 0)
		{
			is_nominal_val_in_dataset.clear();
			vector<bool> tmp;
			for (int i = 1; i < tokens.size(); i++)
			{
				tmp.push_back(atoi(tokens[i].c_str()));
			}
			is_nominal_val_in_dataset.push_back(tmp);
		}
		else if (tokens[0].compare("CatNames:") == 0)
		{
			int j = atoi(tokens[1].c_str());
			cat_to_name[j] = map<int, string>();
			for (int i = 2; i < tokens.size() - 1; i += 2)
			{
				cat_to_name[j][atoi(tokens[i].c_str())] = replaceAll(tokens[i + 1], "[blank]", " ");
			}
		}
		else if (tokens[0].compare("ClassNames:") == 0)
		{
			for (int i = 1; i < tokens.size() - 1; i += 2)
			{
				y_to_name[atoi(tokens[i].c_str())] = replaceAll(tokens[i + 1], "[blank]", " ");
			}
		}
		else if (tokens[0].compare("RootId:") == 0)
		{
			rootId = atoi(tokens[1].c_str());
		}
		else if (tokens[0].compare("Node:") == 0)
		{
			getline(in, line);
			Node *node = new Node();
			node->hasData = false;
			node->possesses_data = false;
			while (line.compare("End Node") != 0)
			{
				tokens = split(line, " ");
				if (tokens[0].compare("id:") == 0)
				{
					node->id = atoi(tokens[1].c_str());
					nodeList[node->id] = node;
					if (node->id == rootId)
					{
						root = node;
					}
				}
				else if (tokens[0].compare("depth:") == 0)
				{
					node->depth = atoi(tokens[1].c_str());
					if (depth < node->depth)
					{
						depth = node->depth;
					}
				}
				else if (tokens[0].compare("parent:") == 0)
				{
					parent[node->id] = atoi(tokens[1].c_str());
				}
				else if (tokens[0].compare("n:") == 0)
				{
					node->total_points = atoi(tokens[1].c_str());
				}
				else if (tokens[0].compare("isLeaf:") == 0)
				{
					node->isLeaf = atoi(tokens[1].c_str());
				}
				else if (tokens[0].compare("isInTree:") == 0)
				{
					node->isInTree = atoi(tokens[1].c_str());
				}
				else if (tokens[0].compare("median:") == 0)
				{
					node->median = atoi(tokens[1].c_str());
				}
				else if (tokens[0].compare("mean:") == 0)
				{
					node->mean = atoi(tokens[1].c_str());
				}
				else if (tokens[0].compare("mae:") == 0)
				{
					node->mae = atoi(tokens[1].c_str());
				}
				else if (tokens[0].compare("mse:") == 0)
				{
					node->mse = atoi(tokens[1].c_str());
				}
				else if (tokens[0].compare("pred:") == 0)
				{
					node->maxVote = atoi(tokens[1].c_str());
				}
				else if (tokens[0].compare("distr:") == 0)
				{
					node->Ni = new double[tokens.size() - 1];
					for (int k = 1; k < tokens.size(); k++)
					{
						node->Ni[k - 1] = atof(tokens[k].c_str());
					}
				}
				else if (tokens[0].compare("imp:") == 0)
				{
					node->impurity = atof(tokens[1].c_str());
				}
				else if (tokens[0].compare("rule:") == 0)
				{
					if (tokens[1].compare("obl:") == 0)
					{
						vector<double> a;
						for (int j = 2; j < tokens.size() - 2; j++)
						{
							a.push_back(atof(tokens[j].c_str()));
						}
						int op = atoi(tokens[tokens.size() - 2].c_str());
						double b = atof(tokens[tokens.size() - 1].c_str());
						node->cond = Condition(a, b, op);
					}
					else if (tokens[1].compare("cross:") == 0)
					{
						int ind1 = atoi(tokens[2].c_str());
						int ind2 = atoi(tokens[3].c_str());
						double x1 = atof(tokens[4].c_str());
						double x2 = atof(tokens[5].c_str());
						int op1 = atoi(tokens[6].c_str());
						int op2 = atoi(tokens[7].c_str());
						node->cond = Condition(ind1, ind2, x1, x2, op1, op2);
					}
					else if (tokens[1].compare("nom:") == 0)
					{
						int nom_feature = atoi(tokens[2].c_str());
						vector<bool> in_part;
						for (int k = 3; k < tokens.size(); k++)
						{
							in_part.push_back(atoi(tokens[k].c_str()));
						}
						vector<vector<bool>> in_dataset;
						getline(in, line);
						int lines = atoi(line.c_str());
						for (int j = 0; j < lines; j++)
						{
							getline(in, line);
							tokens = split(line, " ");
							vector<bool> tmp;
							for (int k = 0; k < tokens.size(); k++)
							{
								tmp.push_back(atoi(tokens[k].c_str()));
							}
							in_dataset.push_back(tmp);
						}
						node->cond = Condition(nom_feature, in_part, in_dataset);
					}
				}
				getline(in, line);
			}
		}
		else if (tokens[0].compare("Leaves:") == 0)
		{
			for (int k = 1; k < tokens.size(); k++)
			{
				int id = atoi(tokens[k].c_str());
				leaves.push_back(nodeList[id]);
			}
		}
		else if (tokens[0].compare("AdjList:") == 0)
		{
			int id = atoi(tokens[1].c_str());
			adjList[id] = vector<Node *>();
			for (int k = 2; k < tokens.size(); k++)
			{
				int id2 = atoi(tokens[k].c_str());
				adjList[id].push_back(nodeList[id2]);
			}
		}
	}
	vector<bool> is_nominal(total_features, false);
	for (int k = 0; k < nominal_indices.size(); k++)
	{
		is_nominal[nominal_indices[k]] = true;
	}
	for (int k = 0; k < names.size(); k++)
	{
		if (is_nominal[k])
		{
			nominalNames.push_back(names[k]);
		}
		else
		{
			numericNames.push_back(names[k]);
		}
	}
	total_numeric_features = total_features - nominal_indices.size();
	total_nominal_features = nominal_indices.size();
	fitted = true;
	nodes_have_data = false;
	root_possesses_data = false;
	map<int, Node *>::iterator it = nodeList.begin();
	while (it != nodeList.end())
	{
		int i = it->first;
		nodeList[i]->cond.names = numericNames;
		nodeList[i]->cond.nominalNames = nominalNames;
		nodeList[i]->cond.cat_to_name = cat_to_name;
		nodeList[i]->maxVote_str = y_to_name[nodeList[i]->maxVote];
		it++;
	}
}
