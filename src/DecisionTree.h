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
#include <stdlib.h>
#include <queue>
#include <string>
#include <map>
#include "Node.h"
#include "ctpl_stl.h"
using namespace std;
class DecisionTree
{	
public:
	int id;
	map<int, vector<Node *>> adjList;
	vector<Node *> leaves;
	map<int, Node *> nodeList;
	Node *root;
	int depth = 0;
	int n_leaf_nodes = 0;
	int total_features = 0;
	int total_numeric_features = 0;
	int total_nominal_features = 0;
	int total_points = 0;
	int n_classes = 0;
	int *Y;
	double *Y_reg;
	bool regression = false;
	double **X;
	double **Xnominal;
	bool fitted = false;
	bool nodes_have_data = true;
	vector<vector<bool>> is_nominal_val_in_dataset;
	vector<int> nominal_indices;
	vector<int> nominal_domain;
	vector<string> names;
	vector<string> nominalNames;
	vector<string> numericNames;
	map<int, map<int, string>> cat_to_name;
	map<int, string> y_to_name;
	string criterion;
	Environment env;
	ctpl::thread_pool *p;
	bool root_possesses_data = true;
	void setParam(string param, string value);
	void setParam(vector<string> params, vector<string> values);
	void setParam(string param, double value);
	void setParam(vector<string> params, vector<double> values);
	DecisionTree();
	~DecisionTree();
	int evaluate(vector<vector<double>> X, vector<int> Y2);
	double evaluate_reg(vector<vector<double>> X, vector<double> Y2);
	vector<int> predict(vector<vector<double>> X);
	vector<double> predict_reg(vector<vector<double>> X);
	void fit(vector<vector<double>> X, vector<int> Y);
	void fit(vector<vector<double>> X, vector<int> Y, vector<int> nominal);
	void fit_reg(vector<vector<double>> X, vector<double> Y);
	void fit_reg(vector<vector<double>> X, vector<double> Y, vector<int> nominal);
	void setAttributeNames(vector<string> names);
	void setNominalNames(int index, vector<int> values, vector<string> categories);
	void setYNames(vector<int> values, vector<string> categories);
	int nodeCount();
	int leafCount();
	string export_graphviz(bool on_edge);
	void removeSubTree(int id);
	void prune(vector<vector<double>> X, vector<int> Y);
	void prune_reg(vector<vector<double>> X, vector<double> Y);
	vector<int> MCCgetLeaves(int id, std::map<int, bool> &T);
	vector<int> MCCgetSubTree(int id, std::map<int, bool> &T, bool inc);
	int MCCgetSubTreeError(vector<int> ids, std::map<int, bool> &T, std::map<int, double> &error);
	int MCCevaluate(vector<vector<double>> &X, vector<int> &Y2, std::map<int, bool> &T);
	vector<int> MCCgetParents(int id, std::map<int, int> &parent);
	void MCCprune(vector<vector<double>> X, vector<int> Y);
	void MCCprune_reg(vector<vector<double>> &X, vector<double> &Y);
	int MCCevaluate_reg(vector<vector<double>> &X, vector<double> &Y2, std::map<int, bool> &T);
	double getError(int id);
	double getSubTreeError(int id, bool collect, bool eval);
	int getDepth();
	string export_graph();
	string save();
	void save(string filename,int mode);
	void save(string filename,bool overwrite);
	void load(string filename);
	void streamload(istream &in);
};
