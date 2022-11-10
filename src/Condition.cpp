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

#include "Condition.h"
#include "General.h"
#include <math.h>
Condition::Condition(){};
Condition::Condition(vector<double> a, double b, int op)
{
	isNominal = false;
	this->a = a;
	this->b = b;
	this->op = op;
}
Condition::Condition(int ind1, int ind2, double x1, double x2, int op1, int op2)
{
	this->ind1 = ind1;
	this->ind2 = ind2;
	this->x1 = x1;
	this->x2 = x2;
	this->op = op1;
	this->op2 = op2;
	isNominal = false;
	isCross = true;
}
Condition::Condition(int nom_feature, vector<bool> in_partition, vector<vector<bool>> is_nominal_val_in_dataset)
{
	isNominal = true;
	this->in_partition = in_partition;
	this->nom_feature = nom_feature;
	this->is_nominal_val_in_dataset = is_nominal_val_in_dataset;
}
bool Condition::evaluate(vector<double> &x, vector<double> &x_nominal)
{
	if (!isNominal)
	{
		if (!isCross)
		{
			double sum = -b;
			for (int j = 0; j < x.size(); j++)
			{
				sum += a[j] * x[j];
			}
			if (op == geq)
			{
				return (sum >= 0);
			}
			else if (op == l)
			{
				return (sum < 0);
			}
			else if (op == leq)
			{
				return (sum <= 0);
			}
			else if (op == g)
			{
				return (sum > 0);
			}
			else
			{
				return true;
			}
		}
		else
		{
			if (op == geq && op2 == geq)
			{
				return x[ind1] >= x1 && x[ind2] >= x2;
			}
			else if (op == geq && op2 == g)
			{
				return x[ind1] >= x1 && x[ind2] > x2;
			}
			else if (op == geq && op2 == leq)
			{
				return x[ind1] >= x1 && x[ind2] <= x2;
			}
			else if (op == geq && op2 == l)
			{
				return x[ind1] >= x1 && x[ind2] < x2;
			}
			else if (op == g && op2 == geq)
			{
				return x[ind1] > x1 && x[ind2] >= x2;
			}
			else if (op == g && op2 == g)
			{
				return x[ind1] > x1 && x[ind2] > x2;
			}
			else if (op == g && op2 == leq)
			{
				return x[ind1] > x1 && x[ind2] <= x2;
			}
			else if (op == g && op2 == l)
			{
				return x[ind1] > x1 && x[ind2] < x2;
			}
			else if (op == leq && op2 == geq)
			{
				return x[ind1] <= x1 && x[ind2] >= x2;
			}
			else if (op == leq && op2 == g)
			{
				return x[ind1] <= x1 && x[ind2] > x2;
			}
			else if (op == leq && op2 == leq)
			{
				return x[ind1] <= x1 && x[ind2] <= x2;
			}
			else if (op == leq && op2 == l)
			{
				return x[ind1] <= x1 && x[ind2] < x2;
			}
			else if (op == l && op2 == geq)
			{
				return x[ind1] < x1 && x[ind2] >= x2;
			}
			else if (op == l && op2 == g)
			{
				return x[ind1] < x1 && x[ind2] > x2;
			}
			else if (op == l && op2 == leq)
			{
				return x[ind1] < x1 && x[ind2] <= x2;
			}
			else
			{
				return x[ind1] < x1 && x[ind2] < x2;
			}
		}
	}
	else
	{
		int val = (int)x_nominal[nom_feature];
		if (!is_nominal_val_in_dataset[nom_feature][val])
		{
			cout << "ERROR: Value " << x_nominal[nom_feature] << " not present in training data\n";
		}
		return in_partition[x_nominal[nom_feature]];
	}
}
string Condition::toString()
{
	if (!isNominal)
	{
		if (!isCross)
		{
			string out = "";
			bool first = true;
			bool is_oblique = false;
			bool is_bivariate = false;
			int count = 0;
			for (int j = 0; j < a.size(); j++)
			{
				if (a[j] != 0)
				{
					count++;
				}
				if (count >= 2)
				{
					is_oblique = true;
				}
			}
			if (count == 2)
			{
				is_bivariate = true;
			}
			double denom = 1;
			bool scale = true;
			double denom2 = std::numeric_limits<double>::max();
			if (scale && is_oblique)
			{
				denom = b * b;
				for (int j = 0; j < a.size(); j++)
				{
					denom += a[j] * a[j];
					if (abs(a[j]) < denom2 && abs(a[j]) > 1e-3)
					{
						denom2 = abs(a[j]);
					}
				}
				denom = sqrt(denom);
				if (is_bivariate)
				{
					denom = denom2;
				}
			}
			for (int j = 0; j < a.size(); j++)
			{
				if (a[j] != 0)
				{
					double val = a[j] / denom;
					if (a[j] < 0)
					{
						val = -a[j] / denom;
					}
					if (!first && a[j] > 0)
					{
						out += "+";
					}
					if (!first && a[j] < 0)
					{
						out += "-";
					}
					if (first && a[j] < 0)
					{
						out += "-";
					}
					double nearest = round(val * 100.0) / 100.0;
					if (nearest == 1.00)
					{
						if (names.size() > 0)
						{
							out += names[j];
						}
						else
						{
							out += "x" + std::string("<sub>") + itos(j) + "</sub>";
						}
					}
					else
					{
						if (names.size() > 0)
						{
							out += dtos(nearest) + "&#183;" + names[j];
						}
						else
						{
							out += dtos(nearest) + "&#183;" + "x" + "<sub>" + itos(j) + "</sub>";
						}
					}
					first = false;
				}
			}
			double nearest = round((b / denom) * 100.0) / 100.0;
			if (nearest == 0)
			{
				nearest = 0;
			}
			if (op == leq)
			{
				out += " &#8804; " + dtos(nearest);
			}
			if (op == geq)
			{
				out += " &#8805; " + dtos(nearest);
			}
			if (op == l)
			{
				out += " &#60; " + dtos(nearest);
			}
			if (op == g)
			{
				out += " &#62; " + dtos(nearest);
			}
			return out;
		}
		else
		{
			string out = "";
			if (names.size() > 0)
			{
				out += names[ind1];
			}
			else
			{
				out += "x" + std::string("<sub>") + itos(ind1) + "</sub>";
			}
			double nearest = round(x1 * 1000.0) / 1000.0;
			if (op == leq)
			{
				out += " &#8804; " + dtos(nearest);
			}
			if (op == geq)
			{
				out += " &#8805; " + dtos(nearest);
			}
			if (op == l)
			{
				out += " &#60; " + dtos(nearest);
			}
			if (op == g)
			{
				out += " &#62; " + dtos(nearest);
			}
			out += " &#8743; ";
			if (names.size() > 0)
			{
				out += names[ind2];
			}
			else
			{
				out += "x" + std::string("<sub>") + itos(ind2) + "</sub>";
			}
			nearest = round(x2 * 1000.0) / 1000.0;
			if (op2 == leq)
			{
				out += " &#8804; " + dtos(nearest);
			}
			if (op2 == geq)
			{
				out += " &#8805; " + dtos(nearest);
			}
			if (op2 == l)
			{
				out += " &#60; " + dtos(nearest);
			}
			if (op2 == g)
			{
				out += " &#62; " + dtos(nearest);
			}
			return out;
		}
	}
	else
	{
		string out = nominalNames[nom_feature] + "&#8712;{";
		for (int k = 0; k < in_partition.size(); k++)
		{
			if (in_partition[k] && is_nominal_val_in_dataset[nom_feature][k])
			{
				if (cat_to_name.find(nom_feature) == cat_to_name.end())
				{
					out += itos(k);
				}
				else
				{
					if (cat_to_name[nom_feature].find(k) == cat_to_name[nom_feature].end())
					{
						out += itos(k);
					}
					else
					{
						out += cat_to_name[nom_feature][k];
					}
				}
				out += ",";
			}
		}
		out = out.substr(0, out.length() - 1);
		out += "}";
		return out;
	}
}
void Condition::invert()
{
	if (isCross || isNominal)
	{
		return;
	}
	if (op == geq)
	{
		op = leq;
	}
	else if (op == l)
	{
		op = g;
	}
	else if (op == leq)
	{
		op = geq;
	}
	else if (op == g)
	{
		op = l;
	}
	for (int j = 0; j < a.size(); j++)
	{
		a[j] = -a[j];
	}
	b = -b;
}
Condition Condition::copy()
{
	Condition cond;
	for (int i = 0; i < a.size(); i++)
	{
		cond.a.push_back(a[i]);
	}
	cond.b = b;
	cond.condition = "" + condition;
	cond.op = op;
	for (int i = 0; i < names.size(); i++)
	{
		cond.names.push_back(names[i]);
	}
	for (int i = 0; i < nominalNames.size(); i++)
	{
		cond.nominalNames.push_back(nominalNames[i]);
	}
	cond.isNominal = isNominal;
	cond.nom_feature = nom_feature;
	for (int i = 0; i < in_partition.size(); i++)
	{
		cond.in_partition.push_back(in_partition[i]);
	}
	for (int i = 0; i < is_nominal_val_in_dataset.size(); i++)
	{
		cond.is_nominal_val_in_dataset.push_back(vector<bool>());
		for (int j = 0; j < is_nominal_val_in_dataset[i].size(); j++)
		{
			cond.is_nominal_val_in_dataset[i].push_back(is_nominal_val_in_dataset[i][j]);
		}
	}
	cond.isCross = isCross;
	cond.ind1 = ind1;
	cond.ind2 = ind2;
	cond.x1 = x1;
	cond.x2 = x2;
	return cond;
}
