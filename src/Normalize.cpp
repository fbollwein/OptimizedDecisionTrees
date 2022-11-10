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

#include <iostream>
#include <limits>
#include <math.h>
#include "Normalize.h"
void normalize(int ***data, double **X2, double **X, int total_points, int d, double *offset, double *denominator, int type)
{
  for (int j = 0; j < d - 1; j++)
  {
    offset[j] = 0;
    denominator[j] = 1;
  }
  for (int j = 0; j < d - 1; j++)
  {
    double min = std::numeric_limits<double>::max();
    double max = std::numeric_limits<double>::min();
    double mean = 0;
    double median;
    double uq;
    double lq;
    double std = 0;
    for (int i = 0; i < total_points; i++)
    {
      mean += X[i][j];
      if (X[i][j] > max)
      {
        max = X[i][j];
      }
      if (X[i][j] < min)
      {
        min = X[i][j];
      }
    }
    mean /= total_points;
    for (int i = 0; i < total_points; i++)
    {
      std += (X[i][j] - min) * (X[i][j] - min);
    }
    std = sqrt(std / (total_points - 1));
    denominator[j] = std;
    if (type == 1)
    {
      if (denominator[j] == 0)
      {
        denominator[j] = 1;
      }
      offset[j] = mean;
    }
    else
    {
      double *X1 = X2[j];
      int **dataX1 = data[j];
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
      if (total_points % 2 == 0)
      {
        median = 0.5 * (X1[dataX1[(total_points) / 2 - 1][0]] + X1[dataX1[(total_points) / 2][0]]);
      }
      else
      {
        median = X1[dataX1[(total_points - 1) / 2][0]];
      }
      denominator[j] = (uq - lq);
      if (denominator[j] == 0)
      {
        denominator[j] = 1;
      }
      offset[j] = median;
    }
  }
  for (int j = 0; j < d - 1; j++)
  {
    for (int i = 0; i < total_points; i++)
    {
      X[i][j] = (X[i][j] - offset[j]) / denominator[j];
    }
  }
}
void normalize_reg(double ***data, double **X2, double **X, int total_points, int d, double *offset, double *denominator, int type)
{
  for (int j = 0; j < d - 1; j++)
  {
    offset[j] = 0;
    denominator[j] = 1;
  }
  for (int j = 0; j < d - 1; j++)
  {
    double min = std::numeric_limits<double>::max();
    double max = std::numeric_limits<double>::min();
    double mean = 0;
    double median;
    double uq;
    double lq;
    double std = 0;
    for (int i = 0; i < total_points; i++)
    {
      mean += X[i][j];
      if (X[i][j] > max)
      {
        max = X[i][j];
      }
      if (X[i][j] < min)
      {
        min = X[i][j];
      }
    }
    mean /= total_points;
    for (int i = 0; i < total_points; i++)
    {
      std += (X[i][j] - min) * (X[i][j] - min);
    }
    std = sqrt(std / (total_points - 1));
    denominator[j] = std;
    if (type == 1)
    {
      if (denominator[j] == 0)
      {
        denominator[j] = 1;
      }
      offset[j] = mean;
    }
    else
    {
      double *X1 = X2[j];
      double **dataX1 = data[j];
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
      if (total_points % 2 == 0)
      {
        median = 0.5 * (X1[(int)dataX1[(total_points) / 2 - 1][0]] + X1[(int)dataX1[(total_points) / 2][0]]);
      }
      else
      {
        median = X1[(int)dataX1[(total_points - 1) / 2][0]];
      }
      denominator[j] = (uq - lq);
      if (denominator[j] == 0)
      {
        denominator[j] = 1;
      }
      offset[j] = median;
    }
  }
  for (int j = 0; j < d - 1; j++)
  {
    for (int i = 0; i < total_points; i++)
    {
      X[i][j] = (X[i][j] - offset[j]) / denominator[j];
    }
  }
}
void denormalize(Result *res, int d, double *offset, double *denominator)
{
  for (int j = 0; j < d - 1; j++)
  {
    res->bestB += (res->a[j] * offset[j]) / denominator[j];
    res->a[j] = res->a[j] / denominator[j];
  }
}