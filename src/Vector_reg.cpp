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

#include "Vector_reg.h"
#include "Environment.h"
using namespace std;
void Vector_reg::use_mae()
{
    keep_vec = true;
    median = true;
};
Vector_reg::Vector_reg()
{
    length = 0;
};
Vector_reg::Vector_reg(const Vector_reg &a)
{
    this->sum_t = a.sum_t;
    this->sum2_t = a.sum2_t;
    this->length = a.length;
    this->sq_error = a.sq_error;
    this->a_error = a.a_error;
    this->is_sorted = a.is_sorted;
    this->keep_vec = a.keep_vec;
    this->median = a.median;
    if (this->keep_vec)
    {
        this->vec = a.vec;
    }
};
Vector_reg::Vector_reg(vector<double> a)
{
    if (keep_vec)
    {
        this->vec = a;
    }
    for (int i = 0; i < a.size(); i++)
    {
        sum_t += a[i];
        sum2_t += a[i] * a[i];
    }
    length = a.size();
    sq_error = -1;
    a_error = -1;
    is_sorted = false;
};
int Vector_reg::size()
{
    return length;
};
double Vector_reg::sum()
{
    return sum_t;
};
double Vector_reg::sum2()
{
    return sum2_t;
};
void Vector_reg::push_back(double b)
{
    sum_t += b;
    sum2_t += b * b;
    if (keep_vec)
    {
        vec.push_back(b);
    }
    length += 1;
    sq_error = -1;
    a_error = -1;
    is_sorted = false;
}
void Vector_reg::merge(Vector_reg b)
{
    sum_t += b.sum();
    sum2_t += b.sum2();
    if (keep_vec)
    {
        vec.insert(vec.end(), b.vec.begin(), b.vec.end());
    }
    length += b.size();
    sq_error = -1;
    a_error = -1;
    is_sorted = false;
};
double Vector_reg::squaredError()
{
    if (sq_error == -1)
    {
        if (length > 0)
        {
            sq_error = sum2_t - 1 * sum_t * (sum_t / length);
        }
        else
        {
            sq_error = 0;
        }
        return sq_error;
    }
    else
    {
        return sq_error;
    }
};
double Vector_reg::absoluteError()
{
    if (!keep_vec)
    {
        cout << "ERROR: Must store values for use with absolute errors.\n";
        exit(0);
    }
    if (a_error == -1)
    {
        a_error = 0;
        if (length > 0)
        {
            double pred = sum_t / (double)length;
            if (median)
            {
                if (is_sorted)
                {
                    if (length % 2 == 0)
                    {
                        pred = (vec[length / 2 - 1] + vec[length / 2]) / 2;
                    }
                    else
                    {
                        pred = vec[length / 2];
                    }
                }
                else
                {
                    vector<double> vec2 = vec;
                    std::sort(vec2.begin(), vec2.end());
                    if (length % 2 == 0)
                    {
                        pred = (vec2[length / 2 - 1] + vec2[length / 2]) / 2;
                    }
                    else
                    {
                        pred = vec2[length / 2];
                    }
                }
            }
            for (int i = 0; i < length; i++)
            {
                a_error += abs(vec[i] - pred);
            }
        }
        else
        {
            a_error = 0;
        }
        return a_error;
    }
    else
    {
        return a_error;
    }
};
double Vector_reg::squaredError(Vector_reg b)
{
    int length2 = length + b.size();
    double sum_t2 = sum_t + b.sum();
    double sum2_t2 = sum2_t + b.sum2();
    return sum2_t2 - 1 * sum_t2 * (sum_t2 / length2);
};
void Vector_reg::copy(Vector_reg *b)
{
    b->sum_t = sum_t;
    b->sum2_t = sum2_t;
    b->length = length;
    b->sq_error = sq_error;
    b->a_error = a_error;
    if (keep_vec)
    {
        b->vec = vec;
    }
    b->is_sorted = is_sorted;
    b->keep_vec = keep_vec;
}
Vector_reg Vector_reg::add(Vector_reg &b)
{
    Vector_reg tmp;
    tmp.sum_t = sum_t + b.sum_t;
    tmp.sum2_t = sum2_t + b.sum2_t;
    tmp.length = length + b.length;
    tmp.sq_error = -1;
    tmp.a_error = -1;
    if (keep_vec)
    {
        tmp.vec.reserve(tmp.length);
        tmp.vec.insert(tmp.vec.end(), vec.begin(), vec.end());
        tmp.vec.insert(tmp.vec.end(), b.vec.begin(), b.vec.end());
    }
    tmp.keep_vec = keep_vec;
    tmp.is_sorted = false;
    return tmp;
};
void Vector_reg::add(Vector_reg *b)
{
    sum_t = sum_t + b->sum_t;
    sum2_t = sum2_t + b->sum2_t;
    length = length + b->length;
    sq_error = -1;
    a_error = -1;
    if (keep_vec)
    {
        vec.insert(vec.end(), b->vec.begin(), b->vec.end());
    }
    keep_vec = keep_vec;
    is_sorted = false;
};
double Vector_reg::error(int criterion)
{
    if (criterion == Environment::mse)
    {
        return squaredError();
    }
    else
    {
        return absoluteError();
    }
}
void Vector_reg::clear()
{
    sum_t = 0;
    sum2_t = 0;
    if (keep_vec)
    {
        vec.clear();
    }
    length = 0;
    sq_error = -1;
    a_error = -1;
    is_sorted = false;
}
void Vector_reg::reserve(int n)
{
    if (keep_vec)
    {
        vec.reserve(n);
    }
}
double Vector_reg::pop_back()
{
    if (!keep_vec)
    {
        cout << "ERROR: Requires to store values.\n";
        exit(0);
    }
    if (length == 0)
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
    double v = vec.back();
    sum_t -= v;
    sum2_t -= v * v;
    length -= 1;
    sq_error = -1;
    a_error = -1;
    vec.pop_back();
    return v;
}
bool Vector_reg::remove(double b)
{
    if (!keep_vec)
    {
        sum_t -= b;
        sum2_t -= b * b;
        length -= 1;
        sq_error = -1;
        a_error = -1;
        return true;
    }
    if (!is_sorted)
    {
        std::vector<double>::iterator it = find(vec.begin(), vec.end(), b);
        if (it == vec.end() || *it != b)
        {
            return false;
        }
        else
        {
            sum_t -= b;
            sum2_t -= b * b;
            length -= 1;
            sq_error = -1;
            a_error = -1;
            vec.erase(it);
            return true;
        }
    }
    else
    {
        std::vector<double>::iterator it = lower_bound(vec.begin(), vec.end(), b);
        if (it != vec.end() && *it == b)
        {
            sum_t -= b;
            sum2_t -= b * b;
            length -= 1;
            sq_error = -1;
            a_error = -1;
            vec.erase(it);
            return true;
        }
        else
        {
            return false;
        }
    }
}
void Vector_reg::sort()
{
    if (!keep_vec)
    {
        is_sorted = false;
        return;
    }
    if (length > 0)
    {
        std::sort(vec.begin(), vec.end());
    }
    is_sorted = true;
}
void Vector_reg::insert_inplace(double b)
{
    if (!keep_vec)
    {
        push_back(b);
        is_sorted = false;
        return;
    }
    if (length == 0)
    {
        push_back(b);
        is_sorted = true;
        return;
    }
    else
    {
        if (!is_sorted)
        {
            push_back(b);
            sort();
            return;
        }
        vec.insert(std::upper_bound(vec.begin(), vec.end(), b), b);
        sum_t += b;
        sum2_t += b * b;
        length += 1;
        sq_error = -1;
        a_error = -1;
    }
}
