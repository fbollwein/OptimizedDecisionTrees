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
#include <deque>
#include <queue>
#include <algorithm>
#include <set>
template <class T>
class VarQueue
{
public:
    static const int PQUE = 1;
    static const int DQUE = 0;
    int ind = 0;
    int n_ques = 0;
    vector<std::priority_queue<T>> pques;
    vector<std::deque<T>> dques;
    std::multiset<double> bounds;
    std::priority_queue<T> pque;
    std::deque<T> dque;
    int type;
    int length = 0;
    VarQueue(int type);
    void insert(T node);
    void insert_front(T node);
    void shuffle();
    void removeBound(double bound);
    T next();
    bool isEmpty();
    int size();
};
template <class T>
void VarQueue<T>::shuffle()
{
    if (type == PQUE)
    {
    }
    else
    {
        if (n_ques > 1)
        {
            for (int i = 0; i < n_ques; i++)
            {
                std::random_shuffle(dques[i].begin(), dques[i].end());
            }
        }
        else
        {
            std::random_shuffle(dque.begin(), dque.end());
        }
    }
}
template <class T>
VarQueue<T>::VarQueue(int type)
{
    this->type = type;
    if (this->type != PQUE && this->type != DQUE)
    {
        this->type = DQUE;
    }
}
template <class T>
void VarQueue<T>::insert(T node)
{
    if (type == PQUE)
    {
        pque.push(node);
        bounds.insert(node.LB);
    }
    else
    {
        dque.push_back(node);
        bounds.insert(node.LB);
    }
    length++;
}
template <class T>
void VarQueue<T>::insert_front(T node)
{
    if (type == PQUE)
    {
        pque.push(node);
    }
    else
    {
        dque.push_front(node);
    }
    length++;
}
template <class T>
T VarQueue<T>::next()
{
    length--;
    if (type == PQUE)
    {
        T node = pque.top();
        pque.pop();
        return node;
    }
    else
    {
        T node = dque.back();
        dque.pop_back();
        return node;
    }
}
template <class T>
void VarQueue<T>::removeBound(double bound)
{
    if (type == PQUE)
    {
        bounds.erase(bounds.find(bound));
    }
    else
    {
        bounds.erase(bounds.find(bound));
    }
}
template <class T>
bool VarQueue<T>::isEmpty()
{
    if (type == PQUE)
    {
        return pque.empty();
    }
    else
    {
        return dque.empty();
    }
}
template <class T>
int VarQueue<T>::size()
{
    if (type == PQUE)
    {
        return pque.size();
    }
    else
    {
        return dque.size();
    }
}