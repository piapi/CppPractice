#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <math.h>
#include <vector>
#include <deque>
#include <iomanip>
#include <queue>
#include <fstream>
using namespace std;
struct Node
{
    vector<double> point;
    Node *left;
    Node *right;

    Node(const vector<double> &point) : point(point), left(nullptr), right(nullptr) {}
};

class KDTree
{
public:
    void build(const vector<vector<double>> &points)
    {
        root = buildTree(points, 0);
    }

    vector<vector<double>> searchKNearestNeighbors(const vector<double> &target, int k)
    {
        priority_queue<pair<double, vector<double>>> pq;
        searchKNearestNeighbors(root, target, k, pq, 0);

        vector<vector<double>> result;
        while (!pq.empty())
        {
            result.push_back(pq.top().second);
            pq.pop();
        }

        return result;
    }

private:
    Node *root;
    // 找到当前维度的中间值的索引
    int FindMidIndex(const vector<vector<double>> &points,int aixs){
        vector<int> temp(points.size(),0);
        for(int i=0;i<points.size();i++){
            temp[i]=points[i][aixs];
        }
        sort(temp.begin(),temp.end());
        int mid=temp[temp.size()/2];
        for(int i=0;i<points.size();i++){
            if(points[i][aixs]==mid){
                return i;
            }
        }
        return 0;
    }
    
    Node *buildTree(const vector<vector<double>> &points, int depth)
    {
        // a样本数量
        int sample_num = points.size();
        if (sample_num == 0)
        {
            return nullptr;
        }
        if (sample_num == 1)
        {
            return new Node(points[0]);
        }
        // 数据维度
        int dim = points[0].size();
        // 当前维度
        int curDim = depth % dim;
        int MIndex=FindMidIndex(points,curDim);
        cout<<depth<<" "<<MIndex<<endl;
        Node *node = new Node(points[MIndex]);
        vector<vector<double>> leftSon;
        vector<vector<double>> rightSon;
        for (size_t i = 0; i < points.size(); i++)
        {
            if(i==MIndex)continue;
            if (points[i][curDim] < points[MIndex][curDim])
            {
                leftSon.push_back(points[i]);
            }
            else
            {
                rightSon.push_back(points[i]);
            }
        }
        node->left = buildTree(leftSon, depth + 1);
        node->right = buildTree(rightSon, depth + 1);

        return node;
    }

    void
    searchKNearestNeighbors(Node *node, const vector<double> &target, int k, priority_queue<pair<double, vector<double>>> &pq, int depth)
    {
        if (node == nullptr)
            return;

        double distance = calculateDistance(node->point, target);
        pq.push(make_pair(distance, node->point));

        if (pq.size() > k)
        {
            pq.pop();
        }

        int axis = depth % target.size(); // 根据维度选择轴

        if (target[axis] < node->point[axis])
        {
            searchKNearestNeighbors(node->left, target, k, pq, depth + 1);
            if (pq.size() < k || abs(target[axis] - node->point[axis]) < pq.top().first)
            {
                searchKNearestNeighbors(node->right, target, k, pq, depth + 1);
            }
        }
        else
        {
            searchKNearestNeighbors(node->right, target, k, pq, depth + 1);
            if (pq.size() < k || abs(target[axis] - node->point[axis]) < pq.top().first)
            {
                searchKNearestNeighbors(node->left, target, k, pq, depth + 1);
            }
        }
    }

    double calculateDistance(const vector<double> &p1, const vector<double> &p2)
    {
        double sum = 0.0;
        for (int i = 0; i < p1.size(); ++i)
        {
            double diff = p1[i] - p2[i];
            sum += diff * diff;
        }
        return sqrt(sum);
    }
};

int main()
{
    vector<vector<double>> points = {
        {{2.0, 3.0, 4.0}},
        {{5.0, 1.0, 9.0}},
        {{7.0, 2.0, 6.0}},
        {{2.3, 3.0, 4.0}},
        {{2.9, 3.0, 4.0}}};

    KDTree kdTree;
    kdTree.build(points);
    vector<double> target({2.1, 3.5, 4.0});
    int k = 3;
    vector<vector<double>> nearestNeighbors = kdTree.searchKNearestNeighbors(target, k);

    cout << "Nearest " << k << " neighbors to the target point (";
    for (int i = 0; i < target.size(); ++i)
    {
        cout << target[i];
        if (i < target.size() - 1)
            cout << ", ";
    }
    cout << "):" << endl;
    for (const auto &point : nearestNeighbors)
    {
        cout << "(";
        for (int i = 0; i < point.size(); ++i)
        {
            cout << point[i];
            if (i < point.size() - 1)
                cout << ", ";
        }
        cout << ")" << endl;
    }

    return 0;
}