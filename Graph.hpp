//
//  Graph.hpp
//  index8(map)
//
//  Created by 孟令凯 on 2021/3/23.
//
#ifndef Graph_hpp
#define Graph_hpp

#include <stdio.h>
#include <iostream>
#include <stdio.h>
#include "string"
#include <vector>
#include <list>
#include <fstream>
#include <queue>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <unordered_set>

using namespace std;
class Graph{
private:
    string str;//图数据
    int n,m,edge_m;
    int max_degree;
    int *max_degree_num;
    int *degree;//每个点邻居数目
    int *pstart,*out_pstart,*in_pstart; //offset of neighbors of nodes 节点相邻点偏移量
    int *edges,*out_edges,*in_edges; //adjacent ids of edges 相邻边
    int *reverse; //the position of reverse edge in edges 反向边在边中的位置
    int *num_c, *num_f;
    float ***core_order_c, ***core_order_f, **neighbor_order_c, **neighbor_order_f;
    
    vector<multimap<float, int, greater<float>>> cores_mp_c;
    vector<multimap<float, int, greater<float>>> nbr_mp_c;
    vector<multimap<float, int, greater<float>>> cores_mp_f;
    vector<multimap<float, int, greater<float>>> nbr_mp_f;
    vector<int> update_degree, update_max_degree_num, update_in_num, update_out_num;
    vector<vector<int>> update_in_edges, update_out_edges;
    vector<vector<float>> nbr_order_old_c, nbr_order_old_f;
public:
    Graph(const char *_dir);
    ~Graph();
    void readGraph();
    void creatIndex();
    void readIndex();
    void readGraph_update();
    void readIndex_update();
    void query(float parameter1,float parameter2, int parameterNum);
    void insert(int u, int v);
    void remove(int u, int v);
    void update_write();
private:
    void getDegree();
    void getTti();
    void getFlowTti();
    void intersection(int *A, int *B, int a1, int a2, int b1, int b2, int *out);
    int binary_search(const int *array, int b, int e, int val) ;
    void getNeighborOrder();
    void getCoreOrder();
    void writeIndex();
    void getAllNei(float parameter1, float parameter2, int v, int *all_nei, int &all_nei_num);
    int BinarySearch(int *a, int value, int n);
    int BinarySearch2(int *a, int value, int b, int e);
    int BinarySearch3(vector<int> a, int value);
    
    void get_add_c_and_f(unordered_map<int,int> &add_c, unordered_map<int,int> &add_f, int u, int v, int &all_c, int &all_f);
    void update_nei1(unordered_map<int,int> &add_c, unordered_map<int,int> &add_f, int all_c, int all_f, int u, int v);//有邻居增加
    void update_nei2(unordered_map<int,int> &add_c, unordered_map<int,int> &add_f, int all_c, int all_f, int u, int v);//无邻居增加
    void update_nei_order_c(int v, int w, float old_sim, float new_sim);
    void update_nei_order_f(int v, int w, float old_sim, float new_sim);
    void update_core_c(int u, float low, float high);//flog = 0:nei_order变长
    void update_core_f(int u, float low, float high);
    
    void delete_update_nei1(unordered_map<int,int> &add_c, unordered_map<int,int> &add_f, int all_c, int all_f, int u, int v);//有邻居减少
    void delete_update_nei2(unordered_map<int,int> &add_c, unordered_map<int,int> &add_f, int all_c, int all_f, int u, int v);//无邻居减少
};
#endif /* Graph_hpp */

