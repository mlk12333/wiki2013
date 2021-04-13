//
//  main.cpp
//  index8(map)
//
//  Created by 孟令凯 on 2021/3/23.
//

#include <iostream>
#include "string"
#include <fstream>
#include "Graph.hpp"
using namespace std;

int main(int argc, const char * argv[]) {
    if(argc<2) return 1;
    // insert code here...
//    string str = "/Users/milk/test_data/zata/index";
//    argv[0] = &str[0];
    clock_t startTime,endTime;
    
    Graph *graph = new Graph(argv[1]);
    
    int num;
    
//    vector<int> a;
//    for(int i = 0;i<10;i++) a.push_back(i);
//
//    a.erase(std::remove(a.begin(),a.end(),3),a.end());
//
//    for(int i = 0;i<a.size();i++) cout<<a[i]<<" ";
    
    
    cin>>num;
    
    if(num == 1){
        startTime = clock();//计时开始
        graph->readGraph();
        endTime = clock();//计时结束
        cout << "The read graph time is: " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
        
        graph->creatIndex();
    }
    else if(num == 2){
        
        startTime = clock();//计时开始
        graph->readIndex();
        endTime = clock();//计时结束
        cout << "The readIndex time is: " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
        
        float a,b;
        int c;
        char d;
        
        cout<<"input?(y/n)";
        
        cin>>d;
        
        while(d == 'y' || d == 'Y'){
            cout<<"input:";
            
            cin>>a>>b>>c;
            
            startTime = clock();//计时开始
            graph->query(a,b,c);
            endTime = clock();//计时结束
            cout << "The query time is: " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
            
            cout<<"input?(y/n)";
            
            cin>>d;
        }
    }
    else{
        startTime = clock();//计时开始
        graph->readGraph_update();
        endTime = clock();//计时结束
        cout << "The readGraph_update time is: " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
        
        startTime = clock();//计时开始
        graph->readIndex_update();
        endTime = clock();//计时结束
        cout << "The readIndex_update time is: " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
        
        startTime = clock();//计时开始
        graph->insert(1,5);
        endTime = clock();//计时结束
        cout << "The insert time is: " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
        
        int a,b,c;
        char d;
        
        cout<<"input?(y/n)";
        
        cin>>d;
        
        while(d == 'y' || d == 'Y'){
            cout<<"input:";
            
            cin>>a>>b>>c;
            
            if(c){
                startTime = clock();//计时开始
                graph->insert(a, b);
                endTime = clock();//计时结束
                cout << "The insert time is: " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
            }else{
                startTime = clock();//计时开始
                graph->remove(a, b);
                endTime = clock();//计时结束
                cout << "The remove time is: " <<(double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
            }
            cout<<"input?(y/n)";
            
            cin>>d;
        }
        
        //graph->readIndex_update();
    }
    
    return 0;
}

