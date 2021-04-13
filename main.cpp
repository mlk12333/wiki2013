//
//  main.cpp
//  changeTxt
//
//  Created by 孟令凯 on 2021/3/10.
//

#include <iostream>
#include "string"
#include "Change.hpp"

using namespace std;
int main(int argc, const char * argv[]) {
    if(argc<3) return 1;
    // insert code here...
    std::cout << "Hello, World!\n";
    
//    string str = "/Users/milk/test_data/zata/p2p-Gnutella31.txt";
//    string output = "/Users/milk/test_data/zata/p2p-Gnutella31";
    
    Change change;
//    change.squenceTxt(output);
    change.resolveTxt(argv[1], argv[2]);
       
    return 0;
}
