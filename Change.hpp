//
//  Change.hpp
//  changeTxt
//
//  Created by 孟令凯 on 2021/3/13.
//

#ifndef Change_hpp
#define Change_hpp

#include <stdio.h>
#include <fstream>
#include "string"
#include "vector"
using namespace std;

class Change{
public:
    void squenceTxt(string str);
    void resolveTxt(const string str, const string output);
};

#endif /* Change_hpp */
