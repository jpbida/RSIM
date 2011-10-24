//
//  example.cpp
//
#include <algorithm>
#include <functional>
#include <numeric>
#include "example.h"

using namespace std;

double average(vector<int> v) {  return (accumulate(v.begin(), v.end(), 0.0))/v.size();}

vector<double> half(const vector<double> & v) {
    
    	vector<double> w(v);
    
    	for( unsigned long i = 0; i < w.size(); i++)
        		w[i] /= 2.0;
    
    	return w; }

std::string TestMod(const std::string &s)   {return string("TestOutput: ") + s;} 

