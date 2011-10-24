#ifndef _PARSE
#define _PARSE
#include <istream>
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include <istream>
#include <sstream>
#include <vector>
void loadDoublets(std::vector<double> *,std::string fname);
double strtofloat(std::string what);
double strtodouble(std::string what);
int strtoint(std::string what);
char strtochar(std::string what);
double round_nplaces(double value, int to);
const double PI = 3.14159265358979323846;
#endif
