// molecule.h -- molecule class and functions
#ifndef FRAGALL_H_
#define FRAGALL_H_
#include <string>
#include "points.h"
#include <vector>

using namespace std;
class Fragall
{
private:
public:
Fragall();
Fragall(string line);

//Data structures
vector<Points> atoms;
int id;

void print();
void loadLine(std::string);
int check();
};

#endif
