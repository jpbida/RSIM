// molecule.h -- molecule class and functions
#ifndef FRAGATOMS_H_
#define FRAGATOMS_H_
#include <string>
#include "points.h"
#include <vector>

using namespace std;
class Fragatoms
{
private:
public:
Fragatoms();
Fragatoms(string line);

vector<Points> atoms;
int id;

void loadLine(string line);
void print();
};


#endif
