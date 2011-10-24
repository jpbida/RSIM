#ifndef GRAPH_H_
#define GRAPH_H_
#include "polymere.h"
#include "element.h"
#include <vector>

using namespace std;
class Graph
{
private:
public:
vector<Element *> elements;
double gen_id(Polymere * pol);
double  search(double);
void sort();
void add(Element * ew);
void print();
};

#endif

