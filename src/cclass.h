// molecule.h -- molecule class and functions
#ifndef CCLASS_H_
#define CCLASS_H_
#include "points.h"
#include <vector>
#include "vclass.h"
class Cclass
{
private:
public:

double z;
double x;
double y;
double r;
std::vector<int> members;
//When creating the fragment classes use the group_id variable to push the pointers to
//the fragments into the members array 
int id;
//Vector of array indexes to each fragment
void translate(Vclass * t);
void print();
double dot(Cclass * t);
void cross(Cclass * t);
void transform(Matrix * m);
};

#endif



