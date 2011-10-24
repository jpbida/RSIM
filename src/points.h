// molecule.h -- molecule class and functions

#ifndef POINTS_H_
#define POINTS_H_
#include "matrix.h"

class Points
{
private:
public:
double z;
double x;
double y;

void translate(Points * t);
void print();
double dot(Points * t);
Points cross(Points * t);
Points transform(Matrix * m);
};

#endif



