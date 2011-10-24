// molecule.h -- molecule class and functions
#include <vector>
#include "points.h"
#ifndef BASE_H_
#define BASE_H_

class Base
{
private:
public:
double z;
double x;
double y;
int group;
float prob;
char seq;
int id;
double bx;
double by;
double bz;
double b2x;
double b2y;
double b2z;

void translate(double,double,double);
void rotate(double a, double b, double g, double x2, double y2 , double z2);
void rotateX(double);
void rotateY(double);
void rotateZ(double);
//Returns center of mass for a given base
};




#endif



