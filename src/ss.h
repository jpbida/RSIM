// molecule.h -- molecule class and functions
#include <vector>
#ifndef SS_H_
#define SS_H_

class SS
{
private:
public:
double z;
double x;
double y;
int group;
int before;
int after;
char seq;
int pos;
int num_mols;
int bp;
void translate(double,double,double);
void rotate(double a, double b, double g, double x2, double y2 , double z2);
void rotateX(double);
void rotateY(double);
void rotateZ(double);
//Returns center of mass for a given base
};




#endif



