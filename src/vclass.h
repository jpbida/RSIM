// molecule.h -- molecule class and functions
#ifndef VCLASS_H_
#define VCLASS_H_
#include "matrix.h"
class Vclass
{
private:
public:
Vclass();
//Fragments are identified by the position of the pos+3 C1' atom relative to the pos-1 C1' atom.
//x,y,z give this relative position.  
//vz is the unit vector from the pos+3 C1' to the pos+2 C1'
//vr is the unit vector from the pos+3 C1' to the pos+4 C1'
double z;
double x;
double y;

double vz;
double vx;
double vy;

double rz;
double rx;
double ry;
int id;
int group;  // Cluster group used for indexing



void translate(Vclass * t);
void print();
double dot(Vclass * t);
double dotv(Vclass * t);
double dotr(Vclass * t);
void cross(Vclass * t);
void transform(Matrix * m);
};

#endif



