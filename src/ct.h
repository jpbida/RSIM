// molecule.h -- molecule class and functions
#include <vector>
#include "ss.h"
#ifndef CT_H_
#define CT_H_

class CT
{
private:
public:
std::vector<SS *> nucs;
int num_mols;
CT();
CT(std::string ss, std::string seq);
void layout();
void print();
int ct2rsim();
void translate(double,double,double);
void rotate(double a, double b, double g, double x2, double y2 , double z2);
void rotateX(double);
void rotateY(double);
void rotateZ(double);
//Returns center of mass for a given base
};




#endif



