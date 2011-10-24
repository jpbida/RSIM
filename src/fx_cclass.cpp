// Montecarlo fractal kinetics simulator //
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "cclass.h"
#include "vclass.h"

void Cclass::translate(Vclass * p){
using namespace std;
//translation in direction of x2,y2,z2
double	x1 = p->x;
double	y1 = p->y;
double	z1 = p->z;
	x = x-x1;
	y = y-y1;
	z = z-z1;
}

double Cclass::dot(Cclass * p)
{
//X-convention rotations Z-rotation, X-rotation, Z-rotaiton 
using namespace std;
double x1=p->x;
double y1=p->y;
double z1=p->z;
double out=x1*x+y1*y+z1*z;
return out;
}

void  Cclass::cross(Cclass * p)
{
//X-convention rotations Z-rotation, X-rotation, Z-rotaiton 
using namespace std;
double b1=p->x;
double b2=p->y;
double b3=p->z;
double a1=x;
double a2=y;
double a3=z;
Cclass * out=new Cclass();
out->x=a2*b3-a3*b2;
out->y=a3*b1-a1*b3;
out->z=a1*b2-a2*b1;
}

void Cclass::print(){
using namespace std;
cout << id << "," << x << "," << y << "," << z << "," << r << endl;
cout << "Members:" << members.size() << endl;
for(int j=0; j<members.size(); j++){
cout << j << ":"  << members[j] << endl;
}

}

void Cclass::transform(Matrix * m){
//Matrix rows are r1,r2,r3

double a11 = m->a11;
double a12 = m->a12;
double a13 = m->a13;
double a21 = m->a21;
double a22 = m->a22;
double a23 = m->a23;
double a31 = m->a31;
double a32 = m->a32;
double a33 = m->a33;

double x2=a11*x + a12*y + a13*z;
double y2=a21*x + a22*y + a23*z;
double z2=a31*x + a32*y + a33*z;

x = x2;
y = y2;
z = z2;

}

