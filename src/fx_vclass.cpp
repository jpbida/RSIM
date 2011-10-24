// Montecarlo fractal kinetics simulator //
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "vclass.h"

Vclass::Vclass(){
id=-1;
x=0;
y=0;
z=0;
vz=0;
vx=0;
vy=0;
rz=0;
rx=0;
ry=0;
group=-1;  // Cluster group used for indexing
}

void Vclass::translate(Vclass * p){
using namespace std;
//translation in direction of x2,y2,z2
double	x1 = p->x;
double	y1 = p->y;
double	z1 = p->z;
	x = x-x1;
	y = y-y1;
	z = z-z1;
}

double Vclass::dot(Vclass * p)
{
//X-convention rotations Z-rotation, X-rotation, Z-rotaiton 
using namespace std;
double x1=p->x;
double y1=p->y;
double z1=p->z;
double out=x1*x+y1*y+z1*z;
return out;
}

double Vclass::dotr(Vclass * p)
{
//X-convention rotations Z-rotation, X-rotation, Z-rotaiton 
using namespace std;
double x1=p->rx;
double y1=p->ry;
double z1=p->rz;
double out=x1*rx+y1*ry+z1*rz;
return out;
}

double Vclass::dotv(Vclass * p)
{
//X-convention rotations Z-rotation, X-rotation, Z-rotaiton 
using namespace std;
double x1=p->vx;
double y1=p->vy;
double z1=p->vz;
double out=x1*vx+y1*vy+z1*vz;
return out;
}

void  Vclass::cross(Vclass * p)
{
//X-convention rotations Z-rotation, X-rotation, Z-rotaiton 
using namespace std;
double b1=p->x;
double b2=p->y;
double b3=p->z;
double a1=x;
double a2=y;
double a3=z;
Vclass * out=new Vclass();
out->x=a2*b3-a3*b2;
out->y=a3*b1-a1*b3;
out->z=a1*b2-a2*b1;
}

void Vclass::print(){
using namespace std;
cout << id << "," << x << "," << y << "," << z << "," << vx << "," << vy << "," << vz << "," << rx << "," << ry << "," << rz << endl;
}

void Vclass::transform(Matrix * m){
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

x2=a11*rx + a12*ry + a13*rz;
y2=a21*rx + a22*ry + a23*rz;
z2=a31*rx + a32*ry + a33*rz;

rx = x2;
ry = y2;
rz = z2;

x2=a11*vx + a12*vy + a13*vz;
y2=a21*vx + a22*vy + a23*vz;
z2=a31*vx + a32*vy + a33*vz;

vx = x2;
vy = y2;
vz = z2;
}

