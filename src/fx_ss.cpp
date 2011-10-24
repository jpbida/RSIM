// Montecarlo fractal kinetics simulator //
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "ss.h"
/*
SS::Base(double x=0, double y=0, double z = 0, double bx=0, double by=0, double bz = 0, double b2x=0, double b2y=0, double b2z = 0){
x=x;
y=y;
z=z;
bx=bx;
by=by;
bz=bz;
b2x=b2x;
b2y=b2y;
b2z=b2z;
} 
*/
void SS::translate(double x2, double y2, double z2){
//translation in direction of x2,y2,z2
	x = x-x2;
	y = y-y2;
	z = z-z2;
}

void SS::rotateX(double ang)
{
//X-convention rotations Z-rotation, X-rotation, Z-rotaiton 
using namespace std;
double yp=y*cos(ang)-z*sin(ang);
double zp=y*sin(ang)+z*cos(ang);
	y = yp;
	z = zp;
}

void SS::rotateY(double ang)
{
//X-convention rotations Z-rotation, X-rotation, Z-rotaiton 
	using namespace std;
double xp=x*cos(ang)+z*sin(ang);
double zp=z*cos(ang)-x*sin(ang);
	x = xp;
	z = zp;
}

void SS::rotateZ(double ang)
{
//X-convention rotations Z-rotation, X-rotation, Z-rotaiton 
using namespace std;
double xp=x*cos(ang)+y*sin(ang);
double yp= -1*x*sin(ang)+y*cos(ang);
x = xp;
y = yp;
}



void SS::rotate(double a,double b,double g,double x2,double y2,double z2)
{
	using namespace std;
	double xp = (cos(a)*cos(g)-sin(a)*cos(b)*sin(g))*(x-x2) + (sin(a)*cos(g)+cos(a)*cos(b)*sin(g))*(y-y2) + (sin(b)*sin(g))*(z-z2) + x2;
	double yp = (-cos(a)*sin(g)-sin(a)*cos(b)*cos(g))*(x-x2) + (-sin(a)*sin(g)+cos(a)*cos(b)*cos(g))*(y-y2) + (sin(b)*cos(g))*(z-z2) + y2;
	double zp = (sin(a)*sin(b))*(x-x2) - (cos(a)*sin(b))*(y-y2) + cos(b)*(z-z2) + z2;
	x = xp;
	y = yp;
	z = zp;
}
