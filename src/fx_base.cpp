// Montecarlo fractal kinetics simulator //
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "base.h"
/*
Base::Base(double x=0, double y=0, double z = 0, double bx=0, double by=0, double bz = 0, double b2x=0, double b2y=0, double b2z = 0){
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
void Base::translate(double x2, double y2, double z2){
//translation in direction of x2,y2,z2
	x = x-x2;
	y = y-y2;
	z = z-z2;
	bx = bx-x2;
	by = by-y2;
	bz = bz-z2;
	b2x = b2x-x2;
	b2y = b2y-y2;
	b2z = b2z-z2;
}

void Base::rotateX(double ang)
{
//X-convention rotations Z-rotation, X-rotation, Z-rotaiton 
using namespace std;
double yp=y*cos(ang)-z*sin(ang);
double zp=y*sin(ang)+z*cos(ang);

double byp=by*cos(ang)-bz*sin(ang);
double bzp=by*sin(ang)+bz*cos(ang);

double b2yp=b2y*cos(ang)-b2z*sin(ang);
double b2zp=b2y*sin(ang)+b2z*cos(ang);

	y = yp;
	z = zp;
	by = byp;
	bz = bzp;
	
	b2y = b2yp;
	b2z = b2zp;
}

void Base::rotateY(double ang)
{
//X-convention rotations Z-rotation, X-rotation, Z-rotaiton 
	using namespace std;
double xp=x*cos(ang)+z*sin(ang);
double zp=z*cos(ang)-x*sin(ang);
double bxp=bx*cos(ang)+bz*sin(ang);
double bzp=bz*cos(ang)-bx*sin(ang);
double b2xp=b2x*cos(ang)+b2z*sin(ang);
double b2zp=b2z*cos(ang)-b2x*sin(ang);
	x = xp;
	z = zp;
	bx = bxp;
	bz = bzp;
	b2x = b2xp;
	b2z = b2zp;
}

void Base::rotateZ(double ang)
{
//X-convention rotations Z-rotation, X-rotation, Z-rotaiton 
using namespace std;
double xp=x*cos(ang)+y*sin(ang);
double yp= -1*x*sin(ang)+y*cos(ang);
double bxp=bx*cos(ang)+by*sin(ang);
double byp= -1*bx*sin(ang)+by*cos(ang);
double b2xp=b2x*cos(ang)+b2y*sin(ang);
double b2yp= -1*b2x*sin(ang)+b2y*cos(ang);
x = xp;
y = yp;
bx = bxp;
by = byp;
b2x = b2xp;
b2y = b2yp;
}



void Base::rotate(double a,double b,double g,double x2,double y2,double z2)
{
	using namespace std;
	double xp = (cos(a)*cos(g)-sin(a)*cos(b)*sin(g))*(x-x2) + (sin(a)*cos(g)+cos(a)*cos(b)*sin(g))*(y-y2) + (sin(b)*sin(g))*(z-z2) + x2;
	double yp = (-cos(a)*sin(g)-sin(a)*cos(b)*cos(g))*(x-x2) + (-sin(a)*sin(g)+cos(a)*cos(b)*cos(g))*(y-y2) + (sin(b)*cos(g))*(z-z2) + y2;
	double zp = (sin(a)*sin(b))*(x-x2) - (cos(a)*sin(b))*(y-y2) + cos(b)*(z-z2) + z2;

	double bxp = (cos(a)*cos(g)-sin(a)*cos(b)*sin(g))*(bx-x2) + (sin(a)*cos(g)+cos(a)*cos(b)*sin(g))*(by-y2) + (sin(b)*sin(g))*(bz-z2) + x2;
	double byp = (-cos(a)*sin(g)-sin(a)*cos(b)*cos(g))*(bx-x2) + (-sin(a)*sin(g)+cos(a)*cos(b)*cos(g))*(by-y2) + (sin(b)*cos(g))*(bz-z2) + y2;
	double bzp = (sin(a)*sin(b))*(bx-x2) - (cos(a)*sin(b))*(by-y2) + cos(b)*(bz-z2) + z2;
	
	double b2xp = (cos(a)*cos(g)-sin(a)*cos(b)*sin(g))*(b2x-x2) + (sin(a)*cos(g)+cos(a)*cos(b)*sin(g))*(b2y-y2) + (sin(b)*sin(g))*(b2z-z2) + x2;
	double b2yp = (-cos(a)*sin(g)-sin(a)*cos(b)*cos(g))*(b2x-x2) + (-sin(a)*sin(g)+cos(a)*cos(b)*cos(g))*(b2y-y2) + (sin(b)*cos(g))*(b2z-z2) + y2;
	double b2zp = (sin(a)*sin(b))*(b2x-x2) - (cos(a)*sin(b))*(b2y-y2) + cos(b)*(b2z-z2) + z2;
	x = xp;
	y = yp;
	z = zp;
	bx = bxp;
	by = byp;
	bz = bzp;
	bx = b2xp;
	by = b2yp;
	bz = b2zp;
}
