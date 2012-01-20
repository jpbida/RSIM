/*This file is part of RSIM.

    RSIM is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    RSIM is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with RSIM.  If not, see <http://www.gnu.org/licenses/>.

*/
// Montecarlo fractal kinetics simulator //
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "points.h"


void Points::translate(Points * p){
using namespace std;
//translation in direction of x2,y2,z2
double	x1 = p->x;
double	y1 = p->y;
double	z1 = p->z;
	x = x-x1;
	y = y-y1;
	z = z-z1;
}

double Points::dot(Points * p)
{
//X-convention rotations Z-rotation, X-rotation, Z-rotaiton 
using namespace std;
double x1=p->x;
double y1=p->y;
double z1=p->z;
double out=x1*x+y1*y+z1*z;
return out;
}

void Points::print(){
using namespace std;
cout << x << ":" << y << ":" << z << endl;
}

Points Points::transform(Matrix * m){
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

