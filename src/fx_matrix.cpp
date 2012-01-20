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
#include "matrix.h"

Matrix Matrix::inverse(){
//Determines inverse of 3x3 matrix

double n11=a22*a33-a23*a32;
double n21=a23*a31-a21*a33;
double n31=a21*a32-a22*a31;

double n12=a13*a32-a12*a33;
double n22=a11*a33-a13*a31;
double n32=a12*a31-a11*a32;

double n13=a12*a23-a13*a22;
double n23=a13*a21-a11*a23;
double n33=a11*a22-a12*a21;

double det=a11*a22*a33 + a12*a23*a31 + a13*a21*a32 - a31*a22*a13 - a32*a23*a11 - a33*a21*a12;
Matrix m;
m.a11=n11/det;
m.a12=n12/det;
m.a13=n13/det;
m.a21=n21/det;
m.a22=n22/det;
m.a23=n23/det;
m.a31=n31/det;
m.a32=n32/det;
m.a33=n33/det;
return m;
}

Matrix Matrix::product(Matrix * Q){
//This X M
Matrix m;

m.a11=a11*Q->a11+a12*Q->a21+a13*Q->a31;
m.a21=a21*Q->a11+a22*Q->a21+a23*Q->a31;
m.a31=a31*Q->a11+a32*Q->a21+a33*Q->a31;

m.a12=a11*Q->a12+a12*Q->a22+a13*Q->a32;
m.a22=a21*Q->a12+a22*Q->a22+a23*Q->a32;
m.a32=a31*Q->a12+a32*Q->a22+a33*Q->a32;

m.a13=a11*Q->a13+a12*Q->a23+a13*Q->a33;
m.a23=a21*Q->a13+a22*Q->a23+a23*Q->a33;
m.a33=a31*Q->a13+a32*Q->a23+a33*Q->a33;

return m;
}
