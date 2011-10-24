// molecule.h -- molecule class and functions

#ifndef MATRIX_H_
#define MATRIX_H_
class Matrix
{
private:
public:
double a11;
double a12;
double a13;
double a21;
double a22;
double a23;
double a31;
double a32;
double a33;

Matrix inverse();
Matrix product(Matrix *);
};




#endif



