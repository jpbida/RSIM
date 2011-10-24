// molecule.h -- molecule class and functions
#include <vector>
#include "polymere.h"
#ifndef DOUBLET_H_
#define DOUBLET_H_

class Doublet
{
private:
public:
//Contains a vector of points for the two bases
//P, C1′, C4′, C2, C4 and C6
float C1_2x;
float C2_2x;
float C4_2x;
float C6_2x;

float C1_2y;
float C2_2y;
float C4_2y;
float C6_2y;

float C1_2z;
float C2_2z;
float C4_2z;
float C6_2z;

float C1_1x;
float C2_1x;
float C4_1x;
float C6_1x;

float C1_1y;
float C2_1y;
float C4_1y;
float C6_1y;

float C1_1z;
float C2_1z;
float C4_1z;
float C6_1z;

//TYPE:  0 = R, 1 = Y
int type1;
int type2;
};




#endif



