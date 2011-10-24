// molecule.h -- molecule class and functions
#ifndef FABASE_H_
#define FABASE_H_
#include <iostream>
#include "points.h"
#include "matrix.h"
#include <cmath>
#include <vector>
#include <istream>
#include <sstream>

class FABase
{
private:
public:
FABase(int type);
int type;
Points N1;
Points N2;
Points N3;
Points N4;
Points N5;
Points N6;
Points N7;
Points N8;
Points N9;

Points C1;
Points C2;
Points C3;
Points C4;
Points C5;
Points C6;
Points C7;
Points C8;
Points C9;

Points O2;
Points O4;
Points O6;

//Hydrogens
Points H21;
Points H22;
Points H1;
Points H2;
Points H8;
Points H61;
Points H62;
Points H5;
Points H6;
Points H41;
Points H42;
Points H3;


//Lone Pairs
Points LP1;
Points LP3;
Points LP7;
Points LP61;
Points LP62;
Points LP21;
Points LP22;
Points LP41;
Points LP42;
//Centers of Mass 
Points CM1C;
Points CM1R;
Points CM2C;
Points CM2R;

//Backbone atoms
Points P;
Points O1P;
Points O2P;
Points O5p;
Points C5p;
Points C4p;
Points O4p;
Points C1p;
Points C2p;
Points O2p;
Points C3p;
Points O3p;

static const int sNbGauss=7;
static const float sProbH[7];
static const float sWeight [7];
static const float sMean [7][3];
static const float sCovarInv [7][3][3];
static const float sCovarDet [7];
//Uses the c1 carbon from the backbone and two atoms from base to create coordinate system 
//for the remaining atoms in the base

//Use ideal base geometry to determine base coordinates 
void alignBase(std::vector <Points> pt);
//void alignAll(std::vector <Points> pt,int frag_id,<Fragatoms *> frags_at);

void print();
void acceptors(std::vector<Points> * pt);
void donors(std::vector<Points> * pt);
void donorsFace(std::vector<int> * pt);
void acceptorsFace(std::vector<int> * pt);
float Eval (float dist, float angle_h, float angle_l, bool decision);
float overlap(FABase);
};

#endif



