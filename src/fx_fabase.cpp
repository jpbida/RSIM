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
#include "math.h"
#include "base.h"
#include "points.h"
#include "fabase.h"

const float FABase::sProbH[7] = {
  1,
  0,
  0,
  0,
  0,
  0,
  0
};


const float FABase::sWeight[7] = {
  0.00752508,
  0.00989524,
  0.0259051,
  0.109535,
  0.120741,
  0.534611,
  0.191654
};

const float FABase::sMean[7][3] = {
  {0.100846, 2.45718, 2.25155},
  {8.78474, 1.13239, -0.0738701},
  {3.28653, 0.449332, 0.474348},
  {5.92268, -0.553811, 0.0362715},
  {5.06479, -0.443828, -0.424931},
  {6.52255, -0.165407, -0.0835928},
  {7.73618, -0.296942, -0.299877}
};
const float FABase::sCovarInv[7][3][3] = {
  {{2.80141, 1.04863, 0.88979},
   {1.04863, 2.37603, -0.59708},
   {0.88979, -0.59708, 2.58077}},
  {{0.173398, 0.292812, 0.0360642},
   {0.292812, 1.02071, 0.192637},
   {0.0360642, 0.192637, 1.75143}},
  {{8.89039, 4.47243, 4.42705},
   {4.47243, 3.16789, 2.61418},
   {4.42705, 2.61418, 3.14729}},
  {{3.1899, 0.842455, 0.862609},
   {0.842455, 0.753278, 0.317195},
   {0.862609, 0.317195, 0.839071}},
  {{11.7232, 13.8285, 11.7912},
   {13.8285, 20.5467, 11.29},
   {11.7912, 11.29, 18.2967}},
  {{0.907316, 0.522994, 0.614166},
   {0.522994, 3.27058, 0.547732},
   {0.614166, 0.547732, 3.37027}},
  {{2.18977, 0.416751, 0.437718},
   {0.416751, 1.10475, 0.0841468},
   {0.437718, 0.0841468, 1.06094}}
};

const float FABase::sCovarDet[7] = {
  0.096653,
  6.40514,
  0.157177,
  0.999763,
  0.00418066,
  0.126178,
  0.45749,
};

float FABase::Eval (float dist, float angle_h, float angle_l, bool decision)
{ 
  float x[3];
//    x[0] = log (pow (dist, 3));
//    x[1] = atanh (cos (angle_h));
//    x[2] = atanh (cos (angle_l));
  x[0] = dist;
  x[1] = angle_h;
  x[2] = angle_l;
  
  float p_x = 0;
  float p_h = 0; 
  for (int i = 0; i < sNbGauss; ++i) {
    float diff[3];
    diff[0] = x[0] - sMean[i][0];
    diff[1] = x[1] - sMean[i][1];
    diff[2] = x[2] - sMean[i][2];
    
    float tmp = exp ((diff[0] * (diff[0] * sCovarInv[i][0][0] +
                                 diff[1] * sCovarInv[i][1][0] +
                                 diff[2] * sCovarInv[i][2][0]) +
                      diff[1] * (diff[0] * sCovarInv[i][0][1] +
                                 diff[1] * sCovarInv[i][1][1] +
                                 diff[2] * sCovarInv[i][2][1]) +
                      diff[2] * (diff[0] * sCovarInv[i][0][2] +
                                 diff[1] * sCovarInv[i][1][2] +
                                 diff[2] * sCovarInv[i][2][2])) * -0.5);
   
    float prob;
    if (std::isnan (tmp) || fabs (sCovarDet[i]) < 0.0005) prob = 0;
    else prob = sWeight[i] * tmp / (pow (2 * M_PI, 1.5) * sqrt (sCovarDet[i]));

    p_x += prob;
    p_h += sProbH[i] * prob;
  }

  if (decision) return (p_h / p_x);
  else return (p_x);
}

void FABase::acceptors(std::vector<Points> * donors){
using namespace std;
//Returns vector of all donors for the base
if(type==1){
donors->push_back(N1);
donors->push_back(H1);
donors->push_back(C8);
donors->push_back(H8);
donors->push_back(N2);
donors->push_back(H21);
donors->push_back(N2);
donors->push_back(H22);
}

if(type==2){
donors->push_back(C2);
donors->push_back(H2);
donors->push_back(C8);
donors->push_back(H8);
donors->push_back(N6);
donors->push_back(H61);
donors->push_back(N6);
donors->push_back(H62);
}

if(type==4){
donors->push_back(C5);
donors->push_back(H5);
donors->push_back(C6);
donors->push_back(H6);
donors->push_back(N4);
donors->push_back(H41);
donors->push_back(N4);
donors->push_back(H42);
}
if(type==3){
donors->push_back(N3);
donors->push_back(H3);
donors->push_back(C5);
donors->push_back(H5);
donors->push_back(C6);
donors->push_back(H6);
}     
}

void FABase::acceptorsFace(std::vector<int> * face){
//Face Type W = 0, H = 1, S = 2, -1 other
if(type==1){
face->push_back(0);
face->push_back(0);
face->push_back(-1);
face->push_back(-1);
face->push_back(0);
face->push_back(0);
face->push_back(2);
face->push_back(2);
}

if(type==2){
face->push_back(-1);
face->push_back(-1);
face->push_back(-1);
face->push_back(-1);
face->push_back(0);
face->push_back(0);
face->push_back(1);
face->push_back(1);
}

if(type==4){
face->push_back(-1);
face->push_back(-1);
face->push_back(-1);
face->push_back(-1);
face->push_back(0);
face->push_back(0);
face->push_back(1);
face->push_back(1);
}
if(type==3){
face->push_back(0);
face->push_back(0);
face->push_back(-1);
face->push_back(-1);
face->push_back(-1);
face->push_back(-1);
}     
}

void FABase::donorsFace(std::vector<int> * face){
using namespace std;
//Returns vector of all donors for the base
if(type==1){
face->push_back(1);
face->push_back(1);
face->push_back(0);
face->push_back(0);
face->push_back(1);
face->push_back(1);
face->push_back(2);
face->push_back(2);
}

if(type==2){
face->push_back(0);
face->push_back(0);
face->push_back(2);
face->push_back(2);
face->push_back(1);
face->push_back(1);
}     
      
if(type==4){
face->push_back(0);
face->push_back(0);
face->push_back(2);
face->push_back(2);
face->push_back(0);
face->push_back(0);
}     
if(type==3){
face->push_back(2);
face->push_back(2);
face->push_back(0);
face->push_back(0);
face->push_back(0);
face->push_back(0);
face->push_back(1);
face->push_back(1);
}

}



void FABase::donors(std::vector<Points> * donors){
using namespace std;
//Returns vector of all donors for the base
if(type==1){
donors->push_back(N7);
donors->push_back(LP7);
donors->push_back(O6);
donors->push_back(LP61);
donors->push_back(O6);
donors->push_back(LP62);
donors->push_back(N3);
donors->push_back(LP3);
}

if(type==2){
donors->push_back(N1);
donors->push_back(LP1);
donors->push_back(N3);
donors->push_back(LP3);
donors->push_back(N7);
donors->push_back(LP7);
}     
      
if(type==4){
donors->push_back(N3);
donors->push_back(LP3);
donors->push_back(O2);
donors->push_back(LP21);
donors->push_back(O2);
donors->push_back(LP22);
}     
if(type==3){
donors->push_back(O2);
donors->push_back(LP21);
donors->push_back(O2);
donors->push_back(LP22);
donors->push_back(O4);
donors->push_back(LP41);
donors->push_back(O4);
donors->push_back(LP42);
}

}

//FAbase::FAbase(int type1,int frag_id,<Fragatoms *> frags)
//Have it file in backbone atoms from fragment database


FABase::FABase(int type1){
type=type1;
//Type 1 = G, 2 = A, 3 = U, 4 = C
//Use ideal values from nucleic acid database
//Guanine
if(type==1){
C1.x=0;C1.y=0;C1.z=0;
N9.x =   0.214 ;N9.y =  0.659 ;N9.z =  1.283;
C4.x =   0.254 ;C4.y =  2.014 ;C4.z =  1.509;
N3.x =   0.034 ;N3.y =  2.979 ;N3.z =  0.591;
C2.x =   0.142 ;C2.y =  4.190 ;C2.z =  1.110;
N2.x =  -0.047 ;N2.y =  5.269 ;N2.z =  0.336;
N1.x =   0.444 ;N1.y =  4.437 ;N1.z =  2.427;
C6.x =   0.676 ;C6.y =  3.459 ;C6.z =  3.389;
O6.x =   0.941 ;O6.y =  3.789 ;O6.z =  4.552;
C5.x =   0.562 ;C5.y =  2.154 ;C5.z =  2.846;
N7.x =   0.712 ;N7.y =  0.912 ;N7.z =  3.448;
C8.x =   0.498 ;C8.y =  0.057 ;C8.z =  2.485;
//Hydrogens
H1.x= 0.5014913;  H1.y=  5.3924780;  H1.z=  2.7164072 ;
H8.x= 0.5398955;  H8.y=  -1.013792;  H8.z=  2.6193450 ;
H21.x= 0.0346323 ; H21.y=6.1855960 ; H21.z=0.7273922  ;
H22.x=-0.2695422 ; H22.y=5.1568580 ; H22.z=-0.6324519 ;
LP3.x= -0.184759; LP3.y=  2.7990803; LP3.z=  -0.3680481; 
LP7.x= 0.9338183; LP7.y=  0.7165434; LP7.z=  4.4032975 ;
LP61.x=0.9886637 ;LP61.y=4.756954;LP61.z=  4.7985626;
LP62.x=1.107458  ;LP62.y=3.087688;LP62.z=  5.245147;
CM1C.x=0.3505; CM1C.y=3.2545; CM1C.z=1.984; 
CM1R.x=0.5515; CM1R.y=5.7315; CM1R.z=2.934; 
CM2C.x=0.48; CM2C.y=0.797; CM2C.z=2.3705; 
CM2R.x=0.9072; CM2R.y=1.0162; CM2R.z=4.1041; 

}

//Adenine
if(type==2){
C1.x=0;C1.y=0;C1.z=0;
N9.x = 0.213; N9.y =   0.660;N9.z =    1.287;
C4.x = 0.250; C4.y =   2.016;C4.z =    1.509;
N3.x = 0.016; N3.y =   2.995;N3.z =    0.619;
C2.x = 0.142; C2.y =   4.189;C2.z =    1.194;
N1.x = 0.451; N1.y =   4.493;N1.z =    2.459;
C6.x = 0.681; C6.y =   3.485;C6.z =    3.329;
N6.x = 0.990; N6.y =   3.787;N6.z =    4.592;
C5.x = 0.579; C5.y =   2.170;C5.z =    2.844;
N7.x = 0.747; N7.y =   0.934;N7.z =    3.454;
C8.x = 0.520; C8.y =   0.074;C8.z =    2.491;
//Hydrogens
 H2.x= -0.03030765; H2.y=  5.03467820;  H2.z= 0.54474007;
 H8.x=   0.5697118; H8.y=  -0.9957162;  H8.z= 2.6311283 ;

H61.x=1.052039    ;H61.y=4.820138   ; H61.z=  4.900509 ;
H62.x=1.177964    ;H62.y=2.998202   ; H62.z=  5.305350 ;
//Lone Pairs
LP1.x=  0.5106161;LP1.y=   5.4463413;LP1.z= 2.7549498;
LP3.x= -0.2199291;LP3.y=   2.8472923;LP3.z= -0.3414790;
LP7.x=  0.9856354;LP7.y=   0.7440680;LP7.z= 4.4063544;
CM1C.x=0.3505; CM1C.y=3.2545; CM1C.z=1.984; 
CM1R.x=0.5515; CM1R.y=5.7315; CM1R.z=2.934; 
CM2C.x=0.48; CM2C.y=0.797; CM2C.z=2.3705; 
CM2R.x=0.9072; CM2R.y=1.0162; CM2R.z=4.1041; 
}


//Uracil
if(type==3){ 
C1.x=0;C1.y=0;C1.z=0;
N1.x=0.212; N1.y=  0.676; N1.z=1.281;
C6.x=0.195; C6.y= -0.023; C6.z=2.466;
C2.x=0.370; C2.y=  2.048; C2.z=1.265;
O2.x=0.390; O2.y=  2.698; O2.z=0.235;
N3.x=0.505; N3.y=  2.629; N3.z=2.502;
C4.x=0.497; C4.y=  1.990; C4.z=3.725;
O4.x=0.629; O4.y=  2.653; O4.z=4.755;
C5.x=0.329; C5.y=  0.571; C5.z=3.657;
//Hydrogens
H3.x= 0.6216583;  H3.y=  3.6220401  ;  H3.z=   2.5181943  ;
H5.x= 0.3104567;  H5.y=  -0.01752211;  H5.z=    4.56237168;
H6.x= 0.0682911;  H6.y=  -1.09503985;  H6.z=    2.43320609;
LP21.x= 0.2890045 ;LP21.y=2.2392270   ;LP21.z= -0.6477951   ;
LP22.x=0.5074143  ;LP22.y=3.6903851   ;LP22.z= 0.2722252    ;
LP41.x= 0.7422079;LP41.y=  3.6450475  ;LP41.z=   4.6999932  ;
LP42.x= 0.6229325;LP42.y=  2.1990896  ;LP42.z=   5.6460267  ;
CM1C.x=0.48; CM1C.y=0.797; CM1C.z=2.3705; 
CM1R.x=0.9072; CM1R.y=1.0162; CM1R.z=4.1041; 
}
//Cystosine
if(type==4){
N1.x=0.212;N1.y=   0.668;N1.z=   1.294;
C6.x=0.193;C6.y=  -0.043;C6.z=   2.462;
C2.x=0.374;C2.y=   2.055;C2.z=   1.315;
O2.x=0.388;O2.y=   2.673;O2.z=   0.240;
N3.x=0.511;N3.y=   2.687;N3.z=   2.504;
C4.x=0.491;C4.y=   1.984;C4.z=   3.638;
N4.x=0.631;N4.y=   2.649;N4.z=   4.788;
C5.x=0.328;C5.y=   0.569;C5.z=   3.645;
H5.x= 0.31385315;  H5.y=   0.01099740 ;  H5.z=     4.56957178;
H6.x= 0.06710056;  H6.y=   -1.11514618;  H6.z=     2.42956710;
H41.x=0.7549068 ; H41.y=  3.7218371  ; H41.z=  4.7797747 ;
H42.x=0.6202849 ; H42.y=  2.1138237  ; H42.z=  5.7260147 ;
LP3.x= 0.627097  ; LP3.y=   3.679826   ; LP3.z=     2.532613  ;
LP21.x=  0.2843344     ;LP21.y=  2.1780018     ;LP21.z=  -0.6226878  ;
LP22.x= 0.5029554   ;LP22.y= 3.6663617   ;LP22.z=  0.2357933   ;
CM1C.x=0.48; CM1C.y=0.797; CM1C.z=2.3705; 
CM1R.x=0.9072; CM1R.y=1.0162; CM1R.z=4.1041; 
C1.x=0;C1.y=0;C1.z=0;
}                 
}



//Align to c1-n-c coordinate system
void FABase::alignBase(std::vector <Points> pts010){
//C1'
double ax=pts010[0].x;
double ay=pts010[0].y;
double az=pts010[0].z;
//N9/N1
double bx=pts010[1].x;
double by=pts010[1].y;
double bz=pts010[1].z;
//C4/C2
double cx=pts010[2].x;
double cy=pts010[2].y;
double cz=pts010[2].z;


//Create a forth point using the cross product of a-b X c-b + b
double A1=ax-bx;
double A2=ay-by;
double A3=az-bz;

double B1=cx-bx;
double B2=cy-by;
double B3=cz-bz;

//Scale B just in case threading is used
double Bd = sqrt(pow(B1,2)+pow(B2,2)+pow(B3,2));
if(type==1 || type==2 ){
B1=B1/Bd * 1.3745;
B2=B2/Bd * 1.3745;
B3=B3/Bd * 1.3745;
}else{
if(type==3 || type==4){
B1=B1/Bd * 1.3965;
B2=B2/Bd * 1.3965;
B3=B3/Bd * 1.3965;
}
}
double Z1=B2*A3-B3*A2;
double Z2=B3*A1-B1*A3;
double Z3=B1*A2-B2*A1;

//Forth point is cross product of A and B
double dx=Z1+bx;
double dy=Z2+by;
double dz=Z3+bz;

//Determine transformation Matrix 

Matrix * Q = new Matrix();
Q->a11=bx-ax;
Q->a21=cx-ax;
Q->a31=dx-ax;
Q->a12=by-ay;
Q->a22=cy-ay;
Q->a32=dy-ay;
Q->a13=bz-az;
Q->a23=cz-az;
Q->a33=dz-az;
/*
   (bx-ax1) (by1-ay1) (bz1-az1) 
   P=  (cx1-ax1) (cy1-ay1) (cz1-az1)
   (dx1-ax1) (dy1-ay1) (dz1-az1)
*/
double ax2=0;
double ay2=0;
double az2=0;
double bx2,by2,bz2,cx2,cy2,cz2;
if(type==1 || type==2){
bx2=N9.x;
by2=N9.y;
bz2=N9.z;
cx2=C4.x;
cy2=C4.y;
cz2=C4.z;
}else{
bx2=N1.x;
by2=N1.y;
bz2=N1.z;
cx2=C2.x;
cy2=C2.y;
cz2=C2.z;
}

//create the  point d
//Create a forth point using the cross product of a-b X c-b + b
 A1=ax2-bx2;
 A2=ay2-by2;
 A3=az2-bz2;

 B1=cx2-bx2;
 B2=cy2-by2;
 B3=cz2-bz2;

 Z1=B2*A3-B3*A2;
 Z2=B3*A1-B1*A3;
 Z3=B1*A2-B2*A1;

//Forth point is cross product of A and B
double dx2=Z1+bx2;
double dy2=Z2+by2;
double dz2=Z3+bz2;

 
Matrix * P = new Matrix();
P->a11=bx2-ax2;
P->a21=cx2-ax2;
P->a31=dx2-ax2;
P->a12=by2-ay2;
P->a22=cy2-ay2;
P->a32=dy2-ay2;
P->a13=bz2-az2;
P->a23=cz2-az2;
P->a33=dz2-az2;
/*

       (bx2-ax2) (by2-ay2) (bz2-az2) 
                                     
   Q = (cx2-ax2) (cy2-ay2) (cz2-az2) 
                                     
       (dx2-ax2) (dy2-ay2) (dz2-az2) 
*/
/*
cout << "Matrix P" << endl;
cout << P->a11 << ":" << P->a12<<":" << P->a13 << endl;
cout << P->a21 << ":" << P->a22<<":" << P->a23 << endl;
cout << P->a31 << ":" << P->a32<<":" << P->a33 << endl;

cout << "Inverse P" << endl;
*/
Matrix Ip = P->inverse();
/*
cout << Ip.a11 << ":" << Ip.a12<<":" << Ip.a13 << endl;
cout << Ip.a21 << ":" << Ip.a22<<":" << Ip.a23 << endl;
cout << Ip.a31 << ":" << Ip.a32<<":" << Ip.a33 << endl;

cout << "Matrix Q" << endl;
cout << Q->a11 << ":" << Q->a12<<":" << Q->a13 << endl;
cout << Q->a21 << ":" << Q->a22<<":" << Q->a23 << endl;
cout << Q->a31 << ":" << Q->a32<<":" << Q->a33 << endl;
*/
Matrix * T = new Matrix();

T->a11=Ip.a11*Q->a11+Ip.a12*Q->a21+Ip.a13*Q->a31;
T->a21=Ip.a21*Q->a11+Ip.a22*Q->a21+Ip.a23*Q->a31;
T->a31=Ip.a31*Q->a11+Ip.a32*Q->a21+Ip.a33*Q->a31;

T->a12=Ip.a11*Q->a12+Ip.a12*Q->a22+Ip.a13*Q->a32;
T->a22=Ip.a21*Q->a12+Ip.a22*Q->a22+Ip.a23*Q->a32;
T->a32=Ip.a31*Q->a12+Ip.a32*Q->a22+Ip.a33*Q->a32;

T->a13=Ip.a11*Q->a13+Ip.a12*Q->a23+Ip.a13*Q->a33;
T->a23=Ip.a21*Q->a13+Ip.a22*Q->a23+Ip.a23*Q->a33;
T->a33=Ip.a31*Q->a13+Ip.a32*Q->a23+Ip.a33*Q->a33;
/*
cout << "Matrix T" << endl;
cout << T->a11 << ":" << T->a12<<":" << T->a13 << endl;
cout << T->a21 << ":" << T->a22<<":" << T->a23 << endl;
cout << T->a31 << ":" << T->a32<<":" << T->a33 << endl;
*/
//Determine the translation vector
double ax3= ax-(ax2*T->a11+ay2*T->a21+az2*T->a31);
double ay3= ay-(ax2*T->a12+ay2*T->a22+az2*T->a32);
double az3= az-(ax2*T->a13+ay2*T->a23+az2*T->a33);

/*
The geometric condition on A,B,C,D not lying in a common plane amounts
to the matrix algebra condition that both P and Q are invertible, so:
         M  = P   * Q 
*/
//Convert each point to new coordinate system
//pts[1-3]=NT1-c1'
//pts[4-6]=NT2-c1'
//pts[7-9]=NT3-c1'
//pts[10-12]=NT1b-c
//pts[13-15]=NT2b-c
//pts[16-18]=NT3b-c
//pts[19-21]=NT1b-n
//pts[22-24]=NT2b-n
//pts[25-27]=NT3b-n
//Update points in base
//Update all points in the base
double x,y,z,w1,w2,w3;
x=C1.x;y=C1.y;z=C1.z;
w1=x*T->a11+y*T->a21+z*T->a31+ax3;
w2=x*T->a12+y*T->a22+z*T->a32+ay3;
w3=x*T->a13+y*T->a23+z*T->a33+az3;
C1.x=w1;C1.y=w2;C1.z=w3;

x=C2.x;y=C2.y;z=C2.z;
w1=x*T->a11+y*T->a21+z*T->a31+ax3;
w2=x*T->a12+y*T->a22+z*T->a32+ay3;
w3=x*T->a13+y*T->a23+z*T->a33+az3;
C2.x=w1;C2.y=w2;C2.z=w3;

x=C3.x;y=C3.y;z=C3.z;
w1=x*T->a11+y*T->a21+z*T->a31+ax3;
w2=x*T->a12+y*T->a22+z*T->a32+ay3;
w3=x*T->a13+y*T->a23+z*T->a33+az3;
C3.x=w1;C3.y=w2;C3.z=w3;

x=C4.x;y=C4.y;z=C4.z;
w1=x*T->a11+y*T->a21+z*T->a31+ax3;
w2=x*T->a12+y*T->a22+z*T->a32+ay3;
w3=x*T->a13+y*T->a23+z*T->a33+az3;
C4.x=w1;C4.y=w2;C4.z=w3;

x=C5.x;y=C5.y;z=C5.z;
w1=x*T->a11+y*T->a21+z*T->a31+ax3;
w2=x*T->a12+y*T->a22+z*T->a32+ay3;
w3=x*T->a13+y*T->a23+z*T->a33+az3;
C5.x=w1;C5.y=w2;C5.z=w3;

x=C6.x;y=C6.y;z=C6.z;
w1=x*T->a11+y*T->a21+z*T->a31+ax3;
w2=x*T->a12+y*T->a22+z*T->a32+ay3;
w3=x*T->a13+y*T->a23+z*T->a33+az3;
C6.x=w1;C6.y=w2;C6.z=w3;

x=C7.x;y=C7.y;z=C7.z;
w1=x*T->a11+y*T->a21+z*T->a31+ax3;
w2=x*T->a12+y*T->a22+z*T->a32+ay3;
w3=x*T->a13+y*T->a23+z*T->a33+az3;
C7.x=w1;C7.y=w2;C7.z=w3;

x=C8.x;y=C8.y;z=C8.z;
w1=x*T->a11+y*T->a21+z*T->a31+ax3;
w2=x*T->a12+y*T->a22+z*T->a32+ay3;
w3=x*T->a13+y*T->a23+z*T->a33+az3;
C8.x=w1;C8.y=w2;C8.z=w3;

x=C9.x;y=C9.y;z=C9.z;
w1=x*T->a11+y*T->a21+z*T->a31+ax3;
w2=x*T->a12+y*T->a22+z*T->a32+ay3;
w3=x*T->a13+y*T->a23+z*T->a33+az3;
C9.x=w1;C9.y=w2;C9.z=w3;

x=N1.x;y=N1.y;z=N1.z;
w1=x*T->a11+y*T->a21+z*T->a31+ax3;
w2=x*T->a12+y*T->a22+z*T->a32+ay3;
w3=x*T->a13+y*T->a23+z*T->a33+az3;
N1.x=w1;N1.y=w2;N1.z=w3;

x=N2.x;y=N2.y;z=N2.z;
w1=x*T->a11+y*T->a21+z*T->a31+ax3;
w2=x*T->a12+y*T->a22+z*T->a32+ay3;
w3=x*T->a13+y*T->a23+z*T->a33+az3;
N2.x=w1;N2.y=w2;N2.z=w3;

x=N3.x;y=N3.y;z=N3.z;
w1=x*T->a11+y*T->a21+z*T->a31+ax3;
w2=x*T->a12+y*T->a22+z*T->a32+ay3;
w3=x*T->a13+y*T->a23+z*T->a33+az3;
N3.x=w1;N3.y=w2;N3.z=w3;

x=N4.x;y=N4.y;z=N4.z;
w1=x*T->a11+y*T->a21+z*T->a31+ax3;
w2=x*T->a12+y*T->a22+z*T->a32+ay3;
w3=x*T->a13+y*T->a23+z*T->a33+az3;
N4.x=w1;N4.y=w2;N4.z=w3;

x=N5.x;y=N5.y;z=N5.z;
w1=x*T->a11+y*T->a21+z*T->a31+ax3;
w2=x*T->a12+y*T->a22+z*T->a32+ay3;
w3=x*T->a13+y*T->a23+z*T->a33+az3;
N5.x=w1;N5.y=w2;N5.z=w3;

x=N6.x;y=N6.y;z=N6.z;
w1=x*T->a11+y*T->a21+z*T->a31+ax3;
w2=x*T->a12+y*T->a22+z*T->a32+ay3;
w3=x*T->a13+y*T->a23+z*T->a33+az3;
N6.x=w1;N6.y=w2;N6.z=w3;

x=N7.x;y=N7.y;z=N7.z;
w1=x*T->a11+y*T->a21+z*T->a31+ax3;
w2=x*T->a12+y*T->a22+z*T->a32+ay3;
w3=x*T->a13+y*T->a23+z*T->a33+az3;
N7.x=w1;N7.y=w2;N7.z=w3;

//I don't think this is needed 
//x=N8.x;y=N8.y;z=N8.z;
//w1=x*T->a11+y*T->a21+z*T->a31+ax3;
//w2=x*T->a12+y*T->a22+z*T->a32+ay3;
//w3=x*T->a13+y*T->a23+z*T->a33+az3;
//N8.x=w1;N8.y=w2;N8.z=w3;

x=N9.x;y=N9.y;z=N9.z;
w1=x*T->a11+y*T->a21+z*T->a31+ax3;
w2=x*T->a12+y*T->a22+z*T->a32+ay3;
w3=x*T->a13+y*T->a23+z*T->a33+az3;
N9.x=w1;N9.y=w2;N9.z=w3;


x=O2.x;y=O2.y;z=O2.z;
w1=x*T->a11+y*T->a21+z*T->a31+ax3;
w2=x*T->a12+y*T->a22+z*T->a32+ay3;
w3=x*T->a13+y*T->a23+z*T->a33+az3;
O2.x=w1;O2.y=w2;O2.z=w3;

x=O6.x;y=O6.y;z=O6.z;
w1=x*T->a11+y*T->a21+z*T->a31+ax3;
w2=x*T->a12+y*T->a22+z*T->a32+ay3;
w3=x*T->a13+y*T->a23+z*T->a33+az3;
O6.x=w1;O6.y=w2;O6.z=w3;

x=O4.x;y=O4.y;z=O4.z;
w1=x*T->a11+y*T->a21+z*T->a31+ax3;
w2=x*T->a12+y*T->a22+z*T->a32+ay3;
w3=x*T->a13+y*T->a23+z*T->a33+az3;
O4.x=w1;O4.y=w2;O4.z=w3;


x=H21.x;y=H21.y;z=H21.z;
w1=x*T->a11+y*T->a21+z*T->a31+ax3;
w2=x*T->a12+y*T->a22+z*T->a32+ay3;
w3=x*T->a13+y*T->a23+z*T->a33+az3;
H21.x=w1;H21.y=w2;H21.z=w3;

x=LP42.x;y=LP42.y;z=LP42.z;
w1=x*T->a11+y*T->a21+z*T->a31+ax3;
w2=x*T->a12+y*T->a22+z*T->a32+ay3;
w3=x*T->a13+y*T->a23+z*T->a33+az3;
LP42.x=w1;LP42.y=w2;LP42.z=w3;

x=LP41.x;y=LP41.y;z=LP41.z;
w1=x*T->a11+y*T->a21+z*T->a31+ax3;
w2=x*T->a12+y*T->a22+z*T->a32+ay3;
w3=x*T->a13+y*T->a23+z*T->a33+az3;
LP41.x=w1;LP41.y=w2;LP41.z=w3;

x=LP22.x;y=LP22.y;z=LP22.z;
w1=x*T->a11+y*T->a21+z*T->a31+ax3;
w2=x*T->a12+y*T->a22+z*T->a32+ay3;
w3=x*T->a13+y*T->a23+z*T->a33+az3;
LP22.x=w1;LP22.y=w2;LP22.z=w3;

x=LP21.x;y=LP21.y;z=LP21.z;
w1=x*T->a11+y*T->a21+z*T->a31+ax3;
w2=x*T->a12+y*T->a22+z*T->a32+ay3;
w3=x*T->a13+y*T->a23+z*T->a33+az3;
LP21.x=w1;LP21.y=w2;LP21.z=w3;

x=LP61.x;y=LP61.y;z=LP61.z;
w1=x*T->a11+y*T->a21+z*T->a31+ax3;
w2=x*T->a12+y*T->a22+z*T->a32+ay3;
w3=x*T->a13+y*T->a23+z*T->a33+az3;
LP61.x=w1;LP61.y=w2;LP61.z=w3;

x=LP62.x;y=LP62.y;z=LP62.z;
w1=x*T->a11+y*T->a21+z*T->a31+ax3;
w2=x*T->a12+y*T->a22+z*T->a32+ay3;
w3=x*T->a13+y*T->a23+z*T->a33+az3;
LP62.x=w1;LP62.y=w2;LP62.z=w3;

x=LP7.x;y=LP7.y;z=LP7.z;
w1=x*T->a11+y*T->a21+z*T->a31+ax3;
w2=x*T->a12+y*T->a22+z*T->a32+ay3;
w3=x*T->a13+y*T->a23+z*T->a33+az3;
LP7.x=w1;LP7.y=w2;LP7.z=w3;

x=LP3.x;y=LP3.y;z=LP3.z;
w1=x*T->a11+y*T->a21+z*T->a31+ax3;
w2=x*T->a12+y*T->a22+z*T->a32+ay3;
w3=x*T->a13+y*T->a23+z*T->a33+az3;
LP3.x=w1;LP3.y=w2;LP3.z=w3;

x=LP1.x;y=LP1.y;z=LP1.z;
w1=x*T->a11+y*T->a21+z*T->a31+ax3;
w2=x*T->a12+y*T->a22+z*T->a32+ay3;
w3=x*T->a13+y*T->a23+z*T->a33+az3;
LP1.x=w1;LP1.y=w2;LP1.z=w3;

x=H3.x;y=H3.y;z=H3.z;
w1=x*T->a11+y*T->a21+z*T->a31+ax3;
w2=x*T->a12+y*T->a22+z*T->a32+ay3;
w3=x*T->a13+y*T->a23+z*T->a33+az3;
H3.x=w1;H3.y=w2;H3.z=w3;

x=H42.x;y=H42.y;z=H42.z;
w1=x*T->a11+y*T->a21+z*T->a31+ax3;
w2=x*T->a12+y*T->a22+z*T->a32+ay3;
w3=x*T->a13+y*T->a23+z*T->a33+az3;
H42.x=w1;H42.y=w2;H42.z=w3;

x=H41.x;y=H41.y;z=H41.z;
w1=x*T->a11+y*T->a21+z*T->a31+ax3;
w2=x*T->a12+y*T->a22+z*T->a32+ay3;
w3=x*T->a13+y*T->a23+z*T->a33+az3;
H41.x=w1;H41.y=w2;H41.z=w3;


x=H6.x;y=H6.y;z=H6.z;
w1=x*T->a11+y*T->a21+z*T->a31+ax3;
w2=x*T->a12+y*T->a22+z*T->a32+ay3;
w3=x*T->a13+y*T->a23+z*T->a33+az3;
H6.x=w1;H6.y=w2;H6.z=w3;

x=H5.x;y=H5.y;z=H5.z;
w1=x*T->a11+y*T->a21+z*T->a31+ax3;
w2=x*T->a12+y*T->a22+z*T->a32+ay3;
w3=x*T->a13+y*T->a23+z*T->a33+az3;
H5.x=w1;H5.y=w2;H5.z=w3;

x=H62.x;y=H62.y;z=H62.z;
w1=x*T->a11+y*T->a21+z*T->a31+ax3;
w2=x*T->a12+y*T->a22+z*T->a32+ay3;
w3=x*T->a13+y*T->a23+z*T->a33+az3;
H62.x=w1;H62.y=w2;H62.z=w3;

x=H61.x;y=H61.y;z=H61.z;
w1=x*T->a11+y*T->a21+z*T->a31+ax3;
w2=x*T->a12+y*T->a22+z*T->a32+ay3;
w3=x*T->a13+y*T->a23+z*T->a33+az3;
H61.x=w1;H61.y=w2;H61.z=w3;

x=H8.x;y=H8.y;z=H8.z;
w1=x*T->a11+y*T->a21+z*T->a31+ax3;
w2=x*T->a12+y*T->a22+z*T->a32+ay3;
w3=x*T->a13+y*T->a23+z*T->a33+az3;
H8.x=w1;H8.y=w2;H8.z=w3;

x=H2.x;y=H2.y;z=H2.z;
w1=x*T->a11+y*T->a21+z*T->a31+ax3;
w2=x*T->a12+y*T->a22+z*T->a32+ay3;
w3=x*T->a13+y*T->a23+z*T->a33+az3;
H2.x=w1;H2.y=w2;H2.z=w3;

x=H1.x;y=H1.y;z=H1.z;
w1=x*T->a11+y*T->a21+z*T->a31+ax3;
w2=x*T->a12+y*T->a22+z*T->a32+ay3;
w3=x*T->a13+y*T->a23+z*T->a33+az3;
H1.x=w1;H1.y=w2;H1.z=w3;

x=H22.x;y=H22.y;z=H22.z;
w1=x*T->a11+y*T->a21+z*T->a31+ax3;
w2=x*T->a12+y*T->a22+z*T->a32+ay3;
w3=x*T->a13+y*T->a23+z*T->a33+az3;
H22.x=w1;H22.y=w2;H22.z=w3;

//Transform all coordinates
delete Q;
delete P;
delete T;
} 


float FABase::overlap(FABase b2){
using namespace std;
//Projects b2 into the base's plane and calculates the area of overlap 
double mx,my,mz,mbx,mby,mbz,mb2x,mb2y,mb2z;
mx=C1.x;
my=C1.y;
mz=C1.z;
if(type < 3){
mbx=N9.x;
mby=N9.y;
mbz=N9.z;
mb2x=C4.x;
mb2y=C4.y;
mb2z=C4.z;
}
if(type > 2){
mbx=N1.x;
mby=N1.y;
mbz=N1.z;
mb2x=C2.x;
mb2y=C2.y;
mb2z=C2.z;
}
double Ax1=mx - mbx;
double Ay1=my - mby;
double Az1=mz - mbz;

double Bx1=mb2x - mbx;
double By1=mb2y - mby;
double Bz1=mb2z - mbz;
mx=b2.C1.x;
my=b2.C1.y;
mz=b2.C1.z;
if(b2.type < 3){
mbx=b2.N9.x;
mby=b2.N9.y;
mbz=b2.N9.z;
mb2x=b2.C4.x;
mb2y=b2.C4.y;
mb2z=b2.C4.z;
}
if(b2.type > 2){
mbx=b2.N1.x;
mby=b2.N1.y;
mbz=b2.N1.z;
mb2x=b2.C2.x;
mb2y=b2.C2.y;
mb2z=b2.C2.z;
}

double Ax2=mx - mbx;
double Ay2=my - mby;
double Az2=mz - mbz;

double Bx2=mb2x - mbx;
double By2=mb2y - mby;
double Bz2=mb2z - mbz;

double Cx1=Ay1*Bz1 - Az1*By1;
double Cy1=Az1*Bx1 - Ax1*Bz1;
double Cz1=Ax1*By1 - Ay1*Bx1;

//Making orthonormal set
double Ex1=Cy1*Az1 - Cz1*Ay1;
double Ey1=Cz1*Ax1 - Cx1*Az1;
double Ez1=Cx1*Ay1 - Cy1*Ax1;

double Cx2 =Ay2*Bz2 - Az2*By2;
double Cy2 =Az2*Bx2 - Ax2*Bz2;
double Cz2 =Ax2*By2 - Ay2*Bx2;

double Ex2=Cy2*Az2 - Cz2*Ay2;
double Ey2=Cz2*Ax2 - Cx2*Az2;
double Ez2=Cx2*Ay2 - Cy2*Ax2;

//A,E,C is a coordinate frame
//Determine theta1,phi1,rho,theta2,phi2 
//plug into energy function
//Caclulate angle between normal vectors
double CC1=sqrt(pow(Cx1,2)+pow(Cy1,2)+pow(Cz1,2));
double CC2=sqrt(pow(Cx2,2)+pow(Cy2,2)+pow(Cz2,2));
double AA1=sqrt(pow(Ax1,2)+pow(Ay1,2)+pow(Az1,2));
double AA2=sqrt(pow(Ax2,2)+pow(Ay2,2)+pow(Az2,2));
double EE1=sqrt(pow(Ex1,2)+pow(Ey1,2)+pow(Ez1,2));
double EE2=sqrt(pow(Ex2,2)+pow(Ey2,2)+pow(Ez2,2));
//Create orthonormal set
Cx1=Cx1/CC1;
Cy1=Cy1/CC1;
Cz1=Cz1/CC1;

Ax1=Ax1/AA1;
Ay1=Ay1/AA1;
Az1=Az1/AA1;

Ex1=Ex1/EE1;
Ey1=Ey1/EE1;
Ez1=Ez1/EE1;

Cx2=Cx2/CC2;
Cy2=Cy2/CC2;
Cz2=Cz2/CC2;

Ax2=Ax2/AA2;
Ay2=Ay2/AA2;
Az2=Az2/AA2;
       
Ex2=Ex2/EE2;
Ey2=Ey2/EE2;
Ez2=Ez2/EE2;
//Transform CM into base coordinates
vector <Points> p1;
vector <Points> p2;
Points c1;
c1.x = b2.CM1C.x;
c1.y = b2.CM1C.y;
c1.z = b2.CM1C.z;

Points c2;
c2.x = b2.CM1R.x;
c2.y = b2.CM1R.y;
c2.z = b2.CM1R.z;

Points c3;
c3.x = b2.CM2C.x;
c3.y = b2.CM2C.y;
c3.z = b2.CM2C.z;

Points c4;
c4.x = b2.CM2R.x;
c4.y = b2.CM2R.y;
c4.z = b2.CM2R.z;

if(b2.type < 3){
double dx_1=(b2.CM1C.x*Ax1+b2.CM1C.y*Ay1+b2.CM1C.z*Az1);
double dy_1=(b2.CM1C.x*Ex1+b2.CM1C.y*Ey1+b2.CM1C.z*Ez1);
double dz_1=(b2.CM1C.x*Cx1+b2.CM1C.y*Cy1+b2.CM1C.z*Cz1);

double dx_2=(b2.CM1R.x*Ax1+b2.CM1R.y*Ay1+b2.CM1R.z*Az1);
double dy_2=(b2.CM1R.x*Ex1+b2.CM1R.y*Ey1+b2.CM1R.z*Ez1);
double dz_2=(b2.CM1R.x*Cx1+b2.CM1R.y*Cy1+b2.CM1R.z*Cz1);

double dx_4=(b2.CM2C.x*Ax1+b2.CM2C.y*Ay1+b2.CM2C.z*Az1);
double dy_4=(b2.CM2C.x*Ex1+b2.CM2C.y*Ey1+b2.CM2C.z*Ez1);
double dz_4=(b2.CM2C.x*Cx1+b2.CM2C.y*Cy1+b2.CM2C.z*Cz1);

double dx_3=(b2.CM2R.x*Ax1+b2.CM2R.y*Ay1+b2.CM2R.z*Az1);
double dy_3=(b2.CM2R.x*Ex1+b2.CM2R.y*Ey1+b2.CM2R.z*Ez1);
double dz_3=(b2.CM2R.x*Cx1+b2.CM2R.y*Cy1+b2.CM2R.z*Cz1);
cout << "R," 
<< dx_1 << ","  << dy_1 << ","  << dz_1 << "," 
<< dx_2 << ","  << dy_2 << ","  << dz_2 << "," 
<< dx_3 << ","  << dy_3 << ","  << dz_3 << "," 
<< dx_4 << ","  << dy_4 << ","  << dz_4 << "," 
<< endl;
}else{
double dx_1=(b2.CM1C.x*Ax1+b2.CM1C.y*Ay1+b2.CM1C.z*Az1);
double dy_1=(b2.CM1C.x*Ex1+b2.CM1C.y*Ey1+b2.CM1C.z*Ez1);
double dz_1=(b2.CM1C.x*Cx1+b2.CM1C.y*Cy1+b2.CM1C.z*Cz1);

double dx_2=(b2.CM1R.x*Ax1+b2.CM1R.y*Ay1+b2.CM1R.z*Az1);
double dy_2=(b2.CM1R.x*Ex1+b2.CM1R.y*Ey1+b2.CM1R.z*Ez1);
double dz_2=(b2.CM1R.x*Cx1+b2.CM1R.y*Cy1+b2.CM1R.z*Cz1);

cout << "Y," 
<< dx_1 << ","  << dy_1 << ","  << dz_1 << "," 
<< dx_2 << ","  << dy_2 << ","  << dz_2 << "," 
<< endl;
}
return 1;
}


void FABase::print(){
using namespace std;
//Type 1 = G, 2 = A, 3 = U, 4 = C
//Use ideal values from nucleic acid database
//Guanine
if(type==1){
cout<<N9.x  << "," <<N9.y   << ","<<N9.z    <<endl;
cout<<C4.x  << "," <<C4.y   << ","<<C4.z    <<endl;
cout<<N3.x  << "," <<N3.y   << ","<<N3.z    <<endl;
cout<<C2.x  << "," <<C2.y   << ","<<C2.z    <<endl;
cout<<N2.x  << "," <<N2.y   << ","<<N2.z    <<endl;
cout<<N1.x  << "," <<N1.y   << ","<<N1.z    <<endl;
cout<<C6.x  << "," <<C6.y   << ","<<C6.z    <<endl;
cout<<O6.x  << "," <<O6.y   << ","<<O6.z    <<endl;
cout<<C5.x  << "," <<C5.y   << ","<<C5.z    <<endl;
cout<<N7.x  << "," <<N7.y   << ","<<N7.z    <<endl;
cout<<C8.x  << "," <<C8.y   << ","<<C8.z    <<endl;
cout<<H1.x  << "," <<H1.y   << ","<<H1.z    <<endl;
cout<<H8.x  << "," <<H8.y   << ","<<H8.z    <<endl;
cout<<H21.x << "," <<H21.y  << ","<<H21.z   <<endl;
cout<<H22.x << "," <<H22.y  << ","<<H22.z   <<endl;
cout<<LP3.x << "," <<LP3.y  << ","<<LP3.z   <<endl;
cout<<LP7.x << "," <<LP7.y  << ","<<LP7.z   <<endl;
cout<<LP61.x<< "," <<LP61.y << ","<<LP61.z  <<endl;
cout<<LP62.x<< "," <<LP62.y << ","<<LP62.z  <<endl;
}

//Adenine
if(type==2){
cout<<N9.x << "," << N9.y << "," << N9.z  <<endl;
cout<<C4.x << "," << C4.y << "," << C4.z  <<endl;
cout<<N3.x << "," << N3.y << "," << N3.z  <<endl;
cout<<C2.x << "," << C2.y << "," << C2.z  <<endl;
cout<<N1.x << "," << N1.y << "," << N1.z  <<endl;
cout<<C6.x << "," << C6.y << "," << C6.z  <<endl;
cout<<N6.x << "," << N6.y << "," << N6.z  <<endl;
cout<<C5.x << "," << C5.y << "," << C5.z  <<endl;
cout<<N7.x << "," << N7.y << "," << N7.z  <<endl;
cout<<C8.x << "," << C8.y << "," << C8.z  <<endl;
cout<<H2.x << "," << H2.y << "," << H2.z  <<endl;
cout<<H8.x << "," << H8.y << "," << H8.z  <<endl;
cout<<H61.x<< "," << H61.y<< "," << H61.z <<endl;
cout<<H62.x<< "," << H62.y<< "," << H62.z <<endl;
cout<<LP1.x<< "," << LP1.y<< "," << LP1.z <<endl;
cout<<LP3.x<< "," << LP3.y<< "," << LP3.z <<endl;
cout<<LP7.x<< "," << LP7.y<< "," << LP7.z <<endl;
}

//Uracil
if(type==3){ 
cout<<N1.x  << "," << N1.y   << "," << N1.z   <<endl;
cout<<C6.x  << "," << C6.y   << "," << C6.z   <<endl;
cout<<C2.x  << "," << C2.y   << "," << C2.z   <<endl;
cout<<O2.x  << "," << O2.y   << "," << O2.z   <<endl;
cout<<N3.x  << "," << N3.y   << "," << N3.z   <<endl;
cout<<C4.x  << "," << C4.y   << "," << C4.z   <<endl;
cout<<O4.x  << "," << O4.y   << "," << O4.z   <<endl;
cout<<C5.x  << "," << C5.y   << "," << C5.z   <<endl;
cout<<H3.x  << "," << H3.y   << "," << H3.z   <<endl;
cout<<H5.x  << "," << H5.y   << "," << H5.z   <<endl;
cout<<H6.x  << "," << H6.y   << "," << H6.z   <<endl;
cout<<LP21.x<< "," << LP21.y << "," << LP21.z <<endl;
cout<<LP22.x<< "," << LP22.y << "," << LP22.z <<endl;
cout<<LP41.x<< "," << LP41.y << "," << LP41.z <<endl;
cout<<LP42.x<< "," << LP42.y << "," << LP42.z <<endl;
}             
             
//Cystosine
if(type==4){
cout<<N1.x  << "," << N1.y   << "," <<N1.z   <<endl;
cout<<C6.x  << "," << C6.y   << "," <<C6.z   <<endl;
cout<<C2.x  << "," << C2.y   << "," <<C2.z   <<endl;
cout<<O2.x  << "," << O2.y   << "," <<O2.z   <<endl;
cout<<N3.x  << "," << N3.y   << "," <<N3.z   <<endl;
cout<<C4.x  << "," << C4.y   << "," <<C4.z   <<endl;
cout<<N4.x  << "," << N4.y   << "," <<N4.z   <<endl;
cout<<C5.x  << "," << C5.y   << "," <<C5.z   <<endl;
cout<<H5.x  << "," << H5.y   << "," <<H5.z   <<endl;
cout<<H6.x  << "," << H6.y   << "," <<H6.z   <<endl;
cout<<H41.x << "," << H41.y  << "," <<H41.z  <<endl;
cout<<H42.x << "," << H42.y  << "," <<H42.z  <<endl;
cout<<LP3.x << "," << LP3.y  << "," <<LP3.z  <<endl;
cout<<LP21.x<< "," << LP21.y << "," <<LP21.z <<endl;
cout<<LP22.x<< "," << LP22.y << "," <<LP22.z <<endl;
}     
}
