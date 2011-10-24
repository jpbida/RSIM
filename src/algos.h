#ifndef _ALGO
#define _ALGO
#include <vector>
#include <string>
#include <cstring>

int nextPath(std::vector<int> * path,int n1,int n2);
void allPath(std::vector<int> * out,int max,int upper);

int nextPick(std::vector<int> * vect);
int makeSmallSide(std::vector<int> * small_side,std::string seq,int pos1,int loop_type,int asym_num);
float scoreSmallSide(std::vector<int> * small_side,std::string seq,int pos1, int loop_type,int s_side);
std::string povBaseCU(std::vector<double> points,std::string color);
std::string povBaseAG(std::vector<double> points,std::string color);
std::string povHeader();
std::string povCamera(double loc_x,double loc_y,double loc_z,double look_x,double look_y,double look_z);
std::string povTriangle(std::vector<double> points,std::string color);
std::string povSphere(std::vector<double> points,double radius,std::string color);
std::string povCylinder(std::vector<double> points,double radius,std::string color);
int forward(std::string s, int j);
int backward(std::string s, int j);

#endif
