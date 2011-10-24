// molecule.h -- molecule class and functions
//Splits space into 1A cubes and keeps one frag from each cube as the simulation progresses
#ifndef DCLUST_H_
#define DCLUST_H_
#include <vector>
#include <string>

using namespace std;
class Dynaclust
{
private:
public:
Dynaclust(int clusts,float res);
int cur_tot;					//Current number of clusters 
int Nc; 					//Final Number of Clusters  
float resolution;				//Minimum distance any two fragments can be to each other 
std::vector<float> pos_nc;			//Position Coordinates of each cluster
						//Functions that create and return the clusters
int add(float x, float y, float z);		//Adds a positon to the cluster list
						//this performs dynamic cluster and updates everything
						//appropriately
};                                       


#endif                                     

