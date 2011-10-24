// Montecarlo fractal kinetics simulator //
#include <istream>
#include <sstream>
#include <vector>
#include "points.h"
#include "fragatoms.h"
#include <iostream>
#include "parse.h"
#include "fragment.h"

using namespace std;
Fragment::Fragment(){
len_1=0;
len_2=0;
len_3=0;
len_4=0;
bond_0=0;
bond_1=0;
bond_2=0;
tor_1=0;
tor_2=0;
}
Fragment::Fragment(string line){
//Parses a line from an input file and loads information into fragment class
string field;
double len_1_a;
double len_2_a;
double len_3_a;
double len_4_a;
double tor_1_a;
double tor_2_a;
double bond_0_a;
double bond_1_a;
double bond_2_a;
int    frag_type_a;
	istringstream iss(line,istringstream::in);	
   		getline(iss,field,'\t');
        len_1_a=strtodouble(field);
                getline(iss,field,'\t');
        len_2_a=strtodouble(field);
                getline(iss,field,'\t');
        len_3_a=strtodouble(field);
                getline(iss,field,'\t');
        len_4_a=strtodouble(field);
                getline(iss,field,'\t');
        tor_1_a=strtodouble(field);
                getline(iss,field,'\t');
        tor_2_a=strtodouble(field);
                getline(iss,field,'\t');
        bond_0_a=strtodouble(field);
                getline(iss,field,'\t');
        bond_1_a=strtodouble(field);
                getline(iss,field,'\t');
        bond_2_a=strtodouble(field);
                getline(iss,field,'\t');
        frag_type_a=strtoint(field);
len_1=len_1_a;
len_2=len_2_a;
len_3=len_3_a;
len_4=len_4_a;
tor_1=tor_1_a;
tor_2=tor_2_a;
bond_0=bond_0_a;
bond_1=bond_1_a;
bond_2=bond_2_a;
frag_type=frag_type_a;
}

void Fragment::print(){
std::printf("Atom torsion,bond angles, and length for fragment type: \%d\nlen1\tlen2\tlen3\tlen4\ttor1\ttor2\tbond0\tbond1\tbond2\n",frag_type);
std::printf("%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\t%3.2f\n",len_1,len_2,len_3,len_4,tor_1,tor_2,bond_0,bond_1,bond_2);
}
