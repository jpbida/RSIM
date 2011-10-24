#include "algos.h"
#include <cmath>
#include <boost/lexical_cast.hpp>
#include <vector>
#include <string>
#include <cstring>
#include <iostream>
#include <stdio.h>
//Simple Povray Api
std::string povHeader(){
std::string shead("#include \"colors.inc\"\n#include \"shapes.inc\"\nlight_source{<-500,1500,-1500> color White*0.9}           // sun light\n\n// sky -------------------------------------------------------------------\nsky_sphere{ pigment{\ncolor rgb<0.14,0.14,0.56>*3\n                      }\n           } // end of sky_sphere\n\n\n\n#declare G_Texture =      \n         texture { pigment{ color Blue }  \n                   finish { phong 1 reflection 0.00}\n                 } // end of texture \n#declare A_Texture =\n         texture { pigment{ color Red }   \n                   finish { phong 1 reflection 0.00}\n                 } // end of texture \n#declare U_Texture = \n         texture { pigment{ color Green } \n                   finish { phong 1 reflection 0.00}\n                 } // end of texture \n#declare C_Texture =\n         texture { pigment{ color Yellow }\n                   finish { phong 1 reflection 0.00}\n                 } // end of texture \n#declare Cylinder_Texture = Silver;\n");
return shead;
}

std::string povCamera(double loc_x,double loc_y,double loc_z,double look_x,double look_y,double look_z){
	std::string scamera("camera {\n\tlocation <");
        scamera+=boost::lexical_cast<std::string>(loc_x);
	scamera.append(",");
        scamera+=boost::lexical_cast<std::string>(loc_y);
	scamera.append(",");
        scamera+=boost::lexical_cast<std::string>(loc_z);
	scamera.append(">\n\tlook_at <");
        scamera+=boost::lexical_cast<std::string>(look_x);
	scamera.append(",");
        scamera+=boost::lexical_cast<std::string>(look_y);
	scamera.append(",");
        scamera+=boost::lexical_cast<std::string>(look_z);
	scamera.append(">\n}\n");
	scamera.append("light_source{ <");
        scamera+=boost::lexical_cast<std::string>(loc_x);
	scamera.append(",");
        scamera+=boost::lexical_cast<std::string>(loc_y);
	scamera.append(",");
        scamera+=boost::lexical_cast<std::string>(loc_z);
	scamera.append("> color rgb<1,1,1>}\n");
return scamera;
}

std::string povTriangle(std::vector<double> points,std::string color){
//vector is three points
std::string triangle("triangle{\n\t<");
        triangle+=boost::lexical_cast<std::string>(points[0]);triangle.append(",");
        triangle+=boost::lexical_cast<std::string>(points[1]);triangle.append(",");
        triangle+=boost::lexical_cast<std::string>(points[2]);triangle.append(">,<");
        triangle+=boost::lexical_cast<std::string>(points[3]);triangle.append(",");
        triangle+=boost::lexical_cast<std::string>(points[4]);triangle.append(",");
        triangle+=boost::lexical_cast<std::string>(points[5]);triangle.append(">,<");
        triangle+=boost::lexical_cast<std::string>(points[6]);triangle.append(",");
        triangle+=boost::lexical_cast<std::string>(points[7]);triangle.append(",");
        triangle+=boost::lexical_cast<std::string>(points[8]);triangle.append(">\n");
triangle.append("texture {\n\t");
triangle.append(color);
triangle.append(" }\n\t}\n");
return triangle;
}

std::string povBaseCU(std::vector<double> points,std::string color){
//tri1
std::vector<double> tri1;
tri1.push_back(points[0]);
tri1.push_back(points[1]);
tri1.push_back(points[2]);

tri1.push_back(points[3]);
tri1.push_back(points[4]);
tri1.push_back(points[5]);

tri1.push_back(points[6]);
tri1.push_back(points[7]);
tri1.push_back(points[8]);
std::string tri1_s=povTriangle(tri1,color);

std::vector<double> tri2;
tri2.push_back(points[6]);
tri2.push_back(points[7]);
tri2.push_back(points[8]);

tri2.push_back(points[9]);
tri2.push_back(points[10]);
tri2.push_back(points[11]);

tri2.push_back(points[0]);
tri2.push_back(points[1]);
tri2.push_back(points[2]);
std::string tri2_s=povTriangle(tri2,color);

std::vector<double> tri3;
tri3.push_back(points[9]);
tri3.push_back(points[10]);
tri3.push_back(points[11]);

tri3.push_back(points[15]);
tri3.push_back(points[16]);
tri3.push_back(points[17]);

tri3.push_back(points[0]);
tri3.push_back(points[1]);
tri3.push_back(points[2]);
std::string tri3_s=povTriangle(tri3,color);

std::vector<double> tri4;
tri4.push_back(points[9]);
tri4.push_back(points[10]);
tri4.push_back(points[11]);

tri4.push_back(points[12]);
tri4.push_back(points[13]);
tri4.push_back(points[14]);

tri4.push_back(points[15]);
tri4.push_back(points[16]);
tri4.push_back(points[17]);

std::string tri4_s=povTriangle(tri4,color);
std::vector<double> ends1,ends2,ends3,ends4,ends5,ends6;
ends1.push_back(points[(1-1)*3+0]);
ends1.push_back(points[(1-1)*3+1]);
ends1.push_back(points[(1-1)*3+2]);

ends1.push_back(points[(2-1)*3+0]);
ends1.push_back(points[(2-1)*3+1]);
ends1.push_back(points[(2-1)*3+2]);

ends2.push_back(points[(2-1)*3+0]);
ends2.push_back(points[(2-1)*3+1]);
ends2.push_back(points[(2-1)*3+2]);

ends2.push_back(points[(3-1)*3+0]);
ends2.push_back(points[(3-1)*3+1]);
ends2.push_back(points[(3-1)*3+2]);

ends3.push_back(points[(3-1)*3+0]);
ends3.push_back(points[(3-1)*3+1]);
ends3.push_back(points[(3-1)*3+2]);

ends3.push_back(points[(4-1)*3+0]);
ends3.push_back(points[(4-1)*3+1]);
ends3.push_back(points[(4-1)*3+2]);

ends4.push_back(points[(4-1)*3+0]);
ends4.push_back(points[(4-1)*3+1]);
ends4.push_back(points[(4-1)*3+2]);

ends4.push_back(points[(5-1)*3+0]);
ends4.push_back(points[(5-1)*3+1]);
ends4.push_back(points[(5-1)*3+2]);

ends5.push_back(points[(5-1)*3+0]);
ends5.push_back(points[(5-1)*3+1]);
ends5.push_back(points[(5-1)*3+2]);

ends5.push_back(points[(6-1)*3+0]);
ends5.push_back(points[(6-1)*3+1]);
ends5.push_back(points[(6-1)*3+2]);

ends6.push_back(points[(6-1)*3+0]);
ends6.push_back(points[(6-1)*3+1]);
ends6.push_back(points[(6-1)*3+2]);

ends6.push_back(points[(1-1)*3+0]);
ends6.push_back(points[(1-1)*3+1]);
ends6.push_back(points[(1-1)*3+2]);

std::string cyl_col("Cylinder_Texture");
std::string line1  = povCylinder(ends1,0.15,cyl_col);
std::string line2  = povCylinder(ends2,0.15,cyl_col);
std::string line3  = povCylinder(ends3,0.15,cyl_col);
std::string line4  = povCylinder(ends4,0.15,cyl_col);
std::string line5  = povCylinder(ends5,0.15,cyl_col);
std::string line6  = povCylinder(ends6,0.15,cyl_col);

std::string tri_out;
tri_out.append(tri1_s);
tri_out.append(tri2_s);
tri_out.append(tri3_s);
tri_out.append(tri4_s);
tri_out.append(line1 );
tri_out.append(line2 );
tri_out.append(line3 );
tri_out.append(line4 );
tri_out.append(line5 );
tri_out.append(line6 );

return tri_out;
}

std::string povBaseAG(std::vector<double> points,std::string color){
std::string cyl_col("Cylinder_Texture");
//tri1
std::vector<double> tri1;
tri1.push_back(points[0]);
tri1.push_back(points[1]);
tri1.push_back(points[2]);

tri1.push_back(points[3]);
tri1.push_back(points[4]);
tri1.push_back(points[5]);



tri1.push_back(points[6]);
tri1.push_back(points[7]);
tri1.push_back(points[8]);
std::string tri1_s=povTriangle(tri1,color);

std::vector<double> tri2;
tri2.push_back(points[6]);
tri2.push_back(points[7]);
tri2.push_back(points[8]);

tri2.push_back(points[9]);
tri2.push_back(points[10]);
tri2.push_back(points[11]);

tri2.push_back(points[0]);
tri2.push_back(points[1]);
tri2.push_back(points[2]);
std::string tri2_s=povTriangle(tri2,color);

std::vector<double> tri3;
tri3.push_back(points[9]);
tri3.push_back(points[10]);
tri3.push_back(points[11]);

tri3.push_back(points[15]);
tri3.push_back(points[16]);
tri3.push_back(points[17]);

tri3.push_back(points[0]);
tri3.push_back(points[1]);
tri3.push_back(points[2]);
std::string tri3_s=povTriangle(tri3,color);

std::vector<double> tri4;
tri4.push_back(points[9]);
tri4.push_back(points[10]);
tri4.push_back(points[11]);

tri4.push_back(points[12]);
tri4.push_back(points[13]);
tri4.push_back(points[14]);

tri4.push_back(points[15]);
tri4.push_back(points[16]);
tri4.push_back(points[17]);
std::string tri4_s=povTriangle(tri4,color);

std::vector<double> tri5;
tri5.push_back(points[0]);
tri5.push_back(points[1]);
tri5.push_back(points[2]);

tri5.push_back(points[(6-1)*3+0]);
tri5.push_back(points[(6-1)*3+1]);
tri5.push_back(points[(6-1)*3+2]);

tri5.push_back(points[(9-1)*3+0]);
tri5.push_back(points[(9-1)*3+1]);
tri5.push_back(points[(9-1)*3+2]);

std::string tri5_s=povTriangle(tri5,color);

std::vector<double> tri6;

tri6.push_back(points[(6-1)*3+0]);
tri6.push_back(points[(6-1)*3+1]);
tri6.push_back(points[(6-1)*3+2]);

tri6.push_back(points[(7-1)*3+0]);
tri6.push_back(points[(7-1)*3+1]);
tri6.push_back(points[(7-1)*3+2]);

tri6.push_back(points[(9-1)*3+0]);
tri6.push_back(points[(9-1)*3+1]);
tri6.push_back(points[(9-1)*3+2]);

std::string tri6_s=povTriangle(tri6,color);

std::vector<double> tri7;

tri7.push_back(points[(7-1)*3+0]);
tri7.push_back(points[(7-1)*3+1]);
tri7.push_back(points[(7-1)*3+2]);

tri7.push_back(points[(8-1)*3+0]);
tri7.push_back(points[(8-1)*3+1]);
tri7.push_back(points[(8-1)*3+2]);

tri7.push_back(points[(9-1)*3+0]);
tri7.push_back(points[(9-1)*3+1]);
tri7.push_back(points[(9-1)*3+2]);

std::string tri7_s=povTriangle(tri7,color);

std::vector<double> ends1,ends2,ends3,ends4,ends5,ends6,ends7,ends8,ends9,ends10;
ends1.push_back(points[(1-1)*3+0]);
ends1.push_back(points[(1-1)*3+1]);
ends1.push_back(points[(1-1)*3+2]);

ends1.push_back(points[(2-1)*3+0]);
ends1.push_back(points[(2-1)*3+1]);
ends1.push_back(points[(2-1)*3+2]);

ends2.push_back(points[(2-1)*3+0]);
ends2.push_back(points[(2-1)*3+1]);
ends2.push_back(points[(2-1)*3+2]);

ends2.push_back(points[(3-1)*3+0]);
ends2.push_back(points[(3-1)*3+1]);
ends2.push_back(points[(3-1)*3+2]);

ends3.push_back(points[(3-1)*3+0]);
ends3.push_back(points[(3-1)*3+1]);
ends3.push_back(points[(3-1)*3+2]);

ends3.push_back(points[(4-1)*3+0]);
ends3.push_back(points[(4-1)*3+1]);
ends3.push_back(points[(4-1)*3+2]);

ends4.push_back(points[(4-1)*3+0]);
ends4.push_back(points[(4-1)*3+1]);
ends4.push_back(points[(4-1)*3+2]);

ends4.push_back(points[(5-1)*3+0]);
ends4.push_back(points[(5-1)*3+1]);
ends4.push_back(points[(5-1)*3+2]);

ends5.push_back(points[(5-1)*3+0]);
ends5.push_back(points[(5-1)*3+1]);
ends5.push_back(points[(5-1)*3+2]);

ends5.push_back(points[(6-1)*3+0]);
ends5.push_back(points[(6-1)*3+1]);
ends5.push_back(points[(6-1)*3+2]);

ends6.push_back(points[(6-1)*3+0]);
ends6.push_back(points[(6-1)*3+1]);
ends6.push_back(points[(6-1)*3+2]);

ends6.push_back(points[(1-1)*3+0]);
ends6.push_back(points[(1-1)*3+1]);
ends6.push_back(points[(1-1)*3+2]);

ends7.push_back(points[(6-1)*3+0]);
ends7.push_back(points[(6-1)*3+1]);
ends7.push_back(points[(6-1)*3+2]);

ends7.push_back(points[(7-1)*3+0]);
ends7.push_back(points[(7-1)*3+1]);
ends7.push_back(points[(7-1)*3+2]);

ends8.push_back(points[(7-1)*3+0]);
ends8.push_back(points[(7-1)*3+1]);
ends8.push_back(points[(7-1)*3+2]);

ends8.push_back(points[(8-1)*3+0]);
ends8.push_back(points[(8-1)*3+1]);
ends8.push_back(points[(8-1)*3+2]);

ends9.push_back(points[(8-1)*3+0]);
ends9.push_back(points[(8-1)*3+1]);
ends9.push_back(points[(8-1)*3+2]);

ends9.push_back(points[(9-1)*3+0]);
ends9.push_back(points[(9-1)*3+1]);
ends9.push_back(points[(9-1)*3+2]);

ends10.push_back(points[(9-1)*3+0]);
ends10.push_back(points[(9-1)*3+1]);
ends10.push_back(points[(9-1)*3+2]);

ends10.push_back(points[(1-1)*3+0]);
ends10.push_back(points[(1-1)*3+1]);
ends10.push_back(points[(1-1)*3+2]);


std::string line1  = povCylinder(ends1,0.15,cyl_col);
std::string line2  = povCylinder(ends2,0.15,cyl_col);
std::string line3  = povCylinder(ends3,0.15,cyl_col);
std::string line4  = povCylinder(ends4,0.15,cyl_col);
std::string line5  = povCylinder(ends5,0.15,cyl_col);
std::string line6  = povCylinder(ends6,0.15,cyl_col);
std::string line7  = povCylinder(ends7,0.15,cyl_col);
std::string line8  = povCylinder(ends8,0.15,cyl_col);
std::string line9  = povCylinder(ends9,0.15,cyl_col);
std::string line10 = povCylinder(ends10,0.15,cyl_col);

std::string tri_out;
tri_out.append(tri1_s);
tri_out.append(tri2_s);
tri_out.append(tri3_s);
tri_out.append(tri4_s);
tri_out.append(tri5_s);
tri_out.append(tri6_s);
tri_out.append(tri7_s);
tri_out.append(tri7_s);
tri_out.append(line1 );
tri_out.append(line2 );
tri_out.append(line3 );
tri_out.append(line4 );
tri_out.append(line5 );
tri_out.append(line6 );
tri_out.append(line7 );
tri_out.append(line8 );
tri_out.append(line9 );
tri_out.append(line10);



return tri_out;
}


std::string povSphere(std::vector<double> points,double radius,std::string color){
//vector is three points
std::string triangle("sphere\n{\n\t<");
        triangle+=boost::lexical_cast<std::string>(points[0]);triangle.append(",");
        triangle+=boost::lexical_cast<std::string>(points[1]);triangle.append(",");
        triangle+=boost::lexical_cast<std::string>(points[2]);triangle.append(">,");
        triangle+=boost::lexical_cast<std::string>(radius);triangle.append("\n");
triangle.append("texture {\n\t");
triangle.append(color);
triangle.append(" }\n\t}\n");
return triangle;
}

std::string povCylinder(std::vector<double> points,double radius,std::string color){
//vector is three points
std::string triangle("cylinder\n{\n\t<");
        triangle+=boost::lexical_cast<std::string>(points[0]);triangle.append(",");
        triangle+=boost::lexical_cast<std::string>(points[1]);triangle.append(",");
        triangle+=boost::lexical_cast<std::string>(points[2]);triangle.append(">,<");
        triangle+=boost::lexical_cast<std::string>(points[3]);triangle.append(",");
        triangle+=boost::lexical_cast<std::string>(points[4]);triangle.append(",");
        triangle+=boost::lexical_cast<std::string>(points[5]);triangle.append(">,");
        triangle+=boost::lexical_cast<std::string>(radius);triangle.append("\n");
triangle.append("texture {\n\tpigment {color ");
triangle.append(color);
triangle.append(" }\n\t}\n\n}\n");
return triangle;
}


float scoreSmallSide(std::vector<int> * small_side,std::string seq,int pos1, int loop_type,int s_side){
//Returns the score for W-C secondary structure created by the asymmetry
std::vector< std::vector< std::pair<int,int> > > cons;
//Check that lenght small_side.size + small_size==1 = e2-s2
float score =0;
int m=0;
int s1=pos1;
int s2=0;
if(loop_type==0){
s2=pos1+7;
}else{
if(loop_type==1){
s2=pos1+8;
}else{
s2=pos1;
}
}

if(s_side==1){
for(int j=0; j<(small_side->size()); j++){
//printf("m:j %d:%d\n",m,j);
std::vector< std::pair<int,int> > bpcons;
	if(small_side->at(j)==0){
		std::pair<int,int> scon((s1-j),(s2+m+0));	
		bpcons.push_back(scon);
	if
		(
	seq.at(s1-j)=='G' && seq.at(s2+m)=='C' ||
	seq.at(s1-j)=='C' && seq.at(s2+m)=='G'
		){score=score+3;}
	if
		(
	seq.at(s1-j)=='U' && seq.at(s2+m)=='G' ||
	seq.at(s1-j)=='U' && seq.at(s2+m)=='A' ||
	seq.at(s1-j)=='A' && seq.at(s2+m)=='U' ||
	seq.at(s1-j)=='G' && seq.at(s2+m)=='U'
		){score=score+2;}
		m++;
		}else{
		std::pair<int,int> scon1((s1-j),(s2+m));	
		bpcons.push_back(scon1);
	float sc1=0;
	float sc2=0;
	if
		(
	seq.at(s1-j)=='G' && seq.at(s2+m)=='C' ||
	seq.at(s1-j)=='C' && seq.at(s2+m)=='G'
		){sc1=3;}
	if
		(
	seq.at(s1-j)=='U' && seq.at(s2+m)=='G' ||
	seq.at(s1-j)=='U' && seq.at(s2+m)=='A' ||
	seq.at(s1-j)=='A' && seq.at(s2+m)=='U' ||
	seq.at(s1-j)=='G' && seq.at(s2+m)=='U'
		){sc1=2;}
		m++;
		std::pair<int,int> scon2((s1-j),(s2+m));	
		bpcons.push_back(scon2);
	if
		(
	seq.at(s1-j)=='G' && seq.at(s2+m)=='C' ||
	seq.at(s1-j)=='C' && seq.at(s2+m)=='G'
		){sc2=3;}
	if
		(
	seq.at(s1-j)=='U' && seq.at(s2+m)=='G' ||
	seq.at(s1-j)=='U' && seq.at(s2+m)=='A' ||
	seq.at(s1-j)=='A' && seq.at(s2+m)=='U' ||
	seq.at(s1-j)=='G' && seq.at(s2+m)=='U'
		){sc2=2;}
	if(sc1 > sc2){score=score+sc1;}else{score=score+sc2;}		
		m++;
		}
cons.push_back(bpcons);
}
}else{
if(s_side==2){
for(int j=0; j<(small_side->size()); j++){
std::vector< std::pair<int,int> > bpcons;
	if(small_side->at(j)==0){
	if
		(
	seq.at(s1-m)=='G' && seq.at(s2+j)=='C' ||
	seq.at(s1-m)=='C' && seq.at(s2+j)=='G'
		){score=score+3;}
	if
		(
	seq.at(s1-m)=='U' && seq.at(s2+j)=='G' ||
	seq.at(s1-m)=='U' && seq.at(s2+j)=='A' ||
	seq.at(s1-m)=='A' && seq.at(s2+j)=='U' ||
	seq.at(s1-m)=='G' && seq.at(s2+j)=='U'
		){score=score+2;}
		std::pair<int,int> scon((s1-m),(s2+j));	
		bpcons.push_back(scon);
		m++;
		}else{
	float sc1=0;
	float sc2=0;
	if
		(
	seq.at(s1-m)=='G' && seq.at(s2+j)=='C' ||
	seq.at(s1-m)=='C' && seq.at(s2+j)=='G'
		){sc1=3;}
	if
		(
	seq.at(s1-m)=='U' && seq.at(s2+j)=='G' ||
	seq.at(s1-m)=='U' && seq.at(s2+j)=='A' ||
	seq.at(s1-m)=='A' && seq.at(s2+j)=='U' ||
	seq.at(s1-m)=='G' && seq.at(s2+j)=='U'
		){sc1=2;}
		std::pair<int,int> scon1((s1-m),(s2+j));	
		bpcons.push_back(scon1);
		m++;
		std::pair<int,int> scon2((s1-m),(s2+j));	
	if
		(
	seq.at(s1-m)=='G' && seq.at(s2+j)=='C' ||
	seq.at(s1-m)=='C' && seq.at(s2+j)=='G'
		){sc2=3;}
	if
		(
	seq.at(s1-m)=='U' && seq.at(s2+j)=='G' ||
	seq.at(s1-m)=='U' && seq.at(s2+j)=='A' ||
	seq.at(s1-m)=='A' && seq.at(s2+j)=='U' ||
	seq.at(s1-m)=='G' && seq.at(s2+j)=='U'
		){sc2=2;}
		bpcons.push_back(scon2);
		m++;
	if(sc1 > sc2){score=score+sc1;}else{score=score+sc2;}		
}
cons.push_back(bpcons);
}
}else{
score = -1;
}




}

//Print out constraints for debugging purpose
//printf("Score: %8.3f\n",score);
for(int j=0; j<cons.size(); j++){
	//printf("pos %d\n",j);
	for(int m=0; m<cons[j].size(); m++){
	//	printf("\t%c:%d-%c:%d\n",seq.at(cons[j][m].first),cons[j][m].first,seq.at(cons[j][m].second),cons[j][m].second);
	}
}

return score;

}

int makeSmallSide(std::vector<int> * small_side,std::string seq,int pos1,int loop_type,int asym_num){

//Small side is a vector that contains the mapping from the small side to the large side
//The sequence is used to get the overall length
//pos1 is the position where the loop starts
//loop_type is even=0, odd=1
//asym_num is the desired asymmetry number

//The program returns:
//	 1 if small side is pos 1,
//	 2 if small side is pos 2,
//	 0 if the asym_num is too high or the ratio is too large


//returns the score for the asymmetric sequence
int len_tot=seq.length();
int s_side=0;
int len1 = pos1+1;
int len2 = 0;
if(loop_type==0){
len2=len_tot-(len1+6);
}else{
if(loop_type==1){
len2=len_tot-(len1+7);
}else{
printf("No loop type given: 0=even, 1=odd\n");
return 0;
}
}
 // Add size for even loops
//Check ratio
float iratio = (float)((float)len1/(float)len2);
//printf("ratio: %8.3f\n",iratio);
if(iratio <= 1.5 && iratio >= 0.67){
//Create small_side vector std::vector<int> small_side;
int min_len=0;
int max_len=0;
if(len1 <= len2){
	min_len=len1;
	max_len=len2;
	s_side=1;
}else{
min_len=len2;
max_len=len1;
s_side=2;
}

small_side->clear();
for(int i=0; i<min_len; i++){
	if(i < (min_len-(max_len-min_len))){
		small_side->push_back(0);
	}else{
		small_side->push_back(1);
	}
}

//Iterate through the asy_num
int go=0;
for(int i=0; i<asym_num; i++){
go=nextPick(small_side);
if(go==0){return 0; break;}
}

/*
printf("small_side: ");
for(int g=0; g<small_side->size(); g++){
printf("%d",small_side->at(g));
}
*/
}else{
s_side=0;
}
return s_side;
}


void allPath(std::vector<int> * out,int max,int upper){
//Base 3 vector of length min(n1,n2)
//Add 1 to the current vector
std::vector<int> path;
for(int m=0; m<max; m++){path.push_back(0);}

int size = max-1;
int good=0;
int max_tot=0;
int count=0;
while(count < upper){
if(count % 10000==0){
for(int m=0; m<out->size(); m++){
printf("%d,%d\n",m,out->at(m));
}
}
max_tot=1;
int rem=1;
for(int i=0; i<max; i++){
if((path[size-i]+rem)<3){
path[size-i]=(path[size-i]+rem) % 3;
rem=0;
}else{
path[size-i]=0;
rem=1;
}
}

//Check to see what the path produces
//moves 0 = (1,1), 1=(0,1)(1,1), 2 = (1,0)(1,1)
//if it reaches (n1,n1) before getting to its length
//Set the vector to the maximum value that can be obtained 
//before getting to that position

int val1=0;
int val2=0;
int cur_val=0;
for(int i=0; i<path.size(); i++){
//Foreach step count the occurances 
if(val1 <= val2 && cur_val == count){out->at(val1)=out->at(val1)++;}
if(path[i]==0){val1++; val2++; max_tot=0;}
if(path[i]==1){val1++; val1++; max_tot=0; val2++;cur_val=pow(3,i)+cur_val;}
if(path[i]==2){val1++; val2++; val2++;cur_val=2*pow(3,i)+cur_val;}
}
count++;
}


}


int nextPath(std::vector<int> * path, int n1,int n2){
//Base 3 vector of length min(n1,n2)
//Add 1 to the current vector
int size = path->size()-1;
int good=0;
int max_tot=0;
while(good==0 && max_tot==0){
max_tot=1;
int rem=1;
for(int i=0; i<path->size(); i++){
if((path->at(size-i)+rem)<3){
path->at(size-i)=(path->at(size-i)+rem) % 3;
rem=0;
}else{
path->at(size-i)=0;
rem=1;
}
}

//Check to see what the path produces
//moves 0 = (1,1), 1=(0,1)(1,1), 2 = (1,0)(1,1)
//if it reaches (n1,n1) before getting to its length
//Set the vector to the maximum value that can be obtained 
//before getting to that position

int val1=0;
int val2=0;
for(int i=0; i<path->size(); i++){
if(val1==n1 && val2==n2){
path->at(i)=2;
good=1;
}else{
if(path->at(i)==0){val1++; val2++; max_tot=0;}
if(path->at(i)==1){val1++; val1++; max_tot=0; val2++;}
if(path->at(i)==2){val1++; val2++; val2++;}
}
}
if(val1==n1 && val2==n2){good=1;}
}

for(int i=0; i<path->size(); i++){
//printf("%d",path->at(i));
}
//printf("\n");

int out=0;
if(max_tot==1){
out=1;
}
return out;
}

int nextPick(std::vector<int> * vect ){
//Returns the next vector in a series that picks all k pairs from n numbers
//Get positions of each 1 
int last=0;
std::vector<int> positions;
for(int i=0; i<vect->size(); i++){
if(vect->at(i)==1){
positions.push_back(i);
}
}

if(positions.size() >0){
last = positions.back();
//printf("Last:%d\n",last);
if(last > (positions.size()-1)){
int minp=0;
int movepos=-1;
for(int i=0; i<(positions.size()-1); i++){
	if(positions[i]>minp){
	movepos=i;
	break;
	}else{
minp++;
}
}

if(movepos==-1){
//move the last position forward one
for(int j=0; j<positions.size(); j++){
vect->at(positions[j])=0;
}
//Move all others infront of the last
for(int j=0; j<positions.size(); j++){
vect->at((positions.back()-1-j))=1;
}

}else{
for(int j=0; j<=movepos; j++){
vect->at(positions[j])=0;
}

for(int j=0; j<=movepos; j++){
vect->at(positions[movepos]-1-j)=1;
}

}
return 1;
}else{
//No more in the list
return 0;
}
}else{
return 0;
}
}

//Functions for manipulating secondary structure
int backward(std::string s, int i){
//i is a postion in the string s
//s is a bracket representation of a secondary structure
//this function finds the matching bracket in the in the string s
int k=1;
int go=1;
int out=-1;
if(s.compare(i,1,")")==0){
i--;
for(int j=i; j>=0; j--){
if(go==1){
        if(s.compare(j,1,"(")==0){
        k--;        
	if(k==0){
                        out=j;
                        go=0;
                }
        }else{
                if(s.compare(j,1,")")==0){
                k++;
        }

        }
}
}
}
return out;
}
 
int forward(std::string s, int i){
//i is a postion in the string s
//s is a bracket representation of a secondary structure
//this function finds the matching bracket in the in the string s
int k=1;
int go=1;
int out=-1;
if(s.compare(i,1,"(")==0){
i++;
for(int j=i; j< s.length(); j++){
if(go==1){
        if(s.compare(j,1,")")==0){
        k--;        
	if(k==0){
                        out=j;
                        go=0;
                }
        }else{
                if(s.compare(j,1,"(")==0){
                k=k+1;
        }

        }
}
}
}
return out;
}






