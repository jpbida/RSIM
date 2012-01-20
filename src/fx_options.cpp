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
#include <string>
#include "stdlib.h"
#include "stdio.h"
#include <cstring>
#include <istream>
#include <sstream>
#include <vector>
#include <iostream>
#include "options.h"
#include "parse.h"
#include "ConfigFile.h"
Options::Options(int argc, char **argv){
using namespace std;
//Currently there are several types of arguments used during simulations
//Score optimization
//Randomization + Simulation
//Degenerate Clustered Fragment Assembly
	char *p0value = NULL;	
	char *p1value = NULL;	
	char *p2value = NULL;	
	char *p3value = NULL;	
	char *p4value = NULL;	
	char *p5value = NULL;	
	char *p6value = NULL;	
	char *p7value = NULL;	
	char *p8value = NULL;	
	char *p9value = NULL;	
	char *p10value = NULL;	
	char *p11value = NULL;	
	char *p12value = NULL;	
	char *p13value = NULL;	
	char *p14value = NULL;	
	char *p15value = NULL;	
	char *p16value = NULL;	
	char *p17value = NULL;	
	char *p18value = NULL;	
	char *p19value = NULL;	
	char *p20value = NULL;	
	int helpflag=0;
	char * pend;
	int c;
      	opterr = 0;
  while ((c = getopt (argc, argv, "l:k:i:f:m:c:r:n:p:g:s:b:a:d:h:v:j:o:e:")) != -1)
    switch (c)
      {
      case 'k': 
	      p20value=optarg;	
	break;
      case 'l': 
	      p19value=optarg;	
	break;
      case 'i':             // Configuration File that contains all of the necessary command line options
	      p18value=optarg;
      case 'e':             // fragment database base coordinates
	      p17value=optarg;
	break;
      case 'o':             // fragment database base coordinates
	      p16value=optarg;
	break;
      case 'j':             // fragment database base coordinates
	      p15value=optarg;
	break;
      case 'v':             // fragment database base coordinates
	      p14value=optarg;
	break;
      case 'a':             // fragment database base coordinates
	      p0value=optarg;
	break;
      case 'f':             // fragment database
	      p1value=optarg;
	break;
      case 'p':             // prefix for output files
	      p2value=optarg;
	break;
      case 'm':             // master structure
	      p3value=optarg;
	break;
      case 'c':             // starting structure 
	      p4value=optarg;
	break;
      case 'n':             // Number of generations  
	      p5value=optarg;
	break;
      case 'g':             // movement group 
	      p7value=optarg;
	break;
      case 's':             // starting position in initial structure file 
	      p8value=optarg;
	break;
      case 'd':             // file containing doublet centers 
	      p9value=optarg;
	break;
      case 'r':             // Number of times to randomize the structure before reducing 
	      p10value=optarg;
	break;
      case 'h':          // help print help message and abort 
		helpflag=1;	
	break;
      case 'x':          // help print help message and abort 
		p11value=optarg;	
	break;
      case 'w':          // help print help message and abort 
		p12value=optarg;	
	break;
      case 'b':          // file contianing coordinates for backbone atoms 
		p13value=optarg;	
	break;
	
      case '?':
        if (isprint (optopt))
{
helpflag=1;

//Error message for simulation binary 
          fprintf (stderr, "Unknown option `-%c'.\n-f fragment database file name\n-c initial starting structure\n-m master structure\n-r number of times to randomize\n-p prefix for output files\n-n number of structures to run\n-s starting position in file \n-g movement group \n-a Base atom fragment database\n-bFull atom fragments\n-d Doublet library params\nExample Call:\n\n ./randmz -f frag_db -a frag_nt -d doublet.lib -g 2 -s 0 -n300 -c p50.struct -m p50_f.struct -p simout -r10 -bfrag_at", optopt);
 }       else{
          fprintf (stderr,
                   "Unknown option character `\\x%x'.\n",
                   optopt);
helpflag=1;
}
      default:
        abort ();
      }
if(helpflag==1){
printf("Closed Moves RNA Monte-Carlo Simulator\n");
}

if(p18value!=NULL)
	{
//Use configuration file to read in values
//Command line values will overwrite if they are given 
string cfile(p18value);
ConfigFile config(cfile);
	string groups;
	config.readInto(groups,"structure");
	printf("Secondary Structure: %s\n",groups.c_str());
//Secondary structure constraints
string gtmp;
int slen=groups.length();
//Remove spaces
string::size_type pos = 0;
bool spacesLeft = true;
while( spacesLeft )
{
pos = groups.find(" ");
if( pos != string::npos )
groups.erase( pos, 1 );
else
spacesLeft = false;
}


for(int j = 0; j< slen; j++){
sec_structure.push_back(-1);
}

for(int j = 0; j< slen; j++){
string pos=groups.substr(j,1);

if(pos.compare("(")==0){
int stop=1;
for(int w=j+1; w<slen; w++){
string pos2=groups.substr(w,1);
if(pos2.compare("(")==0){stop++;}
if(pos2.compare(")")==0){stop=stop-1;}
if(stop==0){
sec_structure[j]=w;
sec_structure[w]=j;
break;
}
}
}
}
printf("\nContraint Length: %d\n", slen);
for(int i =0; i<sec_structure.size(); i++){
//printf("pos: %d - %d\n",i,sec_structure[i]);
}	

//Database Paths
config.readInto(loop_db,"loop_db");
config.readInto(atoms_file,"frag_nt");
config.readInto(gen_confs,"gen_confs");
config.readInto(batoms_file,"frag_at");
config.readInto(fragment_file,"frag_db");
config.readInto(db_prefix,"lib_files");
config.readInto(doublet_params,"doublet_params");
config.readInto(bin_dir,"bin_dir");

//Reading in structure information
config.readInto(sequence,"sequence");
config.readInto(master_struct,"native");
config.readInto(pdb_ids,"pdb_ids");
config.readInto(starting_struct,"starting_pdb");
config.readInto(modelfile,"modelfile");
config.readInto(edge_file,"edge_file");
config.readInto(struct_pos,"model_number");
config.readInto(con_start,"con_start");

//Simulation Parameters
config.readInto(p_o_lower,"p_o_lower");
config.readInto(p_o_upper,"p_o_upper");
config.readInto(iteration_num,"iteration_num");
config.readInto(write_steps,"write_steps");
config.readInto(max_leaves,"max_leaves");
config.readInto(max_clust,"max_clust");
config.readInto(seed,"seed");
config.readInto(constraint_id,"constraint_id");
config.readInto(randomize_num,"randomize_num");
config.readInto(total_simulations,"total_simulations");
config.readInto(movegroup,"movement_group");
config.readInto(debug,"debug");
config.readInto(resolution,"resolution");
config.readInto(small_side,"small_side");
config.readInto(loop_pos,"loop_pos");
config.readInto(loop_type,"loop_type");
config.readInto(asymmetry,"asymmetry");

//output files
config.readInto(prefix,"output_prefix");
config.readInto(input_file,"input_file");
config.readInto(conf_space,"conf_space");
config.readInto(pdb_id,"pdb_id");
out_simname    = prefix;
out_struct     = prefix; 
out_struct.append(".struct"); 
out_scores     = prefix; 
out_scores.append(".scores"); 
out_nnodes     = prefix; 
out_nnodes.append(".nnodes"); 
out_hydro      = prefix; 
out_hydro.append("_hydro.struct"); 
out_doublet    = prefix;
out_doublet.append("_doublet.struct");
out_packing    = prefix;
out_packing.append("_packing.struct");
out_overall    = prefix;
out_overall.append("_overall.struct");
}
	else{
		printf("Configuration file not provided with: -i file_name option.\n  Using default options with command line arguments\n");
//Create file name for nodes that need to be analyzed
//Set Default Values
//Default file locations
//printf("\n\nApplying default values if necessary\n");
if(p20value!=NULL){
string seqtmp(p20value);
sequence=seqtmp;
}
if(p19value!=NULL){
string seqtmp(p19value);
loop_db=seqtmp;
}
if(p15value!=NULL)
	{
	string groups=p15value;
string gtmp;
int slen=groups.length();
//Remove spaces
string::size_type pos = 0;
bool spacesLeft = true;

while( spacesLeft )
{
pos = groups.find(" ");
if( pos != string::npos )
groups.erase( pos, 1 );
else
spacesLeft = false;
}


for(int j = 0; j< slen; j++){
sec_structure.push_back(-1);
}

for(int j = 0; j< slen; j++){
string pos=groups.substr(j,1);

if(pos.compare("(")==0){
int stop=1;
for(int w=j+1; w<slen; w++){
string pos2=groups.substr(w,1);
if(pos2.compare("(")==0){stop++;}
if(pos2.compare(")")==0){stop=stop-1;}
if(stop==0){
sec_structure[j]=w;
sec_structure[w]=j;
break;
}
}
}
}



printf("\nContraint Length: %d\n", slen);
for(int i =0; i<sec_structure.size(); i++){
printf("pos: %d - %d\n",i,sec_structure[i]);
}	
	}
	else{
		//printf("Using default full atom database: ./frag_nt\n"); 
	}
if(p0value!=NULL)
	{
		atoms_file=p0value;
	}
	else{
		atoms_file="frag_nt";
	}

if(p13value!=NULL)
	{
		batoms_file=p13value;
	}
	else{
		batoms_file="frag_at";
		//printf("Using default full atom database: ./frag_at\n"); 
	}

if(p1value!=NULL)
	{
		fragment_file=p1value;
	}
	else{
		fragment_file="frag_db";
		//printf("Using default fragment database: ./frag_db\n"); 
	}

if(p9value!=NULL)
	{
		doublet_params=p9value;
	}
	else{
		doublet_params = "doublet.params";		
//printf("Using default doublet centers file: doublet.params\n"); 
	}

if(p2value!=NULL)
	{
		prefix=p2value;
	}
	else{
		prefix="simout";
		//printf("Using default output prefix: simout\n"); 
	}

if(p3value!=NULL)
	{
		master_struct=p3value;	
	}
	else{
		master_struct="native.struct";
		//printf("Using default native structure: native.struct\n"); 
	}

if(p4value!=NULL)
	{
		starting_struct=p4value;	
	}
	else{
		//starting_struct=NULL;
		//printf("Using default input structure: in.struct\n"); 
	}

if(p17value!=NULL)
	{
		db_prefix=p17value;	
	}
	else{
		db_prefix="../libs/";
		printf("Using default library directory ./libs/\n"); 
	}

if(p5value!=NULL)
	{
		iteration_num= (int)strtod(p5value,&pend);
	}
	else{
		iteration_num=1000;
		//printf("Using default number of iterations: 1000;\n"); 
	}

if(p16value!=NULL)
	{
		seed = (unsigned int)(strtod(p16value,&pend));
	}
	else{
		seed=time(NULL);
		//printf("Using default number of iterations: 1000;\n"); 
	}

if(p14value!=NULL)
	{
//Setting debug = 1 prints out lots of debuging information
		debug= int(strtod(p14value,&pend));
	}
	else{
		debug=0;
		//printf("Using default number of iterations: 1000;\n"); 
	}


if(p7value!=NULL)
	{
		movegroup= int(strtod(p7value,&pend));
		movegroup=movegroup-1;
	}
	else{
		//printf("Using default movement group: 1");
	}


if(p8value!=NULL)
	{
		struct_pos= int(strtod(p8value,&pend));
	}
	else{
		struct_pos = 0;		
//printf("Using default initial structure positon: 0;\n\tIf multiple structures are in the starting structures file, the structure position value determines which is used"); 
	}

if(p10value!=NULL)
	{
		randomize_num= int(strtod(p10value,&pend));
	}
	else{
		randomize_num=50;
		//printf("Using default number of randomizations: 50;\n"); 
	}


if(p11value!=NULL)
	{
		resolution=strtofloat(p11value);
	}
	else{
		resolution = 0.5; 
//printf("Using default cluster resolution: 0.5 Ang\n"); 
	}

if(p12value!=NULL)
	{
		weights=p12value;	
	}
	else{
		weights="weights.dat"; 
//printf("Using default scoring weights: weights.dat\n"); 
	}

out_simname    = prefix;
out_struct     = prefix; 
out_struct.append(".struct"); 
out_scores     = prefix; 
out_scores.append(".scores"); 
out_nnodes     = prefix; 
out_nnodes.append(".nnodes"); 
out_hydro      = prefix; 
out_hydro.append("_hydro.struct"); 
out_doublet    = prefix;
out_doublet.append("_doublet.struct");
out_packing    = prefix;
out_packing.append("_packing.struct");
out_overall    = prefix;
out_overall.append("_overall.struct");

}
}

void Options::print(){
//Prints values for all options

//printf("\n\nCOMMAND LINE OPTIONS\n\nbatoms_file\t%s\natoms_file\t%s\nfragment_file\t%s\noutput prefix\t%s\nnative structure\t%s\nstarting structure\t%s\ndoublet library\t%s\nscoring weights\t%s\niterations\t%d\nrandomization cycles\t%d\nstructure pos\t%d\nmovement group\t%d\ncluster resolution\t%3.3f\nstructures input\t%s\nscores output\t%s\nstructure outputs\t%s\n\n",batoms_file,atoms_file,fragment_file,prefix,master_struct,starting_struct,doublet_params,weights,iteration_num,randomize_num,struct_pos,movegroup,resolution,out_struct,out_scores,out_nnodes);



}



