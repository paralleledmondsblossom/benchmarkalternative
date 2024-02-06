/*
---------------------------------------------------------------------
 This file is part of matchmaker framework
 Copyright (c) 2013,
 By:    Mehmet Deveci,
        Kamer Kaya,
        Bora Ucar,
        Umit V. Catalyurek
---------------------------------------------------------------------
 For license info, please see the LICENSE.txt file in
 the main directory.
---------------------------------------------------------------------
*/

#include "MMArguments.h"
#include <sstream>
#include <string>
#include <cstdlib>
#include <iostream>
using namespace std;

void getArgumentValue(
    string &argumentid, string &argumentValue, char *argumentline){
  
  stringstream stream(stringstream::in | stringstream::out);
  stream << argumentline;
  getline(stream, argumentid, '=');
  if (stream.eof()){
    throw argumentline;
  }
  stream >> argumentValue;
}

void MMArguments::getArguments(int argc, char *argv[]){

  this->input_file = "";
  this->input_type = 1;
  this->initial_matching_type = 0;
  this->rep = -1;
  this->match_type = -1;
  this->blockDim=-1;
  this->threadDim=-1;


  try
  {
    for(int i = 1; i < argc; ++i){
      char *argLine = argv[i];
      string argId = "";
      string argVal = "";
      getArgumentValue(argId, argVal, argLine);
      if (argId == "IT"){
        this->input_type = atoi(argVal.c_str());
        if(this->input_type < 0 || this->input_type >= 2){
          throw argLine;
        }
      } 
      else if(argId == "IF"){
        this->input_file = argVal;
      }
      else if(argId == "MT"){
        this->match_type = atoi(argVal.c_str());
      }
      else if(argId == "IMT"){
        this->initial_matching_type = atoi(argVal.c_str());
      }      
      else if(argId == "R"){
        this->rep = atoi(argVal.c_str());
        if (this->rep < 1){
          throw argLine;
        }
      }
      
      else if(argId == "BD"){
        this->blockDim = atoi(argVal.c_str());
        
        if (this->blockDim < 1){
          throw argLine;
        }
      }
      
      else if(argId == "TD"){
        this->threadDim = atoi(argVal.c_str());
        
        if (this->threadDim < 1){
          throw argLine;
        }
      }
      else {
        throw argLine;
      }
    }
  }
  catch(char * s){
    cout << "Error: " << s << endl;
    this->print_usage();
    exit(1);
  }
  
  if (this->input_file == ""){
    this->print_usage();
    exit(1);
  }
}


void MMArguments::print_usage(){
  cout << "mmaker [argument_list]" << endl;
  cout << "Argument list" << endl; 
  cout << "\tIF: path to the input graph file" <<endl;
  cout << "\tIT: input file format. 0-metis file, 1-matrix market (default)" << endl;
  cout << "\tMT: Maximum Matching Type" << endl;  
  cout << "\t\t0: Sequential DFS" << endl;  
  cout << "\t\t1: Sequential BFS" << endl;  
  cout << "\t\t2: Sequential PF" << endl;  
  cout << "\t\t3: Sequential PFP" << endl;  
  cout << "\t\t4: Sequential HK" << endl;  
  cout << "\t\t5: Sequential HK_DW" << endl;  
  cout << "\t\t6: Sequential ABMP" << endl;  
  cout << "\t\t7: Sequential ABMP_BFS" << endl;  
  cout << "\t\t8: GPU - APFB_GPUBFS" << endl;  
  cout << "\t\t9: GPU - APFB_GPUBFS_WR" << endl;  
  cout << "\t\t10: GPU - APsB_GPUBFS" << endl;  
  cout << "\t\t11: GPU - APsB_GPUBFS_WR" << endl;  
        
  cout << "\tIMT: Initial Matching Type" << endl;    cout << "\t\t0: " << endl;  
  cout << "\t\t0: Cheap Matching" << endl;  
  cout << "\t\t1: SK" << endl;  
  cout << "\t\t2: SK_rand" << endl;  
  cout << "\t\t3: mind_cheap" << endl;  
  cout << "\t\t>=4: no initial matching" << endl;  
  
  
  cout << "\tR: Repetition number. Must be greater than 0." << endl;
  cout << "\tBD: Block Dimension (Grid Size). Must be greater than 0." << endl;
  cout << "\tTD: Thread Dimension (Block Size). Must be greater than 0." << endl;
}


MMArguments::MMArguments(){
}
MMArguments::MMArguments(int argc, char *argv[]){
  this->getArguments(argc, argv);
}
