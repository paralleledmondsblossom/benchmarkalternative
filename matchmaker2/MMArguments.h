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
#include <string>
using namespace std;
class MMArguments{
public:
  int input_type;
  int match_type;
  int initial_matching_type;
  int rep;
  string input_file;
  int blockDim;
  int threadDim;

  void getArguments(int argc, char *argv[]);
  MMArguments();
  MMArguments(int argc, char *argv[]);
  
  void print_usage();
};
