//g++ CustomChannelSetup.cxx -std=c++0x -o CustomChan

#include "tigPZ.h"
#include "grifPZ.h"
#include "tigadc.h"
#include "grifadc.h"
#include "CalibrationParameters.h"
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm> 
#include <vector>   

using namespace std;

int num_options = 3; //number of known options
string custom_options[] = {
  "pz",   // Pole Correction Paramters
  "adc",  // ADC Numbers
  "gain"  // Linear Gains
};
string option_description[] = {
  "Pole Correction Parameters",
  "Fills Host/Collector Pages", 
  "Sets Linear Gains"
};

int search_array(string array[], string search, int len) {
  for (int i = 0; i < len; i++) {
    if (search == array[i]) {
      return i;
    } //if
  } //for
  return -1;
} //search_array

bool CreateCustomChannel(const char * opt, const char * experiment, const char * out = "CustomChan.com") {

  char line[128];
  char var [64];
  char colour[1];
  int gePos[16];
  int j = 0;
 
  // Checks Which Option and Whether TIGRESS or GRIFFIN 
  int option = search_array(custom_options, opt, num_options);
  if ((strcmp(experiment, "tigress") != 0) && (strcmp(experiment, "griffin") != 0)) {
    printf("Experiment Not Recognised, options are 'griffin' or 'tigress' \n");
    return false;
  }

  // Loads Current Array Configuration
  ifstream ge_pos;
  ge_pos.open("gepositions.txt");
  if (ge_pos.is_open()) {
    printf("Loading Ge Positions\n");
    while (!ge_pos.eof()) {
      ge_pos >> gePos[j];
      j++;
    }
  }

  // Definitions
  std::string tempadd, tempmnemonic, tempdig, tempoutsensor;
  double tempgain, tempoffset, tempnonlin, tempdtype;
  int tempchan;
  std::vector < std::string > ADDRESS, MNEMONIC, DIGITIZER;
  std::vector < int > channelnumber;
  std::vector < float > gains, offsets, nonlins;
  ifstream config;
  string conffile = "Conf_File.txt";

  printf("Creating com file: %s\n", out);
  ofstream outfile;
  outfile.open(out);
  if (!outfile.is_open()) {
    printf("Output not opened\n");
  }

  switch (option) {
  case 0:   // Sets Polecorrections
    
    for (int i = 0; i < 64; i++) {
      int pNum = i / 4 + 1;
      int cryNum = i % 4;
      float pZero = 0;
      if (cryNum == 0) sprintf(colour, "B");
      else if (cryNum == 1) sprintf(colour, "G");
      else if (cryNum == 2) sprintf(colour, "R");
      else if (cryNum == 3) sprintf(colour, "W");
      if (strcmp(experiment, "tigress") == 0) {
        pZero = tigpZero[gePos[pNum - 1] - 1][cryNum];
        sprintf(var, "TIG%2.2i%sN00A", pNum, colour);
      } else if (strcmp(experiment, "griffin") == 0) {
        pZero = grifpZero[gePos[pNum - 1] - 1][cryNum];
        sprintf(var, "GRG%2.2i%sN00A", pNum, colour);
      }
      sprintf(line, "create float \"/DAQ/params/grif16/custom/%s/p_polec1\"",
        var);
        outfile << line << "\n";
      sprintf(line, "set \"/DAQ/params/grif16/custom/%s/p_polec1\" '%1.1f'",
        var, pZero);
      if (gePos[pNum - 1] != 0) outfile << line << "\n";
    }
    break;
  case 1:   // Fills Collector Pages
   
    int num_collector;
    if (strcmp(experiment, "tigress") == 0) {
      num_collector = 4;
    }
    if (strcmp(experiment, "griffin") == 0) {
      num_collector = 3;
    }
    for (int j = 0; j < num_collector; j++) {
      for (int i = 0; i < 16; i++) {
        if (strcmp(experiment, "tigress") == 0 && tigadcNum[j][i] != 0) sprintf(var, "grifadc%i.triumf.ca", tigadcNum[j][i]);
        else if (strcmp(experiment, "tigress") == 0 && tigadcNum[j][i] == 0) sprintf(var, "");
        else if (strcmp(experiment, "griffin") == 0 && grifadcNum[j][i] != 0) sprintf(var, "grifadc%i.triumf.ca", grifadcNum[j][i]);
        else if (strcmp(experiment, "griffin") == 0 && grifadcNum[j][i] == 0) sprintf(var, "");
        sprintf(line, "set \"/DAQ/hosts/collector0x%i/digitizers[%i]\" '%s'", j, i,
          var);
        outfile << line << "\n";
      }
    }
    break;
  case 2:   // Sets Linear Gains 
    
    // Loads Configuration File assuming default name Conf_File.txt
    printf("Loading Configuration File\n");
    config.open(conffile.c_str());
    if (config.is_open()) {
      printf("Configuration file opened\n");
      while (!config.eof()) {
        config >> tempchan >> tempadd >> tempmnemonic >> tempgain >> tempoffset >> tempnonlin >> tempdig;
        channelnumber.push_back(tempchan);
        ADDRESS.push_back(tempadd);
        MNEMONIC.push_back(tempmnemonic);
        gains.push_back(tempgain);
        offsets.push_back(tempoffset);
        nonlins.push_back(tempnonlin);
        DIGITIZER.push_back(tempdig);
        j++;
      }
    } else {
      printf("Conf file %s failed to open, abort!\n", conffile.c_str());
      return false;
    }

    for (int i = 0; i < 64; i++) {
      int pNum = i / 4 + 1;
      int cryNum = i % 4;
      float pZero = 0;
      if (cryNum == 0) sprintf(colour, "B");
      else if (cryNum == 1) sprintf(colour, "G");
      else if (cryNum == 2) sprintf(colour, "R");
      else if (cryNum == 3) sprintf(colour, "W");
      if (strcmp(experiment, "tigress") == 0) {
        sprintf(var, "TIG%2.2i%sN00A", pNum, colour);
      } else if (strcmp(experiment, "griffin") == 0) {
        sprintf(var, "GRG%2.2i%sN00A", pNum, colour);
      }

      for (int iii = 0; iii < MNEMONIC.size(); iii++) {
        tempmnemonic = MNEMONIC.at(iii);
        if (tempmnemonic.compare(var) == 0) {
          sprintf(line, "set \"/DAQ/MSC/gain[%i]\" '%f'", iii, gain[i]);
          outfile << line << "\n";
          sprintf(line, "set \"/DAQ/MSC/offset[%i]\" '%f'", iii, offset[i]);
          outfile << line << "\n";
          sprintf(line, "set \"/DAQ/MSC/quadratic[%i]\" '%e'", iii, non_lin[i]);
          outfile << line << "\n";
        }
      }
    }
    break;
  default:
    printf("Option not recognised\nKnown Options are:\n");
    for (int i = 0; i < num_options; i++) printf("%s:\t%s\n", custom_options[i].c_str(), option_description[i].c_str());
    return false;
  }

  outfile.close();
  return true;
}

int main(int argc, char * * argv) {

  if (argc < 3) {
    printf("Too few inputs, give custom option and experiment name\nRecognised options:\n");
    for(int i = 0; i < num_options; i++) printf("%s:\t%s\n",custom_options[i].c_str(), option_description[i].c_str());
    printf("Recognised Experiment:\ntigress\ngriffin\n");
    return 0;
  }

  if (argc == 3) {
    printf("Creating Custom Channel Command File: CustomChan.com\n");
    CreateCustomChannel(argv[1], argv[2]);
  } 

  if (argc == 4) {
    printf("Creating Custom Channel Command File: %s",argv[3]);
    CreateCustomChannel(argv[1], argv[2], argv[3]);
  } 

  if (argc > 4) {
    printf("Too many Arguments! Max 3, Custom Option, Experiment Name(tigress/griffin), Optional Command File Name (Default CustomChan.com)");
  } 
  return 0;

}
