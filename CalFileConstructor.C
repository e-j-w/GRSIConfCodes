// g++ CalFileConstructor.C -std=c++0x -o ConstructCalibrationFile 

#include <fstream>
#include <iostream>
#include <vector>
#include <iomanip> 
#include <stdio.h>
#include <string.h>
using namespace std;

void ConstructCalFile(const char *inp = "Conf_File.txt", const char *out = "CalibrationFile.cal", const char *emma = "NULL"){
  std::string tempadd, tempmnemonic, tempadd2, tempdig;
  double tempgain, tempoffset, tempnonlin, tempdtype;
  int tempchan;

  std::vector < std::string > ADDRESS, MNEMONIC, DIGITIZER;
  std::vector < int > channelnumber;
  std::vector < double > gain, offset, nonlin;
 
  int j = 0;
  int x;
  ifstream infile;
  infile.open(inp);
  if (infile.is_open()) {
    printf("Conf file opened\n");
    while (infile.good()) {
      infile >> tempchan >> tempadd >> tempmnemonic >> tempgain >> tempoffset >> tempnonlin >> tempdig;
      channelnumber.push_back(tempchan);
      ADDRESS.push_back(tempadd);
      MNEMONIC.push_back(tempmnemonic);
      gain.push_back(tempgain);
      offset.push_back(tempoffset);
      nonlin.push_back(tempnonlin);
      DIGITIZER.push_back(tempdig);
      j++;
    }
  } else {
    printf("Conf file %s failed to open, abort!\n", inp);
    return;
  }
  channelnumber.pop_back();
  ADDRESS.pop_back();
  MNEMONIC.pop_back();
  gain.pop_back();
  offset.pop_back();
  nonlin.pop_back();
  DIGITIZER.pop_back();

  // Adds EMMA if a seperate Configuration File is provided, otherwise does nothing. 
  j += 5;
  const char * testval = "NULL";
  if (strcmp(emma, testval) != 0) {
    ifstream inemma;
    inemma.open(emma);
    tempchan = j;
    while (inemma.good()) {
      inemma >> tempadd >> tempmnemonic >> tempgain >> tempoffset >> tempnonlin >> tempdig;
      channelnumber.push_back(tempchan);
      ADDRESS.push_back(tempadd);
      MNEMONIC.push_back(tempmnemonic);
      gain.push_back(tempgain);
      offset.push_back(tempoffset);
      nonlin.push_back(tempnonlin);
      DIGITIZER.push_back(tempdig);
      j++;
      tempchan++;
    }
    channelnumber.pop_back();
    ADDRESS.pop_back();
    MNEMONIC.pop_back();
    gain.pop_back();
    offset.pop_back();
    nonlin.pop_back();
    DIGITIZER.pop_back();
  }
  printf("Config file successfully read, constructing calibration file: %s\n", out);

  ofstream outfile(out, std::ofstream::out);

  for (int i = 0; i < channelnumber.size(); i++) {
    outfile << MNEMONIC.at(i) << " { \n";
    outfile << "Name:\t" << MNEMONIC.at(i) << "\n";
    outfile << "Number:\t" << channelnumber.at(i) << "\n";
    outfile << "Address:\t" << ADDRESS.at(i) << "\n";
    outfile << "Digitizer:\t" << DIGITIZER.at(i) << "\n";
    outfile << "EngCoeff:\t" << offset.at(i) << " " << gain.at(i) << " " << nonlin.at(i) << "\n";
    outfile << "Integration:\t" << 0 << "\n";
    outfile << "ENGChi2:\t" << 0 << "\n";
    outfile << "FileInt:\t" << 0 << "\n";
    outfile << "}\n";
    outfile << "\n";
    outfile << "//====================================//\n";
  }

  cout << "Construction complete!!!\n";

  outfile.close();

}

int main(int argc, char ** argv) {

  printf("Constructing CalFile\n");
  if (argc == 1) {
    printf("Assuming default Configuration File Name: Conf_File.txt\n Output calibration file (default: CalibrationFile.cal)\n");
    ConstructCalFile();
  } else if (argc == 2) {
    printf("Using Configuration File: %s,\n Output calibration file (default: CalibrationFile.cal)\n", argv[1]);
    ConstructCalFile(argv[1]);
  } else if (argc == 3) {
    printf("Using Configuration File: %s,\n Output calibration file %s\n", argv[1], argv[2]);
    ConstructCalFile(argv[1], argv[2]);
  } else if (argc == 4) {
    printf("Using Configuration File: %s,\nOutput calibration file %s\nEMMA Configuration File %s\n", argv[1], argv[2], argv[3]);
    ConstructCalFile(argv[1], argv[2], argv[3]);
  }
  else if (argc > 4) {
    printf("Too Many Arguments! Max 3, Configuration File Name, Output calibration File Name, Optional EMMA Configuration Name\n");
  }
  return 0;

}
