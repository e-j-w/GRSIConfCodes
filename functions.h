#include <stdio.h>
#include <string.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
using namespace std;

void zerogains(float gain[], float offset[], int len);
void zerogains(float gain[], float offset[], float non_lin[], int len);

/* TIGRESS CABLING MAP */
const int tigCollector[16]    = {3,3,3,3, 0,0,1,1, 2,2,2,2, 0,1,1,0}; //map of TIGRESS position to collector
const int tigCollectorPos[16] = {0,1,2,3, 2,0,3,2, 0,1,3,2, 1,1,0,3}; //map of TIGRESS position to collector sub-position (first, 2nd, 3rd, 4th in collector)

/* TIP CABLING MAP */
const int tipPos[128] = {1,2,3,4,5,6,7,8,
9,10,11,12,39,40,41,42,
49,50,51,52,53,54,55,56,
57,58,59,60,61,62,63,64,
65,66,67,68,69,70,71,72,
73,74,75,76,91,92,93,94,
117,118,119,120,95,96,97,98,
99,100,101,102,103,104,105,106,
13,14,15,16,17,18,19,20,
21,22,43,44,45,46,47,48,
23,24,25,26,27,28,29,30,
31,32,33,34,35,36,37,38,
121,122,123,124,125,126,127,128,
109,110,111,112,113,114,115,116,
107,108,77,78,79,80,81,82,
83,84,85,86,87,88,89,90};

//for TIG-10 setup only
//CsI ball test June 2018
//const int tigCollectorTig10[16] = {1,1,2,2, 3,-1,4,4, 3,-1,6,6,7,7,8,8};
//const int tigCollectorPosTig10[16] = {0,1,0,1, 0,0,0,1, 1,1,0,1, 0,1,0,1};
//68Se run Nov 2017
const int tigCollectorTig10[16] =    {1,1,2,2, 3,-1,4,4, 5,-1,5,3, 7,7,8,8};
const int tigCollectorPosTig10[16] = {0,1,0,1, 0,0,0,1,  0,1,1,1,   0,1,0,1};

/* SHARC CABLING MAP */
/* Each entry represents a bank of 8 SHARC channels */
#define SHARC_NUM_BANKS 96
//position, 1-16 (-1 for empty channels)
//positions 1-4 and 13-16 are quads, others are boxes
const int sharcPos[SHARC_NUM_BANKS]     = {5,5,5,5, 5,5,5,5, 5,6,6,6, 6,6,6,6, 6,6,7,7, 7,7,7,7, 7,7,7,8, 8,8,8,8, 8,8,8,8, 10,10,10,10, 10,10,10,10, 10,11,11,11, 11,11,11,11, 11,11,-2,-2, 13,13,13,13, 13,14,-1,-1, 14,14,15,15, 15,15,15,16, 16,16,16,16, 12,12,12,9, 12,12,12,12, 12,12,9,9, 9,9,9,9, 9,9,14,14};
//0=D (closest to target), 1=E, 2=F (furthest from target)
const int sharcDist[SHARC_NUM_BANKS]    = {0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0};
//polarity, 0=P, 1=N
const bool sharcIsN[SHARC_NUM_BANKS]    = {1,1,1,1, 1,1,0,0, 0,0,1,1, 1,1,1,1, 0,0,1,1, 1,1,1,1, 0,0,0,0, 1,1,1,1, 1,1,0,0, 1,1,1,1, 1,1,0,0, 0,0,1,1, 1,1,1,1, 0,0,0,0, 0,0,1,1, 1,1,0,1, 1,1,0,0, 1,1,1,1, 0,0,1,1, 0,0,0,0, 1,1,1,1, 1,1,1,1, 1,1,1,1, 0,0,0,0};
//channel number at start of bank
const int sharcStartCh[SHARC_NUM_BANKS] = {0,8,16,40, 24,32,0,8, 16,16,0,8, 16,40,24,32, 0,8,0,8, 16,40,24,32, 0,8,16,16, 0,8,16,40, 24,32,0,8, 0,8,16,40, 24,32,0,8, 16,16,0,8, 16,40,24,32, 0,8,0,8, 0,8,0,8, 16,16,0,8, 0,8,0,8, 0,8,16,16, 0,8,0,8, 0,8,16,16, 0,8,16,40, 24,32,0,8, 16,40,24,32, 0,8,0,8};

//ATSD adapter remapping arrays. This is how they need to be mapped:
/*
FMC chan | Preamp chan
1 		 | 1
2		 | 3
3 		 | 5
4 		 | 7
5 		 | 9
6 		 | 11
7 		 | 13
8 		 | 15
9 		 | 2
10 		 | 4
11 		 | 6
12		 | 8
13 		 | 10
14 		 | 12
15 		 | 14
16 		 | 16
*/

//we only need to remap those blocks that start at 0 and 24 since those correspond to the blocks of channels coming out of the first 16 channels on the adapter, which are the ones messed up
//the other channels (16-23 and 40-47) are jumpered, and the jumper mapping appears correct
const int FMCRemapSeg0[8] = {0,2,4,6,8,10,12,14};
const int FMCRemapSeg8[8] = {1,3,5,7,9,11,13,15}; //offset from above by 1 because we index at 0 on the detectors
const int FMCRemapSeg24[8] = {24,26,28,30,32,34,36,38};
const int FMCRemapSeg32[8] = {25,27,29,31,33,35,37,39}; //same as above, but just add 24 to every channel to remap 24-39. 


/*SHARC COLLECTOR MAP*/
const int sharcCollector[8] = {0,0,0,0,1,1,1,1};
const int sharcCollectorPos[8] = {0,1,2,3,0,1,2,3};



void write_to_msc(const char * msc, int chancounter, char electronicaddress[], char var[], int DetType, const char *digitizer) {
  char line[128];
  ofstream mscnames;
  if(chancounter == 0) {
    mscnames.open(msc);
  }
  else mscnames.open(msc,ios::app);

  sprintf(line, "set \"/DAQ/PSC/PSC[%i]\" '%s'", chancounter, electronicaddress);
  mscnames << line << "\n";
  sprintf(line, "set \"/DAQ/PSC/chan[%i]\" '%s'", chancounter, var);
  mscnames << line << "\n";
  sprintf(line, "set \"/DAQ/PSC/datatype[%i]\" '%i'", chancounter, DetType);
  mscnames << line << "\n";
  sprintf(line, "set \"/DAQ/PSC/gain[%i]\" '1'", chancounter);
  mscnames << line << "\n";
  sprintf(line, "set \"/DAQ/PSC/offset[%i]\" '0'", chancounter);
  mscnames << line << "\n";
  sprintf(line, "set \"/DAQ/PSC/quadratic[%i]\" '0'", chancounter);
  mscnames << line << "\n";
  sprintf(line, "set \"/DAQ/PSC/digitizer[%i]\" '%s'", chancounter, digitizer);
  mscnames << line << "\n";
  mscnames.close();
}

int makeRF(int chancounter, const char *inp, const char *mscout, int collector, int port, int channel) {
  char line[128];
  char var[64];
  ofstream outfile;
  if(chancounter == 0) {
    outfile.open(inp);
  }
  else outfile.open(inp,ios::app);

  sprintf(var,"RFL00XS00x");
  char electronicaddress[32];
  sprintf(electronicaddress, "0x%01x%01x%02x", collector, port, channel);
  outfile << chancounter << "\t" << electronicaddress << "\t" <<  var << "\t" << 1 << "\t" << 0 << "\t" << 0 << "\tGRF16\n";
  if(strcmp(mscout, "NULL") != 0) {
    write_to_msc(mscout, chancounter, electronicaddress, var, 15, "GRF16");
  }
  outfile.close();
  chancounter++;
  return chancounter;
}

int makeTIGRESS(int first, int last, int portOffset, int chancounter, const char *inp, const char *mscout, float gain[64], float offset[64], float non_lin[64], float seggains[960], float segoffsets[960], std::vector<std::string> MNEMONIC, std::vector<int> customcollector, std::vector<int> customport, std::vector<int> customchannel){
  char line[128];
  char var[64];
  int DetType;
  ofstream outfile;
  if(chancounter == 0) {
    outfile.open(inp);
  }
  else outfile.open(inp,ios::app);

  int num_clover = last - first + 1;
  cout << "Making " << num_clover << " TIGRESS clovers, positions " << first << " through " << last << "." << endl; 

  for (int i = 0; i < num_clover*64; i++) {
    int DetNum = (i / 64) + first;
    if(DetNum > 16)
      continue;
    int cryNum = (i % 64) / 16;
    int channel = i % 16;
    int collector = tigCollector[DetNum-1];
    int port = 0;
    if(cryNum == 0){
      //blue core
      //this is on one of the GRIF-16s with an FMC32 installed,
      //which are assigned to the first 4 ports of each collector
      port = tigCollectorPos[DetNum-1];
    }else{
      //green, red, white cores
      //starting at port 4, cabled as G,R,W,G,R,W,G,R,W,G,R,W
      port = 4 + (cryNum - 1) + (tigCollectorPos[DetNum-1]*3) + portOffset;
    }
    
    char electronicaddress[32];
    sprintf(electronicaddress, "0x%01x%01x%02x", collector, port, channel);
    //cout << "det: " << DetNum << ", address: " << electronicaddress << endl;
    if (channel < 15) {
      char colour[2];
      if (cryNum == 0) sprintf(colour, "B");
      else if (cryNum == 1) sprintf(colour, "G");
      else if (cryNum == 2) sprintf(colour, "R");
      else if (cryNum == 3) sprintf(colour, "W");

      if (channel < 8) {
        DetType = 2;
        sprintf(var, "TIG%2.2i%sP%2.2ix", DetNum, colour, channel + 1);
      } else if (channel == 8) {
        DetType = 1;
        sprintf(var, "TIG%2.2i%sN00B", DetNum, colour);
      } else if (channel == 9) {
        DetType = 0;
        sprintf(var, "TIG%2.2i%sN00A", DetNum, colour);
      } else if (channel < 15 && channel > 9) {
        DetType = 3;
        sprintf(var, "TIS%2.2i%sN%2.2ix", DetNum, colour, channel - 9);
      }
      
      for(int m = 0; m < MNEMONIC.size(); m++) {
        if (strcmp(var,MNEMONIC.at(m).c_str()) == 0) {
          sprintf(electronicaddress, "0x%01x%01x%02x", customcollector.at(m), customport.at(m), customchannel.at(m));
        }
      }
      int aNum = ((DetNum) - 1) * 4 + cryNum;
      int segnum = 8 * aNum + channel;
      if (channel < 8 ) outfile << chancounter << "\t" << electronicaddress << "\t" << var << "\t" << seggains[segnum] << "\t" << segoffsets[segnum] << "\t" << 0 << "\t" << "\tGRF16\n";
      else if (channel == 9) outfile << chancounter << "\t" << electronicaddress << "\t" << var << "\t" << gain[aNum] << "\t" << offset[aNum] << "\t" << non_lin[aNum] << "\t" << "\tGRF16\n";
      else outfile << chancounter << "\t" << electronicaddress << "\t" << var << "\t" << 1 << "\t" << 0 << "\t" << 0 << "\t" << "\tGRF16\n";
      if (strcmp(mscout, "NULL") != 0) {
        write_to_msc(mscout, chancounter, electronicaddress, var, DetType, "GRF16");
      }
      chancounter++;
    }
  }
  outfile.close();
  return chancounter;
}

int makeTIGRESS(int first, int last, int chancounter, const char *inp, const char *mscout, float gain[64], float offset[64], float non_lin[64], float seggains[960], float segoffsets[960], std::vector<std::string> MNEMONIC, std::vector<int> customcollector, std::vector<int> customport, std::vector<int> customchannel){
  return makeTIGRESS(first,last,0,chancounter,inp,mscout,gain,offset,non_lin,seggains,segoffsets,MNEMONIC,customcollector,customport,customchannel);
}

int makeGRIFFINatTIGRESS(int first, int last, int portOffset, int chancounter, const char *inp, const char *mscout, float gain[64], float offset[64], float non_lin[64], std::vector<std::string> MNEMONIC, std::vector<int> customcollector, std::vector<int> customport, std::vector<int> customchannel){

  char line[128];
  char var[64];
  int DetType;
  int cryNum;
  ofstream outfile;
  if(chancounter == 0) {
    outfile.open(inp);
  }
  else outfile.open(inp,ios::app);
  int num_clover = last - first + 1;
  cout << "Making " << num_clover << " GRIFFIN clovers @ TIGRESS, positions " << first << " through " << last << "." << endl;

  for (int i = 0; i < num_clover*32; i++){

    int DetNum = (i / 32) + first;
    int collector = tigCollector[DetNum-1];
    int cardNum = (i % 32) / 16;
    int port = 0;
    if(cardNum == 0){
      //this is on one of the GRIF-16s with an FMC32 installed,
      //which are assigned to the first 4 ports of each collector
      port = tigCollectorPos[DetNum-1];
    }else{
      port = 4 + (cardNum - 1) + (tigCollectorPos[DetNum-1]*3) + portOffset;
    }
    //int port = (i % 256) / 16;
    int channel = i % 16;
    int ab = cardNum%2;
    if(channel == 0) cryNum = 0;
    if(channel == 1) cryNum = 1;
    if(channel == 2) cryNum = 2;
    if(channel == 3) cryNum = 3;
    if(channel < 10 && channel > 4){
      if(ab==0) cryNum = 0;
      else cryNum = 2;
    }
    else if(channel < 15 && channel > 9){
      if(ab==0) cryNum = 1;
      else cryNum = 3;
    }
    char electronicaddress[32];
    sprintf(electronicaddress, "0x%01x%01x%02x", collector, port, channel);
    if (channel < 15 && channel != 4 ) {
      char colour[2];
      if (cryNum == 0) sprintf(colour, "B");
      else if (cryNum == 1) sprintf(colour, "G");
      else if (cryNum == 2) sprintf(colour, "R");
      else if (cryNum == 3) sprintf(colour, "W");

      if (channel < 4 && ab == 0) {
        DetType = 0;
        sprintf(var, "TIG%2.2i%sN00A", DetNum, colour);
      } else if (channel < 4 && ab != 0) {
        DetType = 1;
        sprintf(var, "TIG%2.2i%sN00B", DetNum, colour);
      } else if (channel > 4 && channel < 10) {
        DetType = 3;
        sprintf(var, "TIS%2.2i%sN%2.2iX", DetNum, colour,channel - 4);
      } else if (channel < 15 && channel > 9) {
        DetType = 3;
        sprintf(var, "TIS%2.2i%sN%2.2iX", DetNum, colour, channel - 9);
      }

      for(int m = 0; m < MNEMONIC.size(); m++) {
        if (strcmp(var,MNEMONIC.at(m).c_str()) == 0) {
          sprintf(electronicaddress, "0x%01x%01x%02x", customcollector.at(m), customport.at(m), customchannel.at(m));
        }
      }
      int aNum = ((DetNum) - 1) * 4 + cryNum;
      if (channel < 4) outfile << chancounter << "\t" << electronicaddress << "\t" << var << "\t" << gain[aNum] << "\t" << offset[aNum] << "\t" << non_lin[aNum] << "\t" << "\tGRF16\n";
      else outfile << chancounter << "\t" << electronicaddress << "\t" << var << "\t" << 1 << "\t" << 0 << "\t" << 0 << "\t" << "\tGRF16\n";
      if (strcmp(mscout, "NULL") != 0) {
        write_to_msc(mscout, chancounter, electronicaddress, var, DetType, "GRF16");

      }
      chancounter++;
    }
  }
  outfile.close();
  return chancounter;
}

int makeGRIFFINatTIGRESS(int first, int last, int chancounter, const char *inp, const char *mscout, float gain[64], float offset[64], float non_lin[64], std::vector<std::string> MNEMONIC, std::vector<int> customcollector, std::vector<int> customport, std::vector<int> customchannel){
  return makeGRIFFINatTIGRESS(first,last,0,chancounter,inp,mscout,gain,offset,non_lin,MNEMONIC,customcollector,customport,customchannel);
}

int makeTIGRESSTIG10(int chancounter, const char *inp, const char *mscout, float gain[64], float offset[64], float non_lin[64], float seggains[960], float segoffsets[960], std::vector<std::string> MNEMONIC, std::vector<int> customcollector, std::vector<int> customport, std::vector<int> customchannel){
  char line[128];
  char var[64];
  int DetType, cryNum;
  ofstream outfile;
  if(chancounter == 0) {
    outfile.open(inp);
  }
  else outfile.open(inp,ios::app);

  int num_clover = 16;
  cout << "Making " << num_clover << " TIGRESS clovers (TIG-10)." << endl; 

  for (int i = 0; i < num_clover*60; i++) {
    int DetNum = (i / 60) + 1;
    if(DetNum > 16)
      continue;
    int channel = i % 10;
    int collector = tigCollectorTig10[DetNum-1];
    if(collector < 0){
      continue;
    }
    int detChNum = (i % 60);
    int port = detChNum/10 + (tigCollectorPosTig10[DetNum-1]*6) + 1;
    if(detChNum < 15){
      cryNum = 0;
    }else if(detChNum < 30){
      cryNum = 1;
    }else if(detChNum < 45){
      cryNum = 2;
    }else{
      cryNum = 3;
    }
    
    char electronicaddress[32];
    sprintf(electronicaddress, "0x%01x00%01x%02x", collector, port, channel);
    //cout << "det: " << DetNum << ", address: " << electronicaddress << endl;
    if (channel < 10) {
      char colour[2];
      if (cryNum == 0) sprintf(colour, "B");
      else if (cryNum == 1) sprintf(colour, "G");
      else if (cryNum == 2) sprintf(colour, "R");
      else if (cryNum == 3) sprintf(colour, "W");

      if((detChNum==0)||(detChNum==20)||(detChNum==30)||(detChNum==50)){
        DetType = 0;
        sprintf(var, "TIG%2.2i%sN00A", DetNum, colour);
      }else if((detChNum==9)||(detChNum==29)||(detChNum==39)||(detChNum==59)){
        DetType = 1;
        sprintf(var, "TIG%2.2i%sN00B", DetNum, colour);
      }else if(((detChNum>=10)&&(detChNum<20))||((detChNum>=40)&&(detChNum<50))){
        DetType = 3;
        int suppCh = (channel%5)+1;
        sprintf(var, "TIS%2.2i%sN%2.2ix", DetNum, colour, suppCh);
      }else{
        DetType = 2;
        sprintf(var, "TIG%2.2i%sP%2.2ix", DetNum, colour, channel);
      }
      for(int m = 0; m < MNEMONIC.size(); m++) {
        if (strcmp(var,MNEMONIC.at(m).c_str()) == 0) {
          sprintf(electronicaddress, "0x%01x%01x%02x", customcollector.at(m), customport.at(m), customchannel.at(m));
        }
      }
      int aNum = ((DetNum) - 1) * 4 + cryNum;
      int calChNum = 60*(((collector-1)*2)+tigCollectorPosTig10[DetNum-1]) + detChNum;
      if (DetType==2) outfile << calChNum << "\t" << electronicaddress << "\t" << var << "\t" << 1 << "\t" << 0 << "\t" << 0 << "\t" << "\tTIG10\n";
      else if (DetType==0) outfile << calChNum << "\t" << electronicaddress << "\t" << var << "\t" << gain[aNum] << "\t" << offset[aNum] << "\t" << non_lin[aNum] << "\t" << "\tTIG10\n";
      else outfile << calChNum << "\t" << electronicaddress << "\t" << var << "\t" << 1 << "\t" << 0 << "\t" << 0 << "\t" << "\tTIG10\n";
      if (strcmp(mscout, "NULL") != 0) {
        write_to_msc(mscout, calChNum, electronicaddress, var, DetType, "TIG10");
      }
      chancounter++;
    }
  }
  outfile.close();
  return chancounter;
}


int makeGRIFFIN(int chancounter, const char *inp, const char *mscout, float gain[64], float offset[64], float non_lin[64], std::vector<std::string> MNEMONIC, std::vector<int> customcollector, std::vector<int> customport, std::vector<int> customchannel){

  char line[128];
  char var[64];
  int DetType;
  int cryNum;
  ofstream outfile;
  if(chancounter == 0) {
    outfile.open(inp);
  }
  else outfile.open(inp,ios::app);

  for (int i = 0; i < 512; i++) {
    int DetNum = (i / 32) + 1;
    int collector = (i / 256);
    int port = (i % 256) / 16;
    int channel = i % 16;
    int ab = port%2;
    if(channel == 0) cryNum = 0;
    if(channel == 1) cryNum = 1;
    if(channel == 2) cryNum = 2;
    if(channel == 3) cryNum = 3;
    if(channel < 10 && channel > 4){
      if(ab==0) cryNum = 0;
      else cryNum = 2;
    }
    else if(channel < 15 && channel > 9){
      if(ab==0) cryNum = 1;
      else cryNum = 3;
    }
    char electronicaddress[32];
    sprintf(electronicaddress, "0x%01x%01x%02x", collector, port, channel);
    if (channel < 15 && channel != 4 ) {
      char colour[2];
      if (cryNum == 0) sprintf(colour, "B");
      else if (cryNum == 1) sprintf(colour, "G");
      else if (cryNum == 2) sprintf(colour, "R");
      else if (cryNum == 3) sprintf(colour, "W");

      if (channel < 4 && ab == 0) {
        DetType = 0;
        sprintf(var, "GRG%2.2i%sN00A", DetNum, colour);
      } else if (channel < 4 && ab != 0) {
        DetType = 1;
        sprintf(var, "GRG%2.2i%sN00B", DetNum, colour);
      } else if (channel > 4 && channel < 10) {
        DetType = 7;
        sprintf(var, "GRS%2.2i%sN%2.2iX", DetNum, colour,channel - 4);
      } else if (channel < 15 && channel > 9) {
        DetType = 7;
        sprintf(var, "GRS%2.2i%sN%2.2iX", DetNum, colour, channel - 9);
      }

      for(int m = 0; m < MNEMONIC.size(); m++) {
        if (strcmp(var,MNEMONIC.at(m).c_str()) == 0) {
          sprintf(electronicaddress, "0x%01x%01x%02x", customcollector.at(m), customport.at(m), customchannel.at(m));
        }
      }
      int aNum = ((DetNum) - 1) * 4 + cryNum;
      if (channel < 4) outfile << chancounter << "\t" << electronicaddress << "\t" << var << "\t" << gain[aNum] << "\t" << offset[aNum] << "\t" << non_lin[aNum] << "\t" << "\tGRF16\n";
      else outfile << chancounter << "\t" << electronicaddress << "\t" << var << "\t" << 1 << "\t" << 0 << "\t" << 0 << "\t" << "\tGRF16\n";
      if (strcmp(mscout, "NULL") != 0) {
        write_to_msc(mscout, chancounter, electronicaddress, var, DetType, "GRF16");
      }
      chancounter++;
    }
  }
  outfile.close();
  return chancounter;
}

int makeSCEPTAR(int chancounter, const char *inp, const char *mscout, std::vector<std::string> MNEMONIC, std::vector<int> customcollector, std::vector<int> customport, std::vector<int> customchannel) {
  char line[128];
  char var[64];
  char electronicaddress[32];
  ofstream outfile;
  if(chancounter == 0) {
    outfile.open(inp);
  }
  else outfile.open(inp,ios::app);

  for(int i = 0; i < 32; i++) {
    int collector = 2;
    int port = (i % 256) / 16;
    int channel = i % 16;
    char electronicaddress[32];
    sprintf(electronicaddress, "0x%01x%01x%02x", collector, port, channel);
    if(channel < 10) {
      int DetNum = channel + 1 + port * 10;
      sprintf(var, "SEP%2.2iXN00A", DetNum);
      for(int m = 0; m < MNEMONIC.size(); m++) {
        if (strcmp(var,MNEMONIC.at(m).c_str()) == 0) {
          sprintf(electronicaddress, "0x%01x%01x%02x", customcollector.at(m), customport.at(m), customchannel.at(m));
        }
      }
      outfile << chancounter << "\t" << electronicaddress << "\t" << var << "\t" << 1 << "\t" << 0 << "\t" << 0 << "\tGRF16\n";
      if (strcmp(mscout, "NULL") != 0) {
        write_to_msc(mscout, chancounter, electronicaddress, var, 2, "GRF16");
      }
      chancounter++;
    }
  }
  return chancounter;
}

int makeZDS(int chancounter, const char *inp, const char *mscout, int collector, int port, int channel, std::vector<std::string> MNEMONIC, std::vector<int> customcollector, std::vector<int> customport, std::vector<int> customchannel) {
  char line[128];
  char var[64];
  ofstream outfile;
  if(chancounter == 0) {
    outfile.open(inp);
  }
  else outfile.open(inp,ios::app);

  sprintf(var,"ZDS01XN00A");
  char electronicaddress[32];
  sprintf(electronicaddress, "0x%01x%01x%02x", collector, port, channel);
  for(int m = 0; m < MNEMONIC.size(); m++) {
    if (strcmp(var,MNEMONIC.at(m).c_str()) == 0) {
      sprintf(electronicaddress, "0x%01x%01x%02x", customcollector.at(m), customport.at(m), customchannel.at(m));
    }
  }
  outfile << chancounter << "\t" << electronicaddress << "\t" <<  var << "\t" << 1 << "\t" << 0 << "\t" << 0 << "\tGRF16\n";
  if(strcmp(mscout, "NULL") != 0) {
    write_to_msc(mscout, chancounter, electronicaddress, var, 9, "GRF16");
  }
  outfile.close();
  chancounter++;
  return chancounter;
}
int makePACES(int chancounter, const char *inp, const char *mscout, float gain[], float offset[], float non_lin[], std::vector<std::string> MNEMONIC, std::vector<int> customcollector, std::vector<int> customport, std::vector<int> customchannel) {
  char line[128];
  char var[64];
  int collector = 2;
  int port = 2;
  char electronicaddress[32];

  ofstream outfile;
  if(chancounter == 0) {
    outfile.open(inp);
  }
  else outfile.open(inp,ios::app);

  for(int i = 0; i < 5; i++) {
    int channel = i + 11;
    int DetNum = i + 1;
    sprintf(var, "PAC%2.2iXN00A", DetNum);
    sprintf(electronicaddress, "0x%01x%01x%02x", collector, port, channel);
    for(int m = 0; m < MNEMONIC.size(); m++) {
      if (strcmp(var,MNEMONIC.at(m).c_str()) == 0) {
        sprintf(electronicaddress, "0x%01x%01x%02x", customcollector.at(m), customport.at(m), customchannel.at(m));
      }
    }
    outfile << chancounter << "\t" << electronicaddress << "\t" <<  var << "\t" << gain[i] << "\t" << offset[i] << "\t" << non_lin[i] << "\tGRF16\n";
    if(strcmp(mscout, "NULL") != 0) {
      write_to_msc(mscout, chancounter, electronicaddress, var, 5, "GRF16");
    }
    chancounter++;
  }
  outfile.close();
  return chancounter;
}

int makeLaBr3(int chancounter, const char *inp, const char *mscout, float gain[], float offset[], float non_lin[], std::vector<std::string> MNEMONIC, std::vector<int> customcollector, std::vector<int> customport, std::vector<int> customchannel){
  char line[128];
  char var[64];
  char electronicaddress[32];
  ofstream outfile;
  if(chancounter == 0) {
    outfile.open(inp);
  }
  else outfile.open(inp,ios::app);

  int collector = 2;
  int DetType;
  int DetNum;
  int cryNum;
  for(int i = 0; i < 48; i++) {
    int port = (i % 256) / 16 + 4;
    int channel = i % 16;
    char electronicaddress[32];
    sprintf(electronicaddress, "0x%01x%01x%02x", collector, port, channel);
    if(i < 32) {
      if(channel < 4) {
        DetType = 3;
        DetNum = channel + 1 + (port - 4)*4;
        sprintf(var, "LBL%2.2iXN00A", DetNum);
        for(int m = 0; m < MNEMONIC.size(); m++) {
          if (strcmp(var,MNEMONIC.at(m).c_str()) == 0) {
            sprintf(electronicaddress, "0x%01x%01x%02x", customcollector.at(m), customport.at(m), customchannel.at(m));
          }
        }
        outfile << chancounter << "\t" << electronicaddress << "\t" << var << "\t" << gain[DetNum-1] << "\t" << offset[DetNum-1] << "\t" << non_lin[DetNum-1] << "\tGRF16\n";
        if (strcmp(mscout, "NULL") != 0) {
          write_to_msc(mscout, chancounter, electronicaddress, var, DetType, "GRF16");
        }
        chancounter++;
      }
      else {
        DetType = 8;
        char colour[2];
        DetNum = (channel - 4) / 3 + 1 + 4*(port-4);
        cryNum = (channel - 4) % 3;
        if(cryNum == 0) sprintf(colour, "A");
        else if(cryNum == 1) sprintf(colour, "B");
        else if(cryNum == 2) sprintf(colour, "C");
        sprintf(var, "LBS%2.2i%sN00X", DetNum, colour);
        for(int m = 0; m < MNEMONIC.size(); m++) {
          if (strcmp(var,MNEMONIC.at(m).c_str()) == 0) {
            sprintf(electronicaddress, "0x%01x%01x%02x", customcollector.at(m), customport.at(m), customchannel.at(m));
          }
        }
        outfile << chancounter << "\t" << electronicaddress << "\t" << var << "\t" << 1 << "\t" << 0 << "\t" << 0 << "\tGRF16\n";
        if (strcmp(mscout, "NULL") != 0) {
          write_to_msc(mscout, chancounter, electronicaddress, var, DetType, "GRF16");
        }
        chancounter++;
      }
    }
    else {
      if(channel % 2 == 1) {
        DetType = 4;
        DetNum = channel / 2 + 1;
        sprintf(var, "LBT%2.2iXT00X", DetNum);
        for(int m = 0; m < MNEMONIC.size(); m++) {
          if (strcmp(var,MNEMONIC.at(m).c_str()) == 0) {
            sprintf(electronicaddress, "0x%01x%01x%02x", customcollector.at(m), customport.at(m), customchannel.at(m));
          }
        }
        outfile << chancounter << "\t" << electronicaddress << "\t" << var << "\t" << 1 << "\t" << 0 << "\t" << 0 << "\tGRF16\n";
        if (strcmp(mscout, "NULL") != 0) {
          write_to_msc(mscout, chancounter, electronicaddress, var, DetType, "GRF16");
        }
        chancounter++;
      }
    }
  }
  return chancounter;
}
// EMMA trigger and SSBs
int makeEMMAMisc(int chancounter, const char *inp, const char *mscout) {
  char line[128];
  char var[64];
  char electronicaddress[32];
  ofstream outfile;
  if(chancounter == 0) {
    outfile.open(inp);
  }
  else outfile.open(inp,ios::app);

  int DetType = 8;
  int channel = 15;
  int port = 4;
  int collector = 1;
  sprintf(var,"EMT00XP00x");
  sprintf(electronicaddress, "0x%01x%01x%02x", collector, port, channel);
  outfile << chancounter << "\t" << electronicaddress << "\t" <<  var << "\t" << 1 << "\t" << 0 << "\t" << 0 << "\tGRF16\n";
  if(strcmp(mscout, "NULL") != 0) {
    write_to_msc(mscout, chancounter, electronicaddress, var, DetType, "GRF16");
  }
  chancounter++;

  DetType=9;
  channel=15;
  for (int j = 0; j<2; j++) {
    port=5+j;
    sprintf(var,"ETO%02xXP00x",j);
    sprintf(electronicaddress, "0x%01x%01x%02x", collector, port, channel);
    outfile << chancounter << "\t" << electronicaddress << "\t" <<  var << "\t" << 1 << "\t" << 0 << "\t" << 0 << "\tGRF16\n";
    if(strcmp(mscout, "NULL") != 0) {
      write_to_msc(mscout, chancounter, electronicaddress, var, DetType, "GRF16");
    }
    chancounter++;
  }
  return chancounter;
}

// Upstream S3 in EMMA Chamber
int makeS3Emma(int chancounter, const int collector, const char *inp, const char *mscout, float s3gains[], float s3offsets[], int ring[], int sector[], std::vector<std::string> MNEMONIC, std::vector<int> customcollector, std::vector<int> customport, std::vector<int> customchannel){
  char line[128];
  char var[64];
  char electronicaddress[32];
  ofstream outfile;
  if(chancounter == 0) {
    outfile.open(inp);
  }
  else outfile.open(inp,ios::app);

  for (int i = 0; i < 60; i++) {
    int port = (i % 256)/16 ;
    int channel = i % 16;
    sprintf(electronicaddress, "0x%01x%01x%02x", collector, port, channel);
    if(i<32) {
      int DetType = 4;
      sprintf(var,"ETE01EP%02iX",sector[i]);
      for(int m = 0; m < MNEMONIC.size(); m++) {
        if (strcmp(var,MNEMONIC.at(m).c_str()) == 0) {
          sprintf(electronicaddress, "0x%01x%01x%02x", customcollector.at(m), customport.at(m), customchannel.at(m));
        }
      }
      outfile << chancounter << "\t" << electronicaddress << "\t" <<  var << "\t" << s3gains[i] << "\t" << s3offsets[i] << "\t" << 0 << "\tGRF16\n";
      if(strcmp(mscout, "NULL") != 0) {
        write_to_msc(mscout, chancounter, electronicaddress, var, DetType, "GRF16");
      }
      chancounter++;
    }
    else if (i>31 && i < 56) {
      int DetType = 5;
      sprintf(var,"ETE01EN%02iX",ring[i-32]);
      for(int m = 0; m < MNEMONIC.size(); m++) {
        if (strcmp(var,MNEMONIC.at(m).c_str()) == 0) {
          sprintf(electronicaddress, "0x%01x%01x%02x", customcollector.at(m), customport.at(m), customchannel.at(m));
       }
      }
      outfile << chancounter << "\t" << electronicaddress << "\t" <<  var << "\t" << s3gains[i] << "\t" << s3offsets[i] << "\t" << 0 << "\tGRF16\n";
      if(strcmp(mscout, "NULL") != 0) {
        write_to_msc(mscout, chancounter, electronicaddress, var, DetType, "GRF16");
      }
      chancounter++;
    }
  }
  return chancounter;
}

int makeS3Bambino(int chancounter, const char *inp, const char *mscout, float s3gains[], float s3offsets[], int ring[], int sector[], std::vector<std::string> MNEMONIC, std::vector<int> customcollector, std::vector<int> customport, std::vector<int> customchannel){
  char line[128];
  char var[64];
  char electronicaddress[32];
  ofstream outfile;
  if(chancounter == 0) {
    outfile.open(inp);
  }
  else outfile.open(inp,ios::app);

  for (int i = 0; i < 128; i++) {
    int port = (i % 256)/16 ;
    int channel = i % 16;
    int collector = (i/256) + 4;
    int DetNum =  i / 64 + 1;
    sprintf(electronicaddress, "0x%01x%01x%02x", collector, port, channel);
    if((i % 64)<32) {
      int DetType = 4;
      sprintf(var,"BAE%02iEP%02iX",DetNum, i % 64);
      for(int m = 0; m < MNEMONIC.size(); m++) {
        if (strcmp(var,MNEMONIC.at(m).c_str()) == 0) {
          sprintf(electronicaddress, "0x%01x%01x%02x", customcollector.at(m), customport.at(m), customchannel.at(m));
        }
      }
      outfile << chancounter << "\t" << electronicaddress << "\t" <<  var << "\t" << s3gains[i] << "\t" << s3offsets[i] << "\t" << 0 << "\tGRF16\n";
      if(strcmp(mscout, "NULL") != 0) {
        write_to_msc(mscout, chancounter, electronicaddress, var, DetType, "GRF16");
      }
      chancounter++;
    }
    else if ((i % 64) > 31 && (i % 64) < 56) {
      int DetType = 5;
      sprintf(var,"BAE%02iEN%02iX",DetNum,ring[(i % 64)-32]);
      for(int m = 0; m < MNEMONIC.size(); m++) {
        if (strcmp(var,MNEMONIC.at(m).c_str()) == 0) {
          sprintf(electronicaddress, "0x%01x%01x%02x", customcollector.at(m), customport.at(m), customchannel.at(m));
        }
      }
      outfile << chancounter << "\t" << electronicaddress << "\t" <<  var << "\t" << s3gains[i] << "\t" << s3offsets[i] << "\t" << 0 << "\tGRF16\n";
      if(strcmp(mscout, "NULL") != 0) {
        write_to_msc(mscout, chancounter, electronicaddress, var, DetType, "GRF16");
      }
      chancounter++;
    }
  }
  return chancounter;
}

int makeTIP(int chancounter, const char *inp, const char *mscout, float csigains[], float csioffsets[], std::vector<std::string> MNEMONIC, std::vector<int> customcollector, std::vector<int> customport, std::vector<int> customchannel){
  char line[128];
  char var[64];
  char electronicaddress[32];
  ofstream outfile;
  if(chancounter == 0) {
    outfile.open(inp);
  }
  else outfile.open(inp,ios::app);
  cout << "Making TIP CsI ball." << endl;

  for (int i = 0; i < 128; i++) {
    int port = (i % 64)/16;
    int channel = i % 16;
    int collector = (i/64) + 4;
	  //int collector = (i/256);
    int DetNum =  tipPos[i];
    sprintf(electronicaddress, "0x%01x%01x%02x", collector, port, channel);
    int DetType = 8;
    sprintf(var,"TPC%03iN00X",DetNum);
    for(int m = 0; m < MNEMONIC.size(); m++) {
      if (strcmp(var,MNEMONIC.at(m).c_str()) == 0) {
        sprintf(electronicaddress, "0x%01x%01x%02x", customcollector.at(m), customport.at(m), customchannel.at(m));
      }
    }
    outfile << chancounter << "\t" << electronicaddress << "\t" <<  var << "\t" << csigains[i] << "\t" << csioffsets[i] << "\t" << 0 << "\tGRF16\n";
    if(strcmp(mscout, "NULL") != 0) {
      write_to_msc(mscout, chancounter, electronicaddress, var, 12, "GRF16");
    }
    chancounter++;
  }
  return chancounter;
}

int makeTIPEMMA(int chancounter, const char *inp, const char *mscout, float csigains[], float csioffsets[], std::vector<std::string> MNEMONIC, std::vector<int> customcollector, std::vector<int> customport, std::vector<int> customchannel){
  char line[128];
  char var[64];
  char electronicaddress[32];
  ofstream outfile;
  if(chancounter == 0) {
    outfile.open(inp);
  }
  else outfile.open(inp,ios::app);

  for (int i = 0; i < 128; i++) {
    int port = (i % 256)/16;
    int channel = i % 16;
    int collector = (i/256);
    int DetNum =  i + 1;
    sprintf(electronicaddress, "0x%01x%01x%02x", collector, port, channel);
    int DetType = 8;
    sprintf(var,"TPC%03iN00X",DetNum);
    for(int m = 0; m < MNEMONIC.size(); m++) {
      if (strcmp(var,MNEMONIC.at(m).c_str()) == 0) {
        sprintf(electronicaddress, "0x%01x%01x%02x", customcollector.at(m), customport.at(m), customchannel.at(m));
      }
    }
    outfile << chancounter << "\t" << electronicaddress << "\t" <<  var << "\t" << csigains[i] << "\t" << csioffsets[i] << "\t" << 0 << "\tGRF16\n";
    if(strcmp(mscout, "NULL") != 0) {
      write_to_msc(mscout, chancounter, electronicaddress, var, 12, "GRF16");
    }
    chancounter++;
  }
  return chancounter;
}

// SHARC cabled into an FMC32 on the specified port and collector
// assumed to contain 4 banks of 8 channels, starting at startBank (should be a 0-indexed value
// based on the SHARC cabling map at the top of this file)
int makeSHARC(int chancounter, const int collector, const int port, const int startBank, const char *inp, const char *mscout, const float sharcgains[], const float sharcoffsets[], std::vector<std::string> MNEMONIC, std::vector<int> customcollector, std::vector<int> customport, std::vector<int> customchannel){
  cout << "Making SHARC."<< endl; 
  char line[128];
  char var[64];
  char electronicaddress[32];
  int DetType;
  ofstream outfile;
  if(chancounter == 0) {
    outfile.open(inp);
  }
  else outfile.open(inp,ios::app);
  
  int collectorstartch = 0;

  for (int i = 0; i < 32*8; i++) { //1 FMC32

    int bank = startBank + (i/8); //the 'bank of 8' that this channel is mapped to
    if(bank >= SHARC_NUM_BANKS){
      continue;
    }
    if(port > 15){
      printf("WARNING: Invalid SHARC port: %i\n",port);
      continue;
    }
    int bankCh = i % 8; //the channel number within the bank of 8
    int channel = 16 + (collectorstartch % (32));
    
    //this will adjust the segment mapping to follow funky mapping of the first 16 channels on the ATSD to Samtec adapters
    int adjustedSegment = sharcStartCh[bank]+bankCh; //the existing mapping
    //this if statement should only trigger on the banks that need remapping
    if (0 == sharcStartCh[bank]){ //mappings need to be translated from FMC to preamp: FMC 1 -> segment 0, FMC 2 -> seg 3, FMC 3 -> seg 5,..., FMC 9 -> seg2, FMC 10 -> seg4,..., FMC 16->seg 16
    	adjustedSegment = FMCRemapSeg0[bankCh]; //this should remap the segments 0-15 correctly into the FMC now
    }
    else if (8 == sharcStartCh[bank]){
    	adjustedSegment = FMCRemapSeg8[bankCh]; //this should remap the segments 24-39 correctly into the FMC now
    }
    else if (24 == sharcStartCh[bank]){
    	adjustedSegment = FMCRemapSeg24[bankCh]; //this should remap the segments 24-39 correctly into the FMC now
    }
    else if (32 == sharcStartCh[bank]){
    	adjustedSegment = FMCRemapSeg32[bankCh]; //this should remap the segments 24-39 correctly into the FMC now
    }
    
    //from now on, adjustedSegment should be correct for every channel regardless of it is mapped weirdly by the ATSD preamp or not.
    //I will now replace sharcStartCh[bank]+bankCh with adjustedSegment

    int collectorNum = sharcCollector[i/(32)];
    int collectorPort = sharcCollectorPos[i/(32)];
    
    sprintf(electronicaddress, "0x%01x%01x%02x", collectorNum, collectorPort, channel);
    if(sharcIsN[bank]){
      DetType = 6;
    }else{
      DetType = 7;
    }
    if(sharcPos[bank]>4 && sharcPos[bank]<13){
      if(sharcIsN[bank]){
        switch(sharcDist[bank]){
          case 0:
          default:
            sprintf(var,"SHB%02iDN%02iX",sharcPos[bank],adjustedSegment);
            break;
          case 1:
            sprintf(var,"SHB%02iEN%02iX",sharcPos[bank],adjustedSegment);
            break;
          case 2:
            sprintf(var,"SHB%02iFN%02iX",sharcPos[bank],adjustedSegment);
            break;
        }
      }else{
        switch(sharcDist[bank]){
          case 0:
          default:
            sprintf(var,"SHB%02iDP%02iX",sharcPos[bank],adjustedSegment);
            break;
          case 1:
            sprintf(var,"SHB%02iEP%02iX",sharcPos[bank],adjustedSegment);
            break;
          case 2:
            sprintf(var,"SHB%02iFP%02iX",sharcPos[bank],adjustedSegment);
            break;
        }
      }
    }else if(sharcPos[bank]>0){
      if(sharcIsN[bank]){
        switch(sharcDist[bank]){
          case 0:
          default:
            sprintf(var,"SHQ%02iDN%02iX",sharcPos[bank],adjustedSegment);
            break;
          case 1:
            sprintf(var,"SHQ%02iEN%02iX",sharcPos[bank],adjustedSegment);
            break;
          case 2:
            sprintf(var,"SHQ%02iFN%02iX",sharcPos[bank],adjustedSegment);
            break;
        }
      }else{
        switch(sharcDist[bank]){
          case 0:
          default:
            sprintf(var,"SHQ%02iDP%02iX",sharcPos[bank],adjustedSegment);
            break;
          case 1:
            sprintf(var,"SHQ%02iEP%02iX",sharcPos[bank],adjustedSegment);
            break;
          case 2:
            sprintf(var,"SHQ%02iFP%02iX",sharcPos[bank],adjustedSegment);
            break;
        }
      }
    }else if(sharcPos[bank]==-2){
    	sprintf(var,"SHB%02iEP00X",adjustedSegment);
      /*if(bankCh == 2){
        i += 3;
        collectorstartch += 4;
        continue;
      }else{
        switch(bankCh){
          	continue;
          case 7:
          	std::cout << "\ncase 7: bank chan is " << bankCh;
            sprintf(var,"SHB05EN00X");
            break;
          case 6:
          	std::cout << "\ncase 6: bank chan is " << bankCh;
            sprintf(var,"SHB06EN00X");
            break;
          case 1:
          	std::cout << "\ncase 1: bank chan is " << bankCh;
            sprintf(var,"SHB07EN00X");
            break;
          case 0:
          	std::cout << "\ncase 0: bank chan is " << bankCh;
            sprintf(var,"SHB08EN00X");
            break;
        }
      }*/
    }else{
      collectorstartch++;
      continue;
    }
    
    for(int m = 0; m < MNEMONIC.size(); m++) {
      if(strcmp(var,MNEMONIC.at(m).c_str()) == 0) {
        sprintf(electronicaddress, "0x%01x%01x%02x", customcollector.at(m), customport.at(m), customchannel.at(m));
      }
    }
    outfile << chancounter << "\t" << electronicaddress << "\t" <<  var << "\t" << sharcgains[i] << "\t" << sharcoffsets[i] << "\t" << 0 << "\tGRF16\n";
    if(strcmp(mscout, "NULL") != 0) {
      write_to_msc(mscout, chancounter, electronicaddress, var, DetType, "GRF16");
    }
    chancounter++;
    collectorstartch++;
    
  }
  return chancounter;
}

int makeTRIFIC(int chancounter, const int startcollector, const int startport, const char *inp, const char *mscout, float trificgains[], float trificoffsets[], std::vector<std::string> MNEMONIC, std::vector<int> customcollector, std::vector<int> customport, std::vector<int> customchannel){
  cout << "Making TRIFIC."<< endl; 
  char line[128];
  char var[64];
  char electronicaddress[32];
  ofstream outfile;
  if(chancounter == 0) {
    outfile.open(inp);
  }
  else outfile.open(inp,ios::app);

  int collector = startcollector;
  int collectorstartport = startport;
  int collectorstartch = 0;

  for (int i = 0; i < 64; i++) {

    int port = collectorstartport + collectorstartch/16;
    int channel = collectorstartch % 16;

    //roll over to the next collector if neccesary
    if(port>15){
      collector++;
      port=0;
      channel=0;
      collectorstartport=0;
      collectorstartch=0;
    }

    int DetType = 10;
    switch(i/16){
      case 3:
        sprintf(var,"TFC03YN%02iX",channel);
        DetType = 11;
        break;
      case 2:
        sprintf(var,"TFC05XN%02iX",channel);
        DetType = 11;
        break;
      case 1:
        sprintf(var,"TFC%02iSN00X",(channel*2)+1);
        DetType = 11;
        break;
      case 0:
        sprintf(var,"TFC%02iSP00X",(channel+1)*2);
        DetType = 10;
        break;
      default:
        cout << "Unknown TRIFIC channel: " << i << endl;
        chancounter++;
        collectorstartch++;
        continue;
    }

    sprintf(electronicaddress, "0x%01x%01x%02x", collector, port, channel);
    for(int m = 0; m < MNEMONIC.size(); m++) {
      if (strcmp(var,MNEMONIC.at(m).c_str()) == 0) {
        sprintf(electronicaddress, "0x%01x%01x%02x", customcollector.at(m), customport.at(m), customchannel.at(m));
      }
    }
    outfile << chancounter << "\t" << electronicaddress << "\t" <<  var << "\t" << trificgains[i] << "\t" << trificoffsets[i] << "\t" << 0 << "\tGRF16\n";
    if(strcmp(mscout, "NULL") != 0) {
      write_to_msc(mscout, chancounter, electronicaddress, var, DetType, "GRF16");
    }
    chancounter++;
    collectorstartch++;
  }
  return chancounter;
}


int makeDESCANT(int chancounter, const char *inp, const char *mscout, std::vector<std::string> MNEMONIC, std::vector<int> customcollector, std::vector<int> customport, std::vector<int> customchannel) {
  char line[128];
  char var[64];
  char electronicaddress[32];
  ofstream outfile;
  if(chancounter == 0) {
    outfile.open(inp);
  }
  else outfile.open(inp,ios::app);

  for (int i = 0; i < 70; i++) {
    int collector = 8 + i / 48;
    int DetNum = i + 1;
    int port = ((i / 16 + 1) * 4) % 16;
    int channel = i % 16;

    sprintf(electronicaddress, "0x%01x%01x%02x", collector, port, channel);
    int DetType = 10;
    sprintf(var, "DSC%2.2iXN00X", DetNum);
    for(int m = 0; m < MNEMONIC.size(); m++) {
      if (strcmp(var,MNEMONIC.at(m).c_str()) == 0) {
        sprintf(electronicaddress, "0x%01x%01x%02x", customcollector.at(m), customport.at(m), customchannel.at(m));
      }
    }
    outfile << chancounter << "\t" << electronicaddress << "\t" << var << "\t" << 1 << "\t" << 0 << "\t" << 0 << "\tCAEN\n";
    if (strcmp(mscout, "NULL") != 0) {
      write_to_msc(mscout, chancounter, electronicaddress, var, DetType, "CAEN");
    }
    chancounter++;
  }
  return chancounter;
}

// Loads TIGRESS SegmentCalibrationParmeter
void loadSegmentPar(const char *inp, float gain[], float offset[]){
  if (strcmp(inp, "NULL") != 0) {
    ifstream segpar;
    segpar.open(inp);
    if (segpar.is_open()) {
      printf("Segment parameter file: %s opened!\n", inp);
      int j = 0;
      while (!segpar.eof() && j < 512) {
        segpar >> offset[j] >> gain[j];
        j++;
      }
    } else {
      printf("%s not opened\n", inp);
      return;
    }
    printf("Segment parameters read in successfully\n");
  } else {
    printf("No segment parameter file declared, setting gains, offsets and nonlinear components to zero\n");
    for (int i = 0; i < 512; i++) {
      gain[i] = 1;
      offset[i] = 0;
    }
  }
}

void zerogains(float gain[], float offset[], int len) {
  for(int i = 0; i < len; i++) {
    if (gain[i] == 0 && offset[i] == 0) {
      gain[i] = 1;
      offset[i] = 0;
    }
  }
}

void zerogains(float gain[], float offset[], float non_lin[], int len) {
  for(int i = 0; i < len; i++) {
    if (gain[i] == 0 && offset[i] == 0) {
      gain[i] = 1;
      offset[i] = 0;
      non_lin = 0;
    }
  }
}

