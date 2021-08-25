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
const int tigCollector[16]    = {0,1,2,2,3,0,3,1,3,2,3,2,0,1,1,0}; //map of TIGRESS position to collector
const int tigCollectorPos[16] = {1,3,0,2,0,0,1,2,2,1,3,3,2,1,0,3}; //map of TIGRESS position to collector sub-position (first, 2nd, 3rd, 4th in collector)

void write_to_msc(const char * msc, int chancounter, char electronicaddress[], char var[], int DetType, const char *digitizer) {
  char line[128];
  ofstream mscnames;
  if(chancounter == 0) {
    mscnames.open(msc);
  }
  else mscnames.open(msc,ios::app);

  sprintf(line, "set \"/DAQ/MSC/MSC[%i]\" '%s'", chancounter, electronicaddress);
  mscnames << line << "\n";
  sprintf(line, "set \"/DAQ/MSC/chan[%i]\" '%s'", chancounter, var);
  mscnames << line << "\n";
  sprintf(line, "set \"/DAQ/MSC/datatype[%i]\" '%i'", chancounter, DetType);
  mscnames << line << "\n";
  sprintf(line, "set \"/DAQ/MSC/gain[%i]\" '1'", chancounter);
  mscnames << line << "\n";
  sprintf(line, "set \"/DAQ/MSC/offset[%i]\" '0'", chancounter);
  mscnames << line << "\n";
  sprintf(line, "set \"/DAQ/MSC/quadratic[%i]\" '0'", chancounter);
  mscnames << line << "\n";
  sprintf(line, "set \"/DAQ/MSC/digitizer[%i]\" '%s'", chancounter, digitizer);
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
  //cout << "first: " << first << ", last: " << last << ", num_clover: " << num_clover << endl;

  for (int i = 0; i < num_clover*64; i++) {
    int DetNum = (i / 64) + first;
    if(DetNum > num_clover)
      continue;
    //int port = (i % 256) / 16;
    int channel = i % 16;
    //int collector = (i / 256) + first / 5;
    int collector = tigCollector[DetNum-1];
    int port = (i % 64) / 16 + (tigCollectorPos[DetNum-1]*4) + portOffset;
    int cryNum = (port % 4);
    char electronicaddress[32];
    sprintf(electronicaddress, "0x%01x%01x%02x", collector, port, channel);
    //cout << "det: " << DetNum << ", address: " << electronicaddress << endl;
    if (channel < 15) {
      char colour[1];
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
  makeTIGRESS(first,last,0,chancounter,inp,mscout,gain,offset,non_lin,seggains,segoffsets,MNEMONIC,customcollector,customport,customchannel);
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

  for (int i = 0; i < num_clover*32; i++) {
    int DetNum = (i / 32) + 1;
    int collector = tigCollector[DetNum-1];
    //int collector = (i / 256);
    int port = (i % 32) / 16 + (tigCollectorPos[DetNum-1]*4) + portOffset;
    //int port = (i % 256) / 16;
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
      char colour[1];
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

int makeGRIFFINatTIGRESS(int first, int last, int chancounter, const char *inp, const char *mscout, float gain[64], float offset[64], float non_lin[64], std::vector<std::string> MNEMONIC, std::vector<int> customcollector, std::vector<int> customport, std::vector<int> customchannel){
  makeGRIFFINatTIGRESS(first,last,0,chancounter,inp,mscout,gain,offset,non_lin,MNEMONIC,customcollector,customport,customchannel);
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
      char colour[1];
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
        char colour[1];
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
  int port = 1;
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
    port=2+j;
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
int makeS3Emma(int chancounter, const char *inp, const char *mscout, float s3gains[], float s3offsets[], int ring[], int sector[], std::vector<std::string> MNEMONIC, std::vector<int> customcollector, std::vector<int> customport, std::vector<int> customchannel){
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
    int collector = (i/256);
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

  for (int i = 0; i < 128; i++) {
    int port = (i % 256)/16;
    int channel = i % 16;
    //int collector = (i/256) + 4;
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

