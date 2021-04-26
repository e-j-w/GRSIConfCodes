//g++ Conf_Setup.cxx -std=c++0x -o SetupConfFile

#include "CalibrationParameters.h"
#include "functions.h"

using namespace std;

int num_known_exp = 8; //number of known experiments
string exp_names[] = {
  "tigress",
  "emma",
  "emmas3",
  "griffin",
  "bambino",
  "descant",
  "tip",
  "emmatip",
};
string exp_description[] = {
  "TIGRESS 16 Clover Configuration",
  "TIGRESS 12 Clovers + SSB + EMMA Focal Plane",
  "TIGRESS 12 Clovers + Upstream S3 + SSB + EMMA Focal Plane",
  "GRIFFIN 16 Clover + PACES + SCEPTAR + ZDS + LaBr3",
  "TIGRESS 16 Clover + 2 S3 Detectors Standard Arrangment",
  "GRIFFIN 16 Clover + PACES + SCEPTAR + ZDS + LaBr3 + DESCANT",
  "TIGRESS 16 Clover + TIP",
  "TIGRESS 12 Clovers + TIP + EMMA Focal Plane",
};

int search_array(string array[], string search, int len) {
  for (int i = 0; i < len; i++) {
    if (search == array[i]) {
      return i;
    }//if
  }//for
  return -1;
}//search_array

bool CreateConfFile(const char * experiment, const char * MSC = "NULL", const char * seginp = "NULL") {

  int collector, port, channel;
  int chancounter = 0;
  char line[128];
  const char * outfile = "Conf_File.txt";
  int ring[24] = {23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0}; // Ring Order of S3s
  int sector[32] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31}; // Sector Order of S3s

  std::vector<std::string> MNEMONIC;
  std::string tempmnemonic;
  std::vector<int> customcollector, customport, customchannel;
  ifstream custom;
  custom.open("Source/custom.dat");
  if(custom.is_open()){
    while(!custom.eof()){
      custom >> std::hex >> tempmnemonic >> collector >> port >> channel;
      MNEMONIC.push_back(tempmnemonic);
      customcollector.push_back(collector);
      customport.push_back(port);
      customchannel.push_back(channel);
    }
  }

  int expID = search_array(exp_names, experiment, num_known_exp);

  switch (expID) {
  case 0:
    printf("Creating TIGRESS 16 Clovers\n");
    zerogains(gain, offset, non_lin, (sizeof(gain)/sizeof(gain[0])));
    loadSegmentPar(seginp, seggains, segoffsets);
    chancounter = makeTIGRESS(1, 16, chancounter, outfile, MSC, gain, offset, non_lin, seggains, segoffsets, MNEMONIC, customcollector, customport, customchannel);
    chancounter = makeRF(chancounter, outfile, MSC, 1, 0, 15);
    break;

  case 1:

    printf("Creating TIGRESS-EMMA 12 Clovers\n");
    zerogains(gain, offset, non_lin, (sizeof(gain)/sizeof(gain[0])));
    loadSegmentPar(seginp, seggains, segoffsets);
    chancounter = makeTIGRESS(5, 16, chancounter, outfile, MSC, gain, offset, non_lin, seggains, segoffsets, MNEMONIC, customcollector, customport, customchannel);
    chancounter = makeRF(chancounter, outfile, MSC, 1, 0, 15);
    chancounter = makeEMMAMisc(chancounter, outfile, MSC);
    break;


  case 2:

    printf("Creating TIGRESS-EMMA 12 Clovers + Upstream S3\n");
    zerogains(gain, offset, non_lin, (sizeof(gain)/sizeof(gain[0])));
    loadSegmentPar(seginp, seggains, segoffsets);
    chancounter = makeTIGRESS(5, 16, chancounter, outfile, MSC, gain, offset, non_lin, seggains, segoffsets, MNEMONIC, customcollector, customport, customchannel);

    zerogains(s3gains, s3offsets, (sizeof(s3gains)/sizeof(s3gains[0])));
    chancounter = makeS3Emma(chancounter, outfile, MSC, s3gains, s3offsets, ring, sector, MNEMONIC, customcollector, customport, customchannel);
    chancounter = makeRF(chancounter, outfile, MSC, 1, 0, 15);
    chancounter = makeEMMAMisc(chancounter, outfile, MSC);
    break;

  case 3:

    printf("Creating GRIFFIN 16 Clovers + PACES + ZDS + SCEPTAR + Labr3\n");
    zerogains(gain, offset, non_lin, (sizeof(gain)/sizeof(gain[0])));
    chancounter = makeGRIFFIN(chancounter, outfile, MSC, gain, offset, non_lin, MNEMONIC, customcollector, customport, customchannel);
    chancounter = makeSCEPTAR(chancounter, outfile, MSC, MNEMONIC, customcollector, customport, customchannel);
    chancounter = makeZDS(chancounter, outfile, MSC, 2, 2, 0, MNEMONIC, customcollector, customport, customchannel);

    zerogains(paces_gain, paces_offset, paces_nonlin, (sizeof(paces_gain)/sizeof(paces_gain[0])));
    chancounter = makePACES(chancounter, outfile, MSC, paces_gain, paces_offset, paces_nonlin, MNEMONIC, customcollector, customport, customchannel);

    zerogains(LaBr_gain, LaBr_offset, LaBr_nonlin, (sizeof(LaBr_gain)/sizeof(LaBr_gain[0])));
    chancounter = makeLaBr3(chancounter, outfile, MSC, LaBr_gain, LaBr_offset, LaBr_nonlin, MNEMONIC, customcollector, customport, customchannel);
    break;

  case 4:

    printf("Creating TIGRESS + Standard Bambino \n");
    loadSegmentPar(seginp, seggains, segoffsets);
    zerogains(gain, offset, non_lin, (sizeof(gain)/sizeof(gain[0])));
    chancounter = makeTIGRESS(1, 16, chancounter, outfile, MSC, gain, offset, non_lin, seggains, segoffsets, MNEMONIC, customcollector, customport, customchannel);

    zerogains(s3gains, s3offsets, (sizeof(s3gains)/sizeof(s3gains[0])));
    chancounter = makeS3Bambino(chancounter, outfile, MSC, s3gains, s3offsets, ring, sector, MNEMONIC, customcollector, customport, customchannel);
    chancounter = makeRF(chancounter, outfile, MSC, 1, 0, 15);
    break;

  case 5:

    printf("Creating GRIFFIN 16 Clovers + PACES + ZDS + SCEPTAR + Labr3 + DESCANT\n");
    zerogains(gain, offset, non_lin, (sizeof(gain)/sizeof(gain[0])));
    chancounter = makeGRIFFIN(chancounter, outfile, MSC, gain, offset, non_lin, MNEMONIC, customcollector, customport, customchannel);
    chancounter = makeSCEPTAR(chancounter, outfile, MSC, MNEMONIC, customcollector, customport, customchannel);
    chancounter = makeZDS(chancounter, outfile, MSC, 2, 2, 0, MNEMONIC, customcollector, customport, customchannel);

    zerogains(paces_gain, paces_offset, paces_nonlin, (sizeof(paces_gain)/sizeof(paces_gain[0])));
    chancounter = makePACES(chancounter, outfile, MSC, paces_gain, paces_offset, paces_nonlin, MNEMONIC, customcollector, customport, customchannel);

    zerogains(LaBr_gain, LaBr_offset, LaBr_nonlin, (sizeof(LaBr_gain)/sizeof(LaBr_gain[0])));
    chancounter = makeLaBr3(chancounter, outfile, MSC, LaBr_gain, LaBr_offset, LaBr_nonlin, MNEMONIC, customcollector, customport, customchannel);
    chancounter = makeDESCANT(chancounter, outfile, MSC, MNEMONIC, customcollector, customport, customchannel);
    break;

  case 6:

    printf("Creating TIGRESS + TIP \n");
    zerogains(gain, offset, non_lin, (sizeof(gain)/sizeof(gain[0])));
    loadSegmentPar(seginp, seggains, segoffsets);
    chancounter = makeTIGRESS(1, 16, chancounter, outfile, MSC, gain, offset, non_lin, seggains, segoffsets, MNEMONIC, customcollector, customport, customchannel);
    chancounter = makeTIP(chancounter, outfile, MSC, csigains, csioffsets, MNEMONIC, customcollector, customport, customchannel);
    chancounter = makeRF(chancounter, outfile, MSC, 1, 0, 15);

    break;
  
  case 7:

    printf("Creating TIGRESS-EMMA 12 Clovers + TIP \n");
    zerogains(gain, offset, non_lin, (sizeof(gain)/sizeof(gain[0])));
    loadSegmentPar(seginp, seggains, segoffsets);
    chancounter = makeTIGRESS(5, 16, chancounter, outfile, MSC, gain, offset, non_lin, seggains, segoffsets, MNEMONIC, customcollector, customport, customchannel);
    chancounter = makeTIPEMMA(chancounter, outfile, MSC, csigains, csioffsets, MNEMONIC, customcollector, customport, customchannel);
    chancounter = makeRF(chancounter, outfile, MSC, 1, 0, 15);
    chancounter = makeEMMAMisc(chancounter, outfile, MSC);

    break;

  default:
    printf("Experiment not recognised\nKnown Experiments are:\n");
    for(int i = 0; i < num_known_exp; i++) printf("%s:\t%s\n",exp_names[i].c_str(),exp_description[i].c_str());
    return false;
  }

  if (strcmp(MSC, "NULL") != 0) {
    ofstream mscnames;
    mscnames.open(MSC,ios::app);
    sprintf(line, "trunc \"/DAQ/MSC/MSC\" '%i'", chancounter);
    mscnames << line << "\n";
    sprintf(line, "trunc \"/DAQ/MSC/chan\" '%i'", chancounter);
    mscnames << line << "\n";
    sprintf(line, "trunc \"/DAQ/MSC/datatype\" '%i'", chancounter);
    mscnames << line << "\n";
    sprintf(line, "trunc \"/DAQ/MSC/gain\" '%i'", chancounter);
    mscnames << line << "\n";
    sprintf(line, "trunc \"/DAQ/MSC/offset\" '%i'", chancounter);
    mscnames << line << "\n";
    sprintf(line, "trunc \"/DAQ/MSC/quadratic\" '%i'", chancounter);
    mscnames << line << "\n";
    sprintf(line, "trunc \"/DAQ/MSC/digitizer\" '%i'", chancounter);
    mscnames << line << "\n";

    mscnames.close();
  }
  return true;

}

int main(int argc, char * * argv) {

  bool success = false;
  if (argc == 1) {
    printf("Too few inputs, give experiment name\nKnown Experiments are:\n");
    for(int i = 0; i < num_known_exp; i++) printf("%s:\t%s\n",exp_names[i].c_str(),exp_description[i].c_str());
    return 0;
  }
  if (argc > 4) {
    printf("Too many inputs, max three\n");
    return 0;
  }

  if (argc == 2)
    success = CreateConfFile(argv[1]);

  if (argc == 3)
    success = CreateConfFile(argv[1], argv[2]);

  if (argc == 4)
    success = CreateConfFile(argv[1], argv[2], argv[3]);

  if (success) {
    printf("Config files created successfully!\n");
    return 1;
  } else {
    printf("Config creation failed!\n");
    return 0;
  }

  return 0;

}
