// g++ chan_check.cxx -std=c++0x -I$GRSISYS/include -L$GRSISYS/libraries `grsi-config --cflags --all-libs` `root-config --cflags --libs` -lTreePlayer -lSpectrum -o ChanCheck

#include <iostream>
#include <iomanip>
#include "TCutG.h"
#include "TH1.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TTree.h"
#include "TChain.h"
#include "TH2.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TList.h"
#include "TSpectrum.h"
#include "TChannel.h"
#include "TFragment.h"
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm> 
#include <vector> 
#include "grifadc.h" 
#include "tigadc.h"

using namespace std;

void chan_check(char const* infile, char const* conffile, char const* calfile, char const* outfile){
  // Ignores certain Channels if changed to 0
  // HPGe A, HPGe B, Suppressors, Segments, Ancillary
  int ignore[5] = {1,1,1,1,1};

  // Definitions
  std::string tempadd, tempmnemonic,tempdig, tempoutsensor, grifcN, grif16N, grif16Ch, tempsubsys, tempcharge;
  double tempgain,tempoffset,tempnonlin,tempdtype;
  int tempchan, tempgrifc, tempgrif16, tempgrifch; 
  std::vector<std::string> ADDRESS,MNEMONIC,DIGITIZER;
  std::vector<int> channelnumber;
  std::vector<double> gain,offset,nonlin;
  char var[64];
  int j=0;
  Double_t nBins = 4096;
  Double_t bmin = 0;
  Double_t bmax = 4096;

  // Loads Configuration File
  printf("Loading Configuration File: %s\n",conffile);
  ifstream config;
  config.open(conffile);
  if(config.is_open()){
    printf("Configuration file: %s opened\n",conffile);
    while(!config.eof()){
      config >> tempchan >> tempadd >> tempmnemonic >> tempgain >> tempoffset >> tempnonlin >> tempdig;
      channelnumber.push_back(tempchan);
      ADDRESS.push_back(tempadd);
      MNEMONIC.push_back(tempmnemonic);
      gain.push_back(tempgain);
      offset.push_back(tempoffset);
      nonlin.push_back(tempnonlin);
      DIGITIZER.push_back(tempdig);
      j++;
    }	
  } else{
    printf("Conf file %s failed to open, abort!\n",conffile); 
    return;	
  }

  // Checks whether it is a TIGRESS or GRIFFIN Configuration file 
  string experiment;
  for(int i = 0; i < MNEMONIC.size();i++){
    if( MNEMONIC.at(i).compare(0,2,"TI") == 0 ) {
      experiment = "tigress";
      break;
    }
    if( MNEMONIC.at(i).compare(0,2,"GR") == 0 ) {
      experiment = "griffin";
      break;
    }
  }

  //Define Histograms
  TList * list = new TList;
  TH1D * h[5000];
  char shname[20];
  for (int iii = 0; iii < 2000; iii++) {
    sprintf(shname, "Channel%d", iii);
    h[iii] = new TH1D(shname, Form("Channel d%1.1i", iii), nBins, bmin, bmax);
  }
  TH2F *hh = new TH2F("hh","",2000,0,2000,nBins, bmin, bmax);list->Add(hh);
  TFile * inputfile = new TFile(infile, "READ");
  if (!inputfile->IsOpen()) {
    printf("Opening file failed, aborting\n");
    return;
  }

  //Loads FragmentTree and Calibration File
  TChain * FragmentTree = (TChain * ) inputfile->Get("FragmentTree");
  printf("%i tree files, details:\n", FragmentTree->GetNtrees());
  FragmentTree->ls();
  TTree * tree = (TTree * ) FragmentTree->GetTree();
  Int_t nentries = FragmentTree->GetEntries();
  TFragment * frag = 0;
  FragmentTree->SetBranchAddress("TFragment", & frag);
  TChannel::ReadCalFile(calfile);

  const char *out = "missing.txt";
  ofstream ofile;
  ofile.open(out);
  if (!ofile.is_open()) {
    printf("Output not opened\n");
    return;
  }

  //Fills Histograms
  printf("Begin sort\n");
  int one;
  for (int jentry = 0; jentry < (nentries)/10; jentry++) {
    tree->GetEntry(jentry);
    h[frag->GetChannelNumber()]->Fill(frag->GetCharge());
    hh->Fill(frag->GetChannelNumber(),frag->GetCharge());
    if (jentry % 10000 == 0) cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << nentries << ", " << 100 * jentry / nentries << "% complete" << "\r" << flush;
  }
  cout << "Entry " << nentries << " of " << nentries << ", 100% Complete!\n";
  cout << "Histograms written, sorting complete" << endl;

  // Checks which channels are missing and outputs to a text file the missing channel name and the adc number it is plugged into
  for (int iii = 0; iii < j; iii++) {
    if(h[iii]->Integral(10,nBins) < 100) {
      tempmnemonic = MNEMONIC.at(iii);
      tempadd = ADDRESS.at(iii);
      tempoutsensor = tempmnemonic.back();
      tempsubsys = tempmnemonic.substr(0,3);
      tempcharge = tempmnemonic.substr(6,1);
      grifcN = tempadd.substr(tempadd.size()-4,1);
      grif16N = tempadd.substr(tempadd.size()-3,1);
      grif16Ch = tempadd.substr(tempadd.size()-2,2);
      tempgrifc = atoi(grifcN.c_str());
      tempgrif16 = atoi(grif16N.c_str());
      tempgrifch = strtol(grif16Ch.c_str(),NULL,16);
      if(strcmp(experiment.c_str(),"tigress")==0)sprintf(var,"grifadc%i",tigadcNum[tempgrifc][tempgrif16]);
      else if(strcmp(experiment.c_str(),"griffin")==0)sprintf(var,"grifadc%i",grifadcNum[tempgrifc][tempgrif16]);
      if(ignore[0] == 0 && (tempsubsys.compare("TIG") == 0 || tempsubsys.compare("GRG") == 0) && tempoutsensor.compare("A") == 0) continue;  // Optionally Ignores A cores
      if(ignore[1] == 0 && (tempsubsys.compare("TIG") == 0 || tempsubsys.compare("GRG") == 0) && tempoutsensor.compare("B") == 0) continue;  // Optionally Ignores B Cores
      if(ignore[2] == 0 && (tempsubsys.compare("TIS") == 0 || tempsubsys.compare("GRS") == 0)) continue;                                     // Optionally Ignores Suppressors
      if(ignore[3] == 0 && tempsubsys.compare("TIG") == 0 && tempcharge.compare("P") == 0) continue;                                         // Optionally Ignore Segments (TIGRESS Only)
      if(ignore[4] == 0 && tempsubsys.compare("TIG") != 0 && tempsubsys.compare("GRG") != 0 && tempsubsys.compare("TIS") != 0 && tempsubsys.compare("GRS") != 0) continue;                                         // Optionally Ignore Ancillary Detector
      ofile << tempmnemonic << "\t" << var << "\tChannel\t" << tempgrifch + 1 << "\t" << "\n";
    }
  }
  ofile.close();
  cout << "Writing histograms to " << outfile << endl;
  TFile * myfile = new TFile(outfile, "RECREATE");
  myfile->cd();
  list->Write();
  myfile->Close();

}

int main(int argc, char ** argv) {

  char const * ffile;
  char const * conffile;
  char const * outfile;
  char const * calfile;
  printf("Starting sortcode\n");

  if (argc < 3) {
    cout << "Insufficient arguments, provide fragment tree and Configuration files" << endl;
    return 0;
  } else if (argc == 3) {
    ffile = argv[1];
    conffile = argv[2];
    calfile = "CalibrationFile.cal";
    outfile = "channels.root";
  } else if (argc == 4) {
    ffile = argv[1];
    conffile = argv[2];
    calfile = argv[3];
    outfile = "channels.root";
  } else if (argc == 5) {
    ffile = argv[1];
    conffile = argv[2];
    calfile = argv[3];
    outfile = argv[4];
  } else if (argc > 5) {
    printf("Too many arguments\n");
    return 0;
  }

  printf("Input file:%s \nConfiguration file: %s\nCalibration file: %s\nOutput file: %s\n", ffile, conffile, calfile, outfile);

  chan_check(ffile, conffile,calfile, outfile);

  return 0;
}
