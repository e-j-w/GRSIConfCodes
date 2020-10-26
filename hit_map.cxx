//g++ hit_map.cxx -std=c++0x -I$GRSISYS/include -L$GRSISYS/lib `grsi-config --cflags --all-libs --GRSIData-libs` -I$GRSISYS/GRSIData/include -L$GRSISYS/GRSIData/lib `root-config --cflags --libs` -lTreePlayer -lMathMore -lSpectrum -lMinuit -lPyROOT -o HitMap

// This Draws a hit pattern of S3 Detectors in the standard Bambino Configuratuion (1 upstream and 1 downsteam) 
#include <iostream>
#include <iomanip>
#include "TCutG.h"
#include "TH1.h"
#include "TF1.h"
#include "TTree.h"
#include "TChain.h"
#include "TH2.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TDirectory.h"
#include "TList.h"
#include "TRandom.h"
#include "TS3.h"
#include "TGRSIDetectorHit.h"
#include "TGRSIDetectorInformation.h"
#include "TSpectrum.h"
#include "TChannel.h"
#include "TParserLibrary.h"
#include "TEnv.h"
#include "TVector.h"

using namespace std;

double pi = TMath::Pi();

using namespace std;
void HitMap(char const* infile, float phioffset, char const* calfile, char const* outfile){

  TList * list = new TList;
  TH2F *hitmap_down = new TH2F("hitmap_down","Hitmap Downstream",200,-50,50,200,-50,50);list->Add(hitmap_down);
  TH2F *hitmap_up = new TH2F("hitmap_up","Hitmap Upstream",200,-50,50,200,-50,50);list->Add(hitmap_up);

  TFile * inputfile = new TFile(infile, "READ");
  if (!inputfile->IsOpen()) {
    printf("Opening file failed, aborting\n");
    return;
  }
  TChain * AnalysisTree = (TChain * ) inputfile->Get("AnalysisTree");
  printf("%i tree files, details:\n", AnalysisTree->GetNtrees());
  AnalysisTree->ls();
  TTree * tree = (TTree * ) AnalysisTree->GetTree();
  printf("Reading calibration file: %s\n", calfile);
  TChannel::ReadCalFile(calfile);
  Int_t nentries = AnalysisTree->GetEntries();

  TS3 *s3 = 0;
  TS3Hit *s3hit;

  if (AnalysisTree->FindBranch("TS3")) {
  AnalysisTree->SetBranchAddress("TS3",&s3);
  } 
  else {
    cout << "No TS3 Branch Found" << endl;
    return;
  }

  for(int jentry = 0; jentry < nentries; jentry++){
    AnalysisTree->GetEntry(jentry); 
    s3->SetFrontBackEnergy(0.0001);
    for(int i=0;i<s3->GetPixelMultiplicity();i++){
      s3hit = s3->GetS3Hit(i);
//      if(s3hit->GetIsDownstream()){
        TVector3 pos = s3hit->GetPosition(0,true);
	hitmap_down->Fill(pos.X(),pos.Y());
//      } else{
//        TVector3 pos = s3hit->GetPosition(0,true);
//	hitmap_up->Fill(pos.X(),pos.Y());
//      }
    }
    if(jentry%10000 == 0) cout << setiosflags(ios::fixed) << 100 * jentry/nentries << " % complete" << "\r" << flush;  
  }

  cout << "Entry " << nentries << " of " << nentries << ", 100% Complete!\n";

  TFile * myfile = new TFile(outfile, "RECREATE");
  myfile->cd();
  list->Write();
  myfile->Close();
}
int main(int argc, char ** argv) {

  char const * afile;
  float offset;
  char const * outfile;
  char const * calfile;
  printf("Starting sortcode\n");

  std::string grsi_path = getenv("GRSISYS"); // Finds the GRSISYS path to be used by other parts of the grsisort code
  if(grsi_path.length() > 0) {
	  grsi_path += "/";
  }
  // Read in grsirc in the GRSISYS directory to set user defined options on grsisort startup
  grsi_path += ".grsirc";
  gEnv->ReadFile(grsi_path.c_str(), kEnvChange);

  // Input-chain-file, output-histogram-file
  if (argc == 1) {
	  cout << "Insufficient arguments, provide argument tree files" << endl;
	  return 0;
  } else if (argc == 2) {
	  afile = argv[1];
          offset = 22.0;
	  calfile = "CalibrationFile.cal";
	  outfile = "HitMap.root";
  } else if (argc == 3) {
	  afile = argv[1];
          offset = atof(argv[2]);
	  calfile = "CalibrationFile.cal";
	  outfile = "HitMap.root";
  } else if (argc == 4) {
	  afile = argv[1];
          offset = atof(argv[2]);
	  calfile = argv[3];
	  outfile = "HitMap.root";
  } else if (argc == 5) {
	  afile = argv[1];
          offset = atof(argv[2]);
	  calfile = argv[3];
	  outfile = argv[4];
  } else if (argc > 5) {
	  printf("Too many arguments\n");
	  return 0;
  }

  printf("Input file: %s\nPhi Offset: %f\nCalibration file: %s\nOutput file: %s\n", afile, offset, calfile, outfile);

  TParserLibrary::Get()->Load();

  HitMap(afile, offset, calfile, outfile);

  return 0;
}
