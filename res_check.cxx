//g++ res_check.cxx -std=c++0x -I$GRSISYS/include -L$GRSISYS/lib `grsi-config --cflags --all-libs --GRSIData-libs` -I$GRSISYS/GRSIData/include -L$GRSISYS/GRSIData/lib `root-config --cflags --libs` -lTreePlayer -lMathMore -lSpectrum -lMinuit -lPyROOT -o ResCheck

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
#include "TTigress.h"
#include "TGriffin.h"
#include "TGRSIDetectorInformation.h"
#include "TSpectrum.h"
#include "TChannel.h"
#include "TParserLibrary.h"
#include "TEnv.h"

using namespace std;

void res_check(char const * infile, char const * calfile, char const * outfile) {

  Double_t nBins = 4096;
  Double_t min = 0;
  Double_t max = 2048;

  TList * list = new TList;

  TH1F * gamma_singles[2000];
  char hname[20];
  for (int iii = 0; iii < 65; iii++) {
    sprintf(hname, "h%d", iii+1);
    gamma_singles[iii] = new TH1F(hname, Form("Gamma Singles Crystal %1.1i", iii+1), nBins, min, max);
  } 
  TH2F * gA = new TH2F("gA","",64,0,64,nBins,min,max);list->Add(gA);

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

  int npeaks = 10;
  TTigress * tigress = 0;
  TGriffin * griffin = 0;
  TTigressHit * tigress_hit;
  TGriffinHit * griffin_hit;
  bool tig = true;
  if (AnalysisTree->FindBranch("TTigress")) {
    AnalysisTree->SetBranchAddress("TTigress", & tigress);
     cout << "Checking TIGRESS" << endl;
     tig = true;
  } else {
  if (AnalysisTree->FindBranch("TGriffin")) {
    AnalysisTree->SetBranchAddress("TGriffin", & griffin);
    cout << "Checking GRIFFIN" << endl;
    tig = false;
  } else {
    cout << "No TTigress or TGriffin Branch Found Things will go wrong" << endl;
    }
  }

  printf("Begin sort\n");
  int one;
  for (int jentry = 0; jentry < (nentries - 1); jentry++) {
    tree->GetEntry(jentry);
    if(tig) {
      for (one = 0; one < tigress->GetMultiplicity(); one++) {
        tigress_hit = tigress->GetTigressHit(one);
        gamma_singles[tigress_hit->GetArrayNumber()]->Fill(tigress_hit->GetCharge());
        gA->Fill(tigress_hit->GetArrayNumber(),tigress_hit->GetCharge());
      }
    }
    else {
      for (one = 0; one < griffin->GetMultiplicity(); one++) {
        griffin_hit = griffin->GetGriffinHit(one);
        gamma_singles[griffin_hit->GetArrayNumber() - 1]->Fill(griffin_hit->GetCharge());
        gA->Fill(griffin_hit->GetArrayNumber() - 1, griffin_hit->GetCharge());
      }  
    }
    if (jentry % 10000 == 0) cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << nentries << ", " << 100 * jentry / nentries << "% complete" << "\r" << flush;
  }

  for (int iii = 0; iii < 65; iii++) {
    gamma_singles[iii]->SetBinContent(1, 0);
    if (gamma_singles[iii]->Integral(0, nBins) > 400) list->Add(gamma_singles[iii]);
  }

  cout << "Entry " << nentries << " of " << nentries << ", 100% Complete!\n";

  cout << "Histograms written, sorting complete" << endl;

  double calib_en[2] = {1173.228,1332.492};
  double calib_err[2] = {0.003,0.004};
  double array[64][64];
  double earray[64][64];
  TF1 * gp[64][64];
  ofstream ofile;
  ofile.open("Ge_Resolution.txt");
  if (!ofile.is_open()) {
    printf("Output not opened\n");
    return;
  }

  float avgRes = 0.0f;
  float minRes = 100000.0f;
  float maxRes = 0.0f;
  float res[64];
  int numCry = 0;

  for(int pp = 0; pp < 64; pp++){
    res[pp] = -1.0f;
    int minBin = 1000; //set the minimum bin, to cut out low energy stuff
    if (gamma_singles[pp]->Integral(minBin, nBins) < 400) continue;
	  for(int i=0; i<minBin; i++){
      gamma_singles[pp]->SetBinContent(i,0.0);
    }
    gamma_singles[pp]->GetXaxis()->SetRange(minBin, nBins);
    TSpectrum * s = new TSpectrum(2 * npeaks);
    Int_t nfound = s->Search(gamma_singles[pp], 2, "", 0.3);
    Double_t * xpeaks = s->GetPositionX();
    sort(xpeaks, xpeaks + nfound);
    for (int p = 0; p < nfound; p++) {
      Double_t xp = xpeaks[p];
      Int_t bin = gamma_singles[pp]->GetXaxis()->FindBin(xp);
      Double_t yp = gamma_singles[pp]->GetBinContent(bin);
      gp[p][pp] = new TF1("gp", "[0]*(exp(-((x-[1])^2/(2*[2]^2)))) + [3] + [4]*(0.5*(1-(ROOT::Math::erf(((x-[1])/([2]*2^(0.5)))))))", xp - 20, xp + 20);
      gp[p][pp]->SetParameters(yp, xp, 1, 10, 10, -1);
      gp[p][pp]->SetParLimits(0, 10, 1e9);
      gp[p][pp]->SetParLimits(1, xp - 10, xp + 10);
      gp[p][pp]->SetParLimits(2, 0.5, 15);
      gp[p][pp]->SetParLimits(4, 1, 10000);
      gamma_singles[pp]->Fit(gp[p][pp], "QR+");
      array[p][pp] = gp[p][pp]->GetParameter(1);
      earray[p][pp] = gp[p][pp]->GetParError(1);
      if(p==1){ //1332 keV peak
        res[pp] = (gp[p][pp]->GetParameter(2)/gp[p][pp]->GetParameter(1))*calib_en[p]*2.35482;
        ofile << pp + 1 << ",\t" << gp[p][pp]->GetParameter(1) << ",\t" << gp[p][pp]->GetParameter(2) << ",\t" << res[pp] << endl;
        if(res[pp] >= 3.0) printf( DRED "ArrayNumber:\t%i\tCentroid:\t%f\tResolution (FWHM):\t%f keV\t\n" RESET_COLOR "",pp + 1, gp[p][pp]->GetParameter(1),res[pp]);  
        else printf("ArrayNumber:\t%i\tCentroid:\t%f\tResolution (FWHM):\t%f keV\t\n",pp + 1, gp[p][pp]->GetParameter(1),res[pp]);
        avgRes += res[pp];
        if(res[pp] < minRes){
          minRes = res[pp];
        }
        if(res[pp] > maxRes){
          maxRes = res[pp];
        }
        numCry++;
      }       
    }      
  }

  avgRes /= numCry*1.0f;
  float stdevRes = 0.0f;
  for(int pp = 0; pp < 64; pp++){
    if(res[pp] >= 0.0f){
      stdevRes += (res[pp] - avgRes)*(res[pp] - avgRes);
    }
  }
  stdevRes = sqrtf(stdevRes/(numCry + 1.0f));

  if(avgRes >= 2.8) printf( DRED "Average Resolution (FWHM):\t%f keV\t\n" RESET_COLOR "",avgRes);  
  else printf("Average Resolution (FWHM):\t%f keV\t\n",avgRes);

  if(stdevRes >= 1.0) printf( DRED "   Standard Deviation:\t%f keV\t\n" RESET_COLOR "",stdevRes);  
  else printf("   Standard Deviation:\t\t%f keV\t\n",stdevRes);
  
  if(minRes >= 2.5) printf( DRED "Minimum (Worst) Resolution (FWHM):\t%f keV\t\n" RESET_COLOR "",maxRes);  
  else printf("Minimum (Best) Resolution (FWHM):\t%f keV\t\n",minRes);

  if(maxRes >= 3.2) printf( DRED "Maximum (Worst) Resolution (FWHM):\t%f keV\t\n" RESET_COLOR "",maxRes);  
  else printf("Maximum (Worst) Resolution (FWHM):\t%f keV\t\n",maxRes);

  cout << "Writing histograms to " << outfile << endl;
  TFile * myfile = new TFile(outfile, "RECREATE");
  myfile->cd();
  list->Write();
  myfile->Close();
}

int main(int argc, char ** argv) {

  char const * afile;
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

  cout << "Resolution check code for in-array 60Co data." << endl;
  if (argc == 1) {
	  cout << "Insufficient arguments." << endl;
    cout << "./ResCheck analysis_tree cal_file out_file" << endl;
    cout << "Default cal_file: LabCalFile.cal" << endl;
    cout << "Default out_file: Efficiency.root" << endl;
	  return 0;
  } else if (argc == 2) {
	  afile = argv[1];
	  calfile = "CalibrationFile.cal";
	  outfile = "rCheck.root";
  } else if (argc == 3) {
	  afile = argv[1];
	  calfile = argv[2];
	  outfile = "rCheck.root";
  } else if (argc == 4) {
	  afile = argv[1];
	  calfile = argv[2];
	  outfile = argv[3];
  } else if (argc > 4) {
	  printf("Too many arguments\n");
	  return 0;
  }

  printf("Input file:%s\nCalibration file: %s\nOutput file: %s\n", afile, calfile, outfile);

  TParserLibrary::Get()->Load();

  res_check(afile, calfile, outfile);

  return 0;
}
