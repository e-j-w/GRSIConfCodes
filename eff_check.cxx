//g++ eff_check.cxx -std=c++0x -I$GRSISYS/include -L$GRSISYS/lib `grsi-config --cflags --all-libs --GRSIData-libs` -I$GRSISYS/GRSIData/include -L$GRSISYS/GRSIData/lib `root-config --cflags --libs` -lTreePlayer -lMathMore -lSpectrum -lMinuit -lPyROOT -o EffCheck

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
#include "TGenericDetector.h"
#include "TGRSIDetectorHit.h"
#include "TGRSIDetectorInformation.h"
#include "TSpectrum.h"
#include "TChannel.h"
#include "TRunInfo.h"
#include "TParserLibrary.h"
#include "TEnv.h"

using namespace std;


void eff_check(char const * infile, char const * calfile, char const * outfile) {

  // R-850 Date and initial activity
  int year=2008, month=9, day=15, hour = 12, min = 0, sec = 0;
  double a0 = 35350; // Bq
 
  Double_t nBins = 4096;
  Double_t minBin = 0;
  Double_t maxBin = 2048;
  Double_t binwidth = (maxBin-minBin)/nBins;
  TList * list = new TList;
  TH1F * gamma_singles[6];
  TH2F * gamma_time[6];
  char hname[20];
  char hname2[20];
  for (int iii = 1; iii < 6; iii++) {
    sprintf(hname, "h%d", iii);
    gamma_singles[iii] = new TH1F(hname, Form("Gamma Singles Crystal %1.1i", iii), nBins, minBin, maxBin);
    sprintf(hname2, "t%d", iii);
    gamma_time[iii] = new TH2F(hname2, Form("Gamma Energy vs. Time Crystal %1.1i", iii), 1800, 0, 1800, nBins, minBin, maxBin);
  }

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

  TGenericDetector * gen = 0;
  TDetectorHit * gen_hit;
  if (AnalysisTree->FindBranch("TGenericDetector")) {
    AnalysisTree->SetBranchAddress("TGenericDetector", & gen);
  } else {
    cout << "No TGenericDetector Branch Found" << endl;
    return;
  }

  printf("Begin sort\n");
  for (int jentry = 0; jentry < (nentries - 1); jentry++) {
    tree->GetEntry(jentry);
    for (int i = 0; i < gen->GetMultiplicity(); i++) {
        gen_hit = gen->GetHit(i);
	if(gen_hit->GetDetector()<1 || gen_hit->GetDetector()>5)continue;
        gamma_singles[gen_hit->GetDetector()]->Fill(gen_hit->GetCharge());
        gamma_time[gen_hit->GetDetector()]->Fill(gen_hit->GetTime()/1e9,gen_hit->GetCharge());
      }

    if (jentry % 10000 == 0) cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << nentries << ", " << 100 * jentry / nentries << "% complete" << "\r" << flush;
  }

  for (int iii = 1; iii < 5; iii++) {
    gamma_singles[iii]->SetBinContent(1, 0);
    if (gamma_singles[iii]->Integral(0, nBins) > 400) {
      list->Add(gamma_singles[iii]);
      list->Add(gamma_time[iii]);
    }
  }

  cout << "Entry " << nentries << " of " << nentries << ", 100% Complete!\n";

  TRunInfo* runinfo = new TRunInfo();
  double runstart = runinfo->RunStart();
  double runstop = runinfo->RunStop();
  double runmid = (runstart+runstop)/2;
  double runlength = runinfo->RunLength();

  time_t rawtime;
  struct tm * timeinfo;

  // get current timeinfo:
  time ( &rawtime ); //or: rawtime = time(0);
  timeinfo = localtime ( &rawtime ); 
  timeinfo->tm_year   = year - 1900;
  timeinfo->tm_mon    = month - 1;    //months since January - [0,11]
  timeinfo->tm_mday   = day;          //day of the month - [1,31] 
  timeinfo->tm_hour   = hour;         //hours since midnight - [0,23]
  timeinfo->tm_min    = min;          //minutes after the hour - [0,59]
  timeinfo->tm_sec    = sec;          //seconds after the minute - [0,59]

  double sourcedate = mktime ( timeinfo );
  double days = (runmid - sourcedate) / 86400;
  double sourceActivity = a0 * exp(-0.693 * days/365.25/5.2714); //A = a0*exp(-lamda dt)

  char colour[5][20] = {"","Blue","Green","Red","White"};
  double calib_en[2] = {1173.228,1332.492};
  double calib_err[2] = {0.003,0.004};
  double fwhm[5];
  double efficiency[5];
  TF1 * gp[2][5];
  int npeaks = 10;

  for (int pp = 1; pp < 5; pp++) {
    if (gamma_singles[pp]->Integral(0, nBins) < 400) continue;
    gamma_singles[pp]->GetXaxis()->SetRange(200, nBins);
    TSpectrum * s = new TSpectrum(2 * npeaks);
    Int_t nfound = s->Search(gamma_singles[pp], 2, "", 0.7);
    Double_t * xpeaks = s->GetPositionX();
    sort(xpeaks, xpeaks + nfound);
    for (int p = 0; p < nfound; p++) {
      Double_t xp = xpeaks[p];
      Int_t bin = gamma_singles[pp]->GetXaxis()->FindBin(xp);
      Double_t yp = gamma_singles[pp]->GetBinContent(bin);
      gp[p][pp] = new TF1("gp", "gausn(0) + [3] + [4]*(0.5*(1-(ROOT::Math::erf(((x-[1])/([2]*2^(0.5)))))))", xp - 20, xp + 20);
      gp[p][pp]->SetParameters(yp, xp, 1, 10, 10, -1);
      gp[p][pp]->SetParLimits(0, 10, 1e9);
      gp[p][pp]->SetParLimits(1, xp - 10, xp + 10);
      gp[p][pp]->SetParLimits(2, 0.5, 15);
      gp[p][pp]->SetParLimits(4, 1, 10000);
      gamma_singles[pp]->Fit(gp[p][pp], "QR+");
      if(p==1) {
        fwhm[pp] = 2.35 * (gp[p][pp]->GetParameter(2)/gp[p][pp]->GetParameter(1))*calib_en[p];
        efficiency[pp] = (1/binwidth) * gp[p][pp]->GetParameter(0)/ (runlength * sourceActivity * 0.0012) * 100;
      }   
    }      
  }
  
  cout << "\nResolution Measurement using 60Co (R-850) and tiglab DAQ" << "\n";
  cout << setprecision(0) << "Measurement Time = \t " << runlength << " s\n" << "\n";
  cout << "Resolution (FWHM) @ 1332.5 keV" << "\n";
  for(int i = 1; i < 5; i++)cout <<  setprecision(2) << colour[i] << " Crystal:\t" << fwhm[i] << " keV" << "\n"; 
  cout << "\n\nCorrected Relative Efficiency\n";
  for(int i = 1; i < 5; i++)cout << setprecision(2) << colour[i] << " Crystal:\t" << efficiency[i] << "%" << "\n"; 


  cout << "Writing histograms to " << outfile << endl;
  TFile * myfile = new TFile(outfile, "RECREATE");
  myfile->cd();
  list->Write();
  myfile->Close();
  return;
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

  // Input-chain-file, output-histogram-file
  if (argc == 1) {
	  cout << "Insufficient arguments, provide argument tree files" << endl;
	  return 0;
  } else if (argc == 2) {
	  afile = argv[1];
	  calfile = "LabCalFile.cal";
	  outfile = "Efficiency.root";
  } else if (argc == 3) {
	  afile = argv[1];
	  calfile = argv[2];
	  outfile = "Efficiency.root";
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

  eff_check(afile, calfile, outfile);

  return 0;
}
