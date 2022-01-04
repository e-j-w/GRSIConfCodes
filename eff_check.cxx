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
#define NUM_CRY  12 //number of crystals being tested


void eff_check(char const * infile, char const * calfile, char const * outfile) {

  // R-850 Date and initial activity
  int year=2008, month=9, day=15, hour = 12, min = 0, sec = 0;
  double a0 = 35350; // Bq

  /*// R-1185
  int year=2020, month=5, day=15, hour = 12, min = 0, sec = 0;
  double a0 = 39110; // Bq*/
 
  Double_t nBins = 4096;
  Double_t minBin = 0;
  Double_t maxBin = 2048;
  Double_t binwidth = (maxBin-minBin)/nBins;
  TList * list = new TList;
  TH1F *gamma_singles[NUM_CRY+2], *gamma_singlescal[NUM_CRY+2], *gamma_singles1crycal[NUM_CRY+2], *gamma_singlesallcrycal[NUM_CRY+2], *gamma_gamma_time;
  TH2F *gamma_time[NUM_CRY+2], *gamma_gamma_cal, *gammacal_numcry, *gammacal1crycal_numcry[NUM_CRY+2];
  char hname[20];
  char hname2[20];
  for (int iii = 1; iii < NUM_CRY+2; iii++) {
    sprintf(hname, "h%d", iii);
    gamma_singles[iii] = new TH1F(hname, Form("Gamma Singles Crystal %1.1i", iii), nBins, minBin, maxBin);
    sprintf(hname, "hc%d", iii);
    gamma_singlescal[iii] = new TH1F(hname, Form("Calibrated Gamma Singles Crystal %1.1i", iii), nBins, minBin, maxBin);
    sprintf(hname, "h1c%d", iii);
    gamma_singles1crycal[iii] = new TH1F(hname, Form("Calibrated Gamma Singles Crystal %1.1i (single crystal hit)", iii), nBins, minBin, maxBin);
    sprintf(hname, "hac%d", iii);
    gamma_singlesallcrycal[iii] = new TH1F(hname, Form("Calibrated Gamma Singles Crystal %1.1i (all crystal sum)", iii), nBins, minBin, maxBin);
    sprintf(hname2, "h1cnc%d", iii);
    gammacal1crycal_numcry[iii] = new TH2F(hname2, Form("Calibrated Gamma Singles Crystal %1.1i No Summing vs. Num Crystals", iii), nBins, minBin, maxBin, 4, 1, 4);
    gammacal1crycal_numcry[iii]->GetXaxis()->SetTitle("Single crystal energy (keV)");
  	gammacal1crycal_numcry[iii]->GetYaxis()->SetTitle("Num crystals hit");
    sprintf(hname2, "t%d", iii);
    gamma_time[iii] = new TH2F(hname2, Form("Gamma Energy vs. Time Crystal %1.1i", iii), 1800, 0, 1800, nBins, minBin, maxBin);
  }
  
  gamma_gamma_time = new TH1F("Gamma Time vs. Gamma Time","Gamma Time vs. Gamma Time", 4096, -4096, 4096);
  gamma_gamma_cal = new TH2F("Calibrated Gamma Energy vs. Gamma Energy", "Gamma Energy vs. Gamma Energy", nBins, minBin, maxBin, nBins, minBin, maxBin);
  gammacal_numcry = new TH2F("Calibrated Gamma Energy vs. Num Crystals", "Calibrated Gamma Energy vs. Num Crystals", nBins, minBin, maxBin, 4, 1, 4);
  gammacal_numcry->GetXaxis()->SetTitle("Summed energy (keV)");
  gammacal_numcry->GetYaxis()->SetTitle("Num crystals hit");

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
  TDetectorHit *gen_hit, *gen_hit2;
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
				if(gen_hit->GetDetector()<1 || gen_hit->GetDetector()>NUM_CRY+1)continue;
        gamma_singles[gen_hit->GetDetector()]->Fill(gen_hit->GetCharge());
        gamma_time[gen_hit->GetDetector()]->Fill(gen_hit->GetTime()/1e9,gen_hit->GetCharge());
        for(int j = i+1; j < gen->GetMultiplicity(); j++){
        	gen_hit2 = gen->GetHit(j);
        	gamma_gamma_time->Fill(gen_hit->GetTime()-gen_hit2->GetTime());
        }
    }

    if (jentry % 10000 == 0) cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << nentries << ", " << 100 * jentry / nentries << "% complete" << "\r" << flush;
  }

  for (int iii = 1; iii < NUM_CRY+1; iii++) {
    gamma_singles[iii]->SetBinContent(1, 0);
    if (gamma_singles[iii]->Integral(0, nBins) > 400) {
      list->Add(gamma_singles[iii]);
      list->Add(gamma_singlescal[iii]);
      list->Add(gamma_singles1crycal[iii]);
      list->Add(gamma_singlesallcrycal[iii]);
      list->Add(gammacal1crycal_numcry[iii]);
      list->Add(gamma_time[iii]);
    }
  }
  list->Add(gamma_gamma_time);
  list->Add(gamma_gamma_cal);
  list->Add(gammacal_numcry);

  cout << "Entry " << nentries << " of " << nentries << ", 100% Complete!\n";

  TRunInfo* runinfo = new TRunInfo();
  double runstart = runinfo->RunStart();
  double runstop = runinfo->RunStop();
  double runmid = (runstart+runstop)/2;
  double runlength = runinfo->RunLength();
  //double runlength = 1800.;

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
  cout << "Run performed " << days << " days since source activity measured." << endl;
  double sourceActivity = a0 * exp(-0.693 * days/365.25/5.2714); //A = a0*exp(-lamda dt)

  char colour[NUM_CRY+1][256] = {"","Blue","Green","Red","White","cry5","cry6","cry7","cry8","cry9","cry10","cry11","cry12"};
  double calib_en[2] = {1173.228,1332.492};
  double calib_err[2] = {0.003,0.004};
  double fwhm[NUM_CRY+1];
  double efficiency[NUM_CRY+1];
  double cal_par[NUM_CRY+1][2];
  TF1 * gp[2][NUM_CRY+1];
  int npeaks = 10;

  for (int pp = 1; pp < NUM_CRY+1; pp++) {
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
        //Germanium detector efficiencies are measured relative to that 
        //of a 3" x 3" NaI crystal which has a absolute full energy peak efficiency 
        //of 0.0012 at 25 cm for 1332 keV.
        //Then crystal efficiency expressed as a percentage = (yield (c/s) / source strength (Bq) / 0.0012)*100 .
        efficiency[pp] = (1/binwidth) * gp[p][pp]->GetParameter(0)/ (runlength * sourceActivity * 0.0012) * 100;
      }
    }
    if(nfound == 2){
    	//store cal parameters
    	cal_par[pp][1] = (calib_en[1]-calib_en[0]) / abs(gp[1][pp]->GetParameter(1) - gp[0][pp]->GetParameter(1));
    	if(gp[1][pp]->GetParameter(1) > gp[0][pp]->GetParameter(1)){
    		cal_par[pp][0] = calib_en[1] - cal_par[pp][1]*gp[1][pp]->GetParameter(1);
    	}else{
    		cal_par[pp][0] = calib_en[1] - cal_par[pp][1]*gp[0][pp]->GetParameter(1);
    	}
    }else{
    	cal_par[pp][0] = 0;
    	cal_par[pp][1] = 1;
    }
  }
  
  float allCryCharge;
  for (int jentry = 0; jentry < (nentries - 1); jentry++) {
    tree->GetEntry(jentry);
    allCryCharge = 0.;
    if(gen->GetMultiplicity()==1){
    	gamma_singles1crycal[gen->GetHit(0)->GetDetector()]->Fill(gen->GetHit(0)->GetCharge()*cal_par[gen->GetHit(0)->GetDetector()][1] + cal_par[gen->GetHit(0)->GetDetector()][0]);
    }
    for (int i = 0; i < gen->GetMultiplicity(); i++) {
      gen_hit = gen->GetHit(i);
			if(gen_hit->GetDetector()<1 || gen_hit->GetDetector()>NUM_CRY+1)
				continue;
      gamma_singlescal[gen_hit->GetDetector()]->Fill(gen_hit->GetCharge()*cal_par[gen_hit->GetDetector()][1] + cal_par[gen_hit->GetDetector()][0]);
      gammacal1crycal_numcry[gen_hit->GetDetector()]->Fill(gen_hit->GetCharge()*cal_par[gen_hit->GetDetector()][1] + cal_par[gen_hit->GetDetector()][0],gen->GetMultiplicity());
			if(((gen_hit->GetTime() - gen->GetHit(0)->GetTime()) < 100)&&((gen_hit->GetTime() - gen->GetHit(0)->GetTime()) > -100)){
				allCryCharge += gen_hit->GetCharge()*cal_par[gen_hit->GetDetector()][1] + cal_par[gen_hit->GetDetector()][0];
			}
			for(int j = i+1; j < gen->GetMultiplicity(); j++){
      	gen_hit2 = gen->GetHit(j);
      	gamma_gamma_cal->Fill(gen_hit2->GetCharge()*cal_par[gen_hit->GetDetector()][1] + cal_par[gen_hit->GetDetector()][0],gen_hit->GetCharge()*cal_par[gen_hit->GetDetector()][1] + cal_par[gen_hit->GetDetector()][0]);
      	gamma_gamma_cal->Fill(gen_hit->GetCharge()*cal_par[gen_hit->GetDetector()][1] + cal_par[gen_hit->GetDetector()][0],gen_hit2->GetCharge()*cal_par[gen_hit->GetDetector()][1] + cal_par[gen_hit->GetDetector()][0]); //symmetrized
      }
    }
    gamma_singlesallcrycal[gen->GetHit(0)->GetDetector()]->Fill(allCryCharge);
    gammacal_numcry->Fill(allCryCharge,gen->GetMultiplicity());

    if (jentry % 10000 == 0) cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << nentries << ", " << 100 * jentry / nentries << "% complete" << "\r" << flush;
  }

  
  cout << "\nResolution Measurement using 60Co (R-850) and tiglab DAQ" << "\n";
  //cout << "\nResolution Measurement using 60Co (R-1185) and tiglab DAQ" << "\n";
  cout << setprecision(0) << "Measurement Time = \t " << runlength << " s\n" << "\n";
  cout << "Resolution (FWHM) @ 1332.5 keV" << "\n";
  for(int i = 1; i < NUM_CRY+1; i++)cout <<  setprecision(2) << colour[i] << " Crystal:\t" << fwhm[i] << " keV" << "\n"; 
  cout << "\n\nCorrected Relative Efficiency\n";
  for(int i = 1; i < NUM_CRY+1; i++)cout << setprecision(2) << colour[i] << " Crystal:\t" << efficiency[i] << "%" << "\n"; 


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
	  cout << "Insufficient arguments." << endl;
    cout << "./EffCheck analysis_tree cal_file out_file" << endl;
    cout << "Default cal_file: LabCalFile.cal" << endl;
    cout << "Default out_file: Efficiency.root" << endl;
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
