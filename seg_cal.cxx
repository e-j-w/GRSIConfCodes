//g++ seg_cal.C -std=c++0x -I$GRSISYS/include -L$GRSISYS/lib `grsi-config --cflags --all-libs --GRSIData-libs` -I$GRSISYS/GRSIData/include -L$GRSISYS/GRSIData/lib `root-config --cflags --libs` -lTreePlayer -lMathMore -lSpectrum -lMinuit -lPyROOT -o SegCal

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
#include "TTip.h"
#include "TTigress.h"
#include "TSharc.h"
#include "TRF.h"
#include "TSpectrum.h"
#include "TChannel.h"
#include "TParserLibrary.h"
#include "TEnv.h"

using namespace std;

void seg_cal(char const * infile, char const * calfile, char const * outfile) {


  Double_t nBins = 2048;
  Double_t min = 0;
  Double_t max = 2048;

  Double_t nBins2 = 2048;
  Double_t min2 = 0;
  Double_t max2 = 2048;

  TList * list = new TList;

  TH2F * gamma_singles[512];
  char hname[20];
  for (int iii = 0; iii < 512; iii++) {
    sprintf(hname, "h%d", iii);
    gamma_singles[iii] = new TH2F(hname, Form("Segment Number %1.1i", iii), nBins2, min2, max2, nBins, min, max);
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

  int npeaks = 10;

  TTigress * tigr = 0;
  AnalysisTree->SetBranchAddress("TTigress", & tigr);
  TTigressHit * tigress_hit;
  //TGRSIDetectorHit *segment_hit;
  printf("Begin sort\n");
  int one;
  for (int jentry = 0; jentry < (nentries - 1); jentry++) {
    tree->GetEntry(jentry);
    for (one = 0; one < tigr->GetMultiplicity(); one++) {
      tigress_hit = tigr->GetTigressHit(one);
      for (int j = 0; j < tigress_hit->GetSegmentMultiplicity(); j++) {
        TGRSIDetectorHit segment_hit = tigress_hit->GetSegmentHit(j);
        if (segment_hit.GetCharge() < 10) continue;
        int seg_index = 8 * (tigress_hit->GetArrayNumber()) + segment_hit.GetSegment() - 1;
        gamma_singles[seg_index]->Fill(segment_hit.GetCharge(), tigress_hit->GetEnergy());
//	cout << j << endl;
      }
    }
    if (jentry % 10000 == 0) cout << setiosflags(ios::fixed) << "TIGRESS Entry " << jentry << " of " << nentries << ", " << 100 * jentry / nentries << "% complete" << "\r" << flush;
  }

  for (int iii = 0; iii < 512; iii++) {
    if (gamma_singles[iii]->Integral(min, max, min2, max2) > 50) list->Add(gamma_singles[iii]);
  }

  cout << "Entry " << nentries << " of " << nentries << ", 100% Complete!\n";

  cout << "Histograms written, sorting complete" << endl;

  cout << "Drawing TGraphs" << endl;

  const char *out = "SegmentParameter.txt";
  ofstream ofile;
  ofile.open(out);
  if (!ofile.is_open()) {
    printf("Output not opened\n");
    return;
  }
  TGraph * gr1[512];
  TF1 * ccurve[512];
  for (int k = 0; k < 512; k++) {
    if (gamma_singles[k]->Integral(min2, max2, min, max) < 50) ofile <<"0\t0" << endl;
    if (gamma_singles[k]->Integral(min2, max2, min, max) < 50) continue;
    int counter = 0;
    gr1[k] = new TGraph();
    gr1[k]->SetTitle(Form("Segment_Index_%i", k));
    gr1[k]->SetName(Form("Segment_Index_%i", k));
    for (int i = 0; i < gamma_singles[k]->GetSize(); i++) {
      int x, y, z;
      gamma_singles[k]->GetBinXYZ(i, x, y, z);
      double X = gamma_singles[k]->GetXaxis()->GetBinCenter(x);
      double Y = gamma_singles[k]->GetYaxis()->GetBinCenter(y);
     if (gamma_singles[k]->GetBinContent(x, y) > 10) {
        if (x < 10 || y < 10) continue;
        if (X > 1500) continue;
        gr1[k]->SetPoint(counter, X, Y);
        counter++;
      }
    }
    gr1[k]->Draw("A*");
    ccurve[k] = new TF1("curve", "pol1");
    gr1[k]->Fit(ccurve[k],"Q");
    list->Add(gr1[k]);
    if (ccurve[k]->GetParameter(1)>0.45 && ccurve[k]->GetParameter(1)<0.55)ofile << ccurve[k]->GetParameter(0) << "\t" << ccurve[k]->GetParameter(1) << endl;
    if (ccurve[k]->GetParameter(1)<0.45 || ccurve[k]->GetParameter(1)>0.55) {
      // Honestly I can't remember why I put this check in, I'm sure there was a reason, but I can't remember. I think it may mean the fit failed or the segment has a problem
      int anum = k/8;
      int segnum = k % 8 + 1;
      cout << "Something has gone wrong? check Array Number " << anum << " Segment" << segnum << endl;
      ofile << ccurve[k]->GetParameter(0) << "\t" << ccurve[k]->GetParameter(1) << endl;
    }

  }
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

  // Input-chain-file, output-histogram-file
  if (argc == 1) {
	  cout << "Insufficient arguments, provide argument tree files" << endl;
	  return 0;
  } else if (argc == 2) {
	  afile = argv[1];
	  calfile = "CalibrationFile.cal";
	  outfile = "segCheck.root";
  } else if (argc == 3) {
	  afile = argv[1];
	  calfile = argv[2];
	  outfile = "segCheck.root";
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

  seg_cal(afile, calfile, outfile);

  return 0;
}
