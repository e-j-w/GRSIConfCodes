//file path to default calibration file
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <dirent.h>

Double_t num_bins = 8192; //number of bins to read in from histograms
Double_t min_bin = 0; //first bin in histogram
Double_t max_bin = 4096; //last bin in histogram
Int_t num_cores = 64;
bool multirun = true; // Uses 1 subrun or all available

//function to load histograms and add them to TLists - Uses Multiple Subruns
void load_sum_histograms(const char analysis_filepath[], char cal_filepath[], TH1F *hist[], Int_t source_count, Int_t num_cores) {

  TChain *analysis = new TChain("AnalysisTree");
  string infile = analysis_filepath;
  int size = infile.find_last_of("/");
  string directory = infile.substr(0,size);
  int num = infile.find_last_of("_");
  string runnumber = infile.substr(num-5,5);

  //Finds all subruns for a specific run
  vector<string> runlist;
  DIR * pDIR;
  struct dirent * entry;
  if ((pDIR = opendir(directory.c_str()))) {
    while ((entry = readdir(pDIR))) {
      if (strstr(entry->d_name, runnumber.c_str())) {
        string file = directory;
	file.append("/");
	file.append(entry->d_name);
	if(strstr(file.c_str(),"analysis")) runlist.push_back(file);
      }
    }
    closedir(pDIR);
  }

  sort(runlist.begin(),runlist.end()); // Puts subruns in order
  for(int i = 0; i<runlist.size(); i++) {
    analysis->Add(runlist.at(i).c_str());
  }

  printf("%i tree files, details:\n", analysis->GetNtrees());
  TTree *tree = (TTree *) analysis->GetTree();

  //opening calibration file that organizes data into channels
  TChannel::ReadCalFile(cal_filepath);

  Int_t num_entries = analysis->GetEntries();
  TTigress * tigress = 0;
  TGriffin * griffin = 0;
  bool tig = true;
  if (analysis->FindBranch("TTigress")) {
    analysis->SetBranchAddress("TTigress", & tigress);
    cout << "Creating TIGRESS Graphs" << endl;
    tig = true;
  } else {
    if (analysis->FindBranch("TGriffin")) {
      analysis->SetBranchAddress("TGriffin", & griffin);
      cout << "Creating GRIFFIN Graphs" << endl;
      tig = false;
    } else {
      cout << "No TTigress or TGriffin Branch Found Things will go wrong" << endl;
    }
  }
  char hname[20];
  for(int i = 0; i < num_cores; i++) {
    sprintf(hname,"hist%i_%i",source_count, i);
    hist[i] = new TH1F(hname, hname, num_bins, min_bin, max_bin);
  }

  cout << "Histograms created." << endl;

  //filling histograms with data from analysis root file
  for (int i = 0; i < num_entries - 1; i++) {
    analysis->GetEntry(i);
    if(tig) {
      for (int j = 0; j < tigress->GetMultiplicity(); j++) {
        TTigressHit *tigress_hit = tigress->GetTigressHit(j);
        hist[tigress_hit->GetArrayNumber()]->Fill(tigress_hit->GetCharge());
      }//for
    } else {
      for (int j = 0; j < griffin->GetMultiplicity(); j++) {
        TGriffinHit *griffin_hit = griffin->GetGriffinHit(j);
        hist[griffin_hit->GetArrayNumber()-1]->Fill(griffin_hit->GetCharge());
      }//for
    }
    if (i % 10000 == 0) {
      cout << setiosflags(ios::fixed) << "Entry " << i << " of " << num_entries << ", " << 100 * (i) / num_entries << "% complete" << "\r" << flush;
    }//if
  }//for
  cout << setiosflags(ios::fixed) << "Entry " << num_entries << " of " << num_entries << ", 100% complete" << "\r" << flush;	
}//load_sum_histograms

//Loads histogram from one subrun only
void load_histograms(const char analysis_filepath[], char cal_filepath[], TH1F *hist[], Int_t source_count, Int_t num_cores) {

  //opening analysis root file
  TFile *input_file = new TFile(analysis_filepath, "READ");
  if (!input_file->IsOpen()) {
    cout << "Cannot open input file!" << endl;
    return;
  }//if
  TChain *analysis = (TChain *) input_file->Get("AnalysisTree");
  TTree *tree = (TTree *) analysis->GetTree();

  //opening calibration file that organizes data into channels
  TChannel::ReadCalFile(cal_filepath);

  Int_t num_entries = analysis->GetEntries();
  TTigress * tigress = 0;
  TGriffin * griffin = 0;
  bool tig = true;
  if (analysis->FindBranch("TTigress")) {
    analysis->SetBranchAddress("TTigress", & tigress);
    cout << "Creating TIGRESS Graphs" << endl;
    tig = true;
  } else {
    if (analysis->FindBranch("TGriffin")) {
      analysis->SetBranchAddress("TGriffin", & griffin);
      cout << "Creating GRIFFIN Graphs" << endl;
      tig = false;
    } else {
      cout << "No TTigress or TGriffin Branch Found Things will go wrong" << endl;
    }
  }

  char hname[20];
  for(int i = 0; i < num_cores; i++) {
    sprintf(hname,"hist%i_%i",source_count, i);
    hist[i] = new TH1F(hname, hname, num_bins, min_bin, max_bin);
  }

  cout << "Histograms created." << endl;

  //filling histograms with data from analysis root file
  for (int i = 0; i < num_entries - 1; i++) {
    tree->GetEntry(i);
    if(tig) {
      for (int j = 0; j < tigress->GetMultiplicity(); j++) {
        TTigressHit *tigress_hit = tigress->GetTigressHit(j);
        hist[tigress_hit->GetArrayNumber()]->Fill(tigress_hit->GetCharge());
      }//for
    } else {
      for (int j = 0; j < griffin->GetMultiplicity(); j++) {
        TGriffinHit *griffin_hit = griffin->GetGriffinHit(j);
        hist[griffin_hit->GetArrayNumber()-1]->Fill(griffin_hit->GetCharge());
      }//for
    }
    if (i % 10000 == 0) {
      cout << setiosflags(ios::fixed) << "Entry " << i << " of " << num_entries << ", " << 100 * (i) / num_entries << "% complete" << "\r" << flush;
    }//if
  }//for
  cout << setiosflags(ios::fixed) << "Entry " << num_entries << " of " << num_entries << ", 100% complete" << "\r" << flush;

}
//main function
void make_ge_histograms() {

  TList *list = new TList;

  Int_t num_sources = 0; //number of sources used in this calibration
  string source[5];
  string rootfile[5];
  int iii = 0;
  string sourcefile = "sourcelist.dat";
  ifstream fp;
  fp.open(sourcefile.c_str());
  while (fp.good()) {
    fp >> source[iii] >> rootfile[iii];
    iii++;
  }
  char empty_gains_cal[70] = "CalibrationFile.cal";
  num_sources = iii-1;
  if (num_sources == 0) {
    cout << "No sources inputted. Exiting program..." << endl;
    return;
  }//

  TH1F *cal_hist[num_sources][num_cores];
  for (int k = 0; k < num_sources; k++) {
    if(multirun)load_sum_histograms(rootfile[k].c_str(), empty_gains_cal, cal_hist[k], k, num_cores);
    else load_histograms(rootfile[k].c_str(), empty_gains_cal, cal_hist[k], k, num_cores);
    for (int j = 0; j < num_cores; j++) list->Add(cal_hist[k][j]);
  }
  //writes out the histograms in the TList
  TFile *out = new TFile("CalHist.root", "RECREATE");
  out->cd();
  list->Write();
  out->Close();//closes files

}//histogram_make
