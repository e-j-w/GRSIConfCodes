#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <dirent.h>

Double_t num_bins = 8192; //number of bins to read in from histograms
Double_t min_bin = 0; //first bin in histogram
Double_t max_bin = 4096; //last bin in histogram
Double_t r2d = 180 / TMath::Pi();
Double_t distance = 145.0 // Detector distance from source 110.0 (high efficiency) or 145.0 (high peak to total)

// Gets Back to back detectors for addback
bool p180[16][16] = {
{false,false,false,false,false,false,false,false,false,false,false,false,false,false,true,false},
{false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,true},
{false,false,false,false,false,false,false,false,false,false,false,false,true,false,false,false},
{false,false,false,false,false,false,false,false,false,false,false,false,false,true,false,false},
{false,false,false,false,false,false,false,false,true,false,false,false,false,false,false,false},
{false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false},
{false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false},
{false,false,false,false,false,false,false,false,false,false,false,true,false,false,false,false},
{false,false,false,false,true,false,false,false,false,false,false,false,false,false,false,false},
{false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false},
{false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false},
{false,false,false,false,false,false,false,true,false,false,false,false,false,false,false,false},
{false,false,true,false,false,false,false,false,false,false,false,false,false,false,false,false},
{false,false,false,true,false,false,false,false,false,false,false,false,false,false,false,false},
{true,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false},
{false,true,false,false,false,false,false,false,false,false,false,false,false,false,false,false}
};

//function to load histograms and add them to TLists
void load_histograms(const char analysis_filepath[], char cal_filepath[], TH1F *singles[], TH1F *addback[], TH2F *singles180[], TH2F *addback180[], Int_t source_count, int &RunStart, int &RunStop) {
  TChain *analysis = new TChain("AnalysisTree");
  string infile = analysis_filepath;
  int size = infile.find_last_of("/");
  string directory = infile.substr(0,size);
  int num = infile.find_last_of("_");
  string runnumber = infile.substr(num-5,5);

  //Finds all subruns for a specific run in the same directory
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

  double runstart = 0;
  double runstop = 0;

  for(int i = 0; i<runlist.size(); i++) {
    TFile *temp_file = new TFile(runlist.at(i).c_str());
    TRunInfo* runinfo = new TRunInfo();
    if(i==0) runstart = runinfo->RunStart();
    if(i==runlist.size()-1) runstop = runinfo->RunStop();
    temp_file->Close();
    analysis->Add(runlist.at(i).c_str());
  }
  RunStop = runstop;
  RunStart = runstart;

  printf("%i tree files, details:\n", analysis->GetNtrees());
  TTree *tree = (TTree *) analysis->GetTree();
  Int_t num_entries = analysis->GetEntries();

  //opening calibration file that organizes data into channels
  TChannel::ReadCalFile(cal_filepath);


  // Checks for any crystals to be exlcluded
  string  badfile = "badcrystal.dat";
  vector<int> ignore;
  string crynum;
  ifstream fp;
  fp.open(badfile.c_str());
  while (getline(fp,crynum)) {
    ignore.push_back(stoi(crynum));
  }

  // Pointerrs for things.
  TTigress * tigress = 0;
  TGriffin * griffin = 0;
  TGriffinBgo * fbgo = 0;

  bool tig = true; // Switch betweem filling TIGRESS or GRIFFIN histogram
  bool suppTIG, suppTIG2; // TIGRESS singles Suppression
  bool suppADD, suppADD2; // TIGRESS Addback Suppression

  // Checks if analysis tree contains TIGRESS or GRIFFIN data
  if (analysis->FindBranch("TTigress")) {
    analysis->SetBranchAddress("TTigress", & tigress);
    cout << "Creating TIGRESS Graphs" << endl;
    tig = true;
  } else {
    if (analysis->FindBranch("TGriffin")) {
      analysis->SetBranchAddress("TGriffin", & griffin);
      cout << "Creating GRIFFIN Graphs" << endl;
      tig = false;
      if(analysis->FindBranch("TGriffinBgo")) {
        analysis->SetBranchAddress("TGriffinBgo", & fbgo);
      }
    } else {
      cout << "No TTigress or TGriffin Branch Found Things will go wrong" << endl;
    }
  }

  // Defines Histograms
  char hname[20];
  sprintf(hname,"singles_%i",source_count);
  char aname[20];
  sprintf(aname,"addback_%i",source_count);
  char hname180[20];
  sprintf(hname180,"singles180_%i",source_count);
  char aname180[20];
  sprintf(aname180,"addback180_%i",source_count);

  singles[0] = new TH1F(hname, "Ge Singles", num_bins, min_bin, max_bin);
  addback[0] = new TH1F(aname, "Ge Addback", num_bins, min_bin, max_bin);
  singles180[0] = new TH2F(hname180, "Ge Singles 180", num_bins, min_bin, max_bin, num_bins, min_bin, max_bin);
  addback180[0] = new TH2F(aname180, "Ge Addback 180", num_bins, min_bin, max_bin, num_bins, min_bin, max_bin);

  cout << "Histograms created." << endl;

  //filling histograms with data from analysis root file
  for (int i = 0; i < num_entries - 1; i++) {
    analysis->GetEntry(i);
    if(tig) { // If TIGRESS data
      for (int j = 0; j < tigress->GetMultiplicity(); j++) {
        TTigressHit *tigress_hit = tigress->GetTigressHit(j);
        suppTIG = tigress_hit->BGOFired();
        if(suppTIG) continue; // Suppression
        if(find(ignore.begin(),ignore.end(),tigress_hit->GetArrayNumber()) != ignore.end()) continue; // Skips crystals in bad-crystal.dat
        singles[0]->Fill(tigress_hit->GetEnergy());
	if(tigress->GetMultiplicity()>1) {
	  for(int k = 0; k < tigress->GetMultiplicity(); k++) {
	    if( j != k ) {
	      TTigressHit *tigress_hit2 = tigress->GetTigressHit(k);
	      suppTIG2 = tigress_hit2->BGOFired();
              if(suppTIG2) continue; // Suppression
              TVector3 pos1 = tigress_hit->GetPosition(distance);
              TVector3 pos2 = tigress_hit2->GetPosition(distance);
//              TVector3 pos1 = tigress_hit->GetCorePosition(distance); // Core Only position if segments present, does not yet exist in GRSISort
//              TVector3 pos2 = tigress_hit2->GetCorePosition(distance); //  Core Only position if segments present, does not yet exist in GRSISort
	      if(abs(tigress_hit->GetTime() - tigress_hit2->GetTime()) > 100) continue;
	      if(find(ignore.begin(),ignore.end(),tigress_hit2->GetArrayNumber()) == ignore.end() ) {
   	        if(pos1.Angle(pos2)*r2d < 182 && pos1.Angle(pos2)*r2d > 178)singles180[0]->Fill(tigress_hit->GetEnergy(),tigress_hit2->GetEnergy());
	      }
 	    }
	  }
	}
      }//for
      for (int j = 0; j < tigress->GetAddbackMultiplicity(); j++) {
        TTigressHit *tigress_add_hit = tigress->GetAddbackHit(j);
        suppADD = tigress_add_hit->BGOFired();
        if(suppADD) continue; // Suppression
        if(find(ignore.begin(),ignore.end(),tigress_add_hit->GetArrayNumber()) != ignore.end()) continue;
	addback[0]->Fill(tigress_add_hit->GetEnergy());
        if(tigress->GetAddbackMultiplicity() > 1) {
          for(int k = 0; k < tigress->GetAddbackMultiplicity(); k++) {
            if( j != k ) {
              TTigressHit *tigress_add_hit2 = tigress->GetAddbackHit(k);
              suppADD2 = tigress_add_hit2->BGOFired();
              if(suppADD2) continue; // Suppression
              if(abs(tigress_add_hit->GetTime() - tigress_add_hit2->GetTime()) > 300 ) continue;
              if(find(ignore.begin(),ignore.end(),tigress_add_hit2->GetArrayNumber()) == ignore.end() ) {
	        if(p180[tigress_add_hit->GetDetector()-1][tigress_add_hit2->GetDetector()-1])addback180[0]->Fill(tigress_add_hit->GetEnergy(),tigress_add_hit2->GetEnergy()); // Back to back clovers for addback
              }
            }
          }
        }
      }  //for
    } else { // If GRIFFIN data
      for (int j = 0; j < griffin->GetSuppressedMultiplicity(fbgo); j++) {
        TGriffinHit *griffin_hit = griffin->GetSuppressedHit(j);
        if(find(ignore.begin(),ignore.end(),griffin_hit->GetArrayNumber()) != ignore.end()) continue;
	singles[0]->Fill(griffin_hit->GetEnergy());
        if(griffin->GetSuppressedMultiplicity(fbgo)>1) {
          for(int k = 0; k < griffin->GetSuppressedMultiplicity(fbgo); k++) {
            if( j != k ) {
              TGriffinHit *griffin_hit2 = griffin->GetSuppressedHit(k);
              if(abs(griffin_hit->GetTime() - griffin_hit2->GetTime()) > 100 ) continue;
              if( find(ignore.begin(),ignore.end(),griffin_hit->GetArrayNumber()) == ignore.end() && find(ignore.begin(),ignore.end(),griffin_hit2->GetArrayNumber()) == ignore.end() ) {
                TVector3 pos1 = griffin_hit->GetPosition(distance);
                TVector3 pos2 = griffin_hit2->GetPosition(distance);
                if(pos1.Angle(pos2)*r2d < 182 && pos1.Angle(pos2)*r2d > 178)singles180[0]->Fill(griffin_hit->GetEnergy(),griffin_hit2->GetEnergy());
              }
            }
          }
        }
      } //for
      for (int j = 0; j < griffin->GetSuppressedAddbackMultiplicity(fbgo); j++) {
        TGriffinHit *griffin_add_hit = griffin->GetSuppressedAddbackHit(j);
        if(find(ignore.begin(),ignore.end(),griffin_add_hit->GetArrayNumber()) != ignore.end()) continue;
	addback[0]->Fill(griffin_add_hit->GetEnergy());
        if(griffin->GetSuppressedAddbackMultiplicity(fbgo)>1) {
          for(int k = 0; k < griffin->GetSuppressedAddbackMultiplicity(fbgo); k++) {
            if( j != k ) {
              TGriffinHit *griffin_add_hit2 = griffin->GetSuppressedAddbackHit(k);
              if(abs(griffin_add_hit->GetTime() - griffin_add_hit2->GetTime()) > 300 ) continue;
              if( find(ignore.begin(),ignore.end(),griffin_add_hit->GetArrayNumber()) == ignore.end() && find(ignore.begin(),ignore.end(),griffin_add_hit2->GetArrayNumber()) == ignore.end() ) {
                if(p180[griffin_add_hit->GetDetector()-1][griffin_add_hit2->GetDetector()-1])addback180[0]->Fill(griffin_add_hit->GetEnergy(),griffin_add_hit2->GetEnergy()); // Back to back clovers for addback
              }
            }
          }
        }
      }  //for
    }
    if (i % 10000 == 0) {
      cout << setiosflags(ios::fixed) << "Entry " << i << " of " << num_entries << ", " << 100 * (i) / num_entries << "% complete" << "\r" << flush;
    }//if
  }//for
  cout << "Entry " << num_entries << " of " << num_entries << ", 100% complete" << endl;

  // Makes sure errors are set correctly
  for(int i = 0; i < num_bins; i++) {
    double bc = singles[0]->GetBinContent(i);
    double bca = addback[0]->GetBinContent(i);
    if(bc==0) {
      singles[0]->SetBinError(i,1);
    } else {
      singles[0]->SetBinError(i, pow(bc,0.5));
    }
    if(bca==0) {
      addback[0]->SetBinError(i,1);
    } else {
      addback[0]->SetBinError(i, pow(bca,0.5));
    }
  }
}//load_histograms

//main function
void histogram_make() {

  TList *list = new TList;
  Int_t num_sources = 0; //number of sources used in this calibration
  string sourceID[10];
  string rootfile[10];
  int iii = 0;
  // Reads in run information
  string sourcefile = "runlist.dat";
  ifstream fp;
  fp.open(sourcefile.c_str());
  while (fp.good()) {
    fp >> sourceID[iii] >> rootfile[iii];
    iii++;
  }
  char empty_gains_cal[70] = "CalibrationFile.cal";
  num_sources = iii-1;
  if (num_sources == 0) {
    cout << "No sources inputted. Exiting program..." << endl;
    return;
  }//

  // Declare histograms
  TH1F *singles_hist[num_sources][1];
  TH1F *addback_hist[num_sources][1];
  TH2F *singles_180[num_sources][1];
  TH2F *addback_180[num_sources][1];
  TH2F *singadd_180[num_sources][1];

  int runstart[num_sources];
  int runstop[num_sources];
  for (int k = 0; k < num_sources; k++) {
    load_histograms(rootfile[k].c_str(), empty_gains_cal, singles_hist[k], addback_hist[k], singles_180[k], addback_180[k], singadd_180[k], k, runstart[k], runstop[k]);
    list->Add(singles_hist[k][0]);
    list->Add(addback_hist[k][0]);
    list->Add(singles_180[k][0]);
    list->Add(addback_180[k][0]);
    list->Add(singadd_180[k][0]);

  }
  //writes out the histograms in the TList
  TFile *out = new TFile("EffHist.root", "RECREATE");
  out->cd();
  list->Write();
  out->Close();//closes files

  // Writes information to text file for next code
  // File contains source information and run lengths
  ofstream outfile("histinfo.dat");
  for(int i = 0; i < num_sources; i++ ) outfile << sourceID[i] << "\t" << runstart[i] << "\t" << runstop[i] << endl;
  outfile.close();
}//histogram_make
