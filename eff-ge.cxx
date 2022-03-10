// g++ eff-ge.cxx -std=c++0x  -o Efficiency
#include <iostream>
#include <iomanip>
#include <string.h>
#include <stdio.h>
#include <fstream>

using namespace std;
void efficiency(string sourceid[], string rootfiles[], int NumSources ) {
  string histname = "EffHist.root";
  string temp;
  bool goodroot = true;

  for(int i = 0; i < NumSources; i++ ) {
    ifstream f(rootfiles[i].c_str());
    if(!f.good()) {
      cout << "File: " << rootfiles[i] << " Does Not Exists" << endl;
      goodroot = false;
    }
  }
  if(goodroot){
    ofstream outfile("runlist.dat");
    for(int i = 0; i < NumSources; i++ ) outfile << sourceid[i] << "\t" << rootfiles[i] << endl;
    outfile.close();
    ifstream h(histname.c_str());
    if(h.good()){
      cout << "Efficiency histograms have already been sorted" << endl;
      cout << "Resort y/n?" << endl;
      cin >> temp;
      if(strcmp(temp.c_str(),"y") == 0 || strcmp(temp.c_str(),"Y") == 0) {
        cout << "Resorting efficiency histograms" << endl;
        system("grsisort -lq histogram_make.c");
      } else {
        cout << "Using existing histograms" << endl;
      }
    } else {
      system("grsisort -lq histogram_make.c");
    }
    system("grsisort -lq fit_efficiency.c");
  }
}
int main(int argc, char **argv){

  string sourceid[10];
  string rootfiles[10];
  int numsources = ((argc-1)/2);
  int check = ((argc-1)%2);
  int count = 0;
  if(argc > 2 && check==0) {
    for(int i = 1; i < argc; i++) {
      if(i%2 == 1){
        sourceid[count] = argv[i];
      } else if(i%2 == 0){
	rootfiles[count] = argv[i];
	count++;
      }
    }
    efficiency(sourceid, rootfiles, numsources);
  } else {
    cout << "Wrong Number of Arguments, need to be in form;" << "\n" << "SourceID1 AnalysisTree1 ... SourceIDX AnalysisTreeX" << endl;
  }
  return 0;
}
