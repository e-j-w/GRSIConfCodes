Double_t num_bins = 8192; //number of bins to read in from histograms
Double_t min_bin = 0; //first bin in histogram
Double_t max_bin = 4096; //last bin in histogram
const int num_cores = 64; //number of cores in TIGRESS detector

//sorts an array into ascending order
//preserves correspondance between the main array and the 'carry' array
void insertion_sort(Double_t array[], Double_t carry[], int len) {
  int i;
  for (i = 1; i < len; i++) {
    Double_t val = array[i];
    Double_t carry_val = carry[i];
    for (int j = i - 1; j >= 0; j--) {
      if (array[j+1] < array[j]) {
        array[j+1] = array[j];
        carry[j+1] = carry[j];
        array[j] = val;
        carry[j] = carry_val;
      } else {
        break;
      }//else
    }//for
  }//for
}//insertion_sort

//Loads histogram from one subrun only
void load_histograms(const char analysis_filepath[], char cal_filepath[], TH1F *hist[], Int_t source_count, Int_t num_cores) {

  //opening analysis root file
  TFile *input_file = new TFile(analysis_filepath, "READ");
  if (!input_file->IsOpen()) {
    cout << "Cannot open input file!" << endl;
    return 0;
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
      cout << "ERROR: No TTigress or TGriffin Branch Found." << endl;
      exit(-1);
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

//using gamma ray spectra and actual energy peaks, calculates a linear equation to calibrate the detector's outputs
void get_centroids(TH1F *hist[], Double_t centroids[num_cores], Double_t centroid_errors[num_cores]) {

  //arrays to temporarily store the centroids of the peaks
  Double_t c[num_cores][2];
  Double_t c_er[num_cores][2];

  for (int i = 0; i < num_cores; i++) {
    int max_bin = hist[i]->GetXaxis()->GetXmax();
    int min_bin = hist[i]->GetXaxis()->GetXmin();
    hist[i]->SetAxisRange(100, max_bin-100,"X");
    hist[i]->SetBinContent(1, 0);
    Double_t intr = hist[i]->Integral(min_bin, max_bin, "width");
    if (intr < 1000) {
      cout << i << " FAILED!" << endl;
      for (int j = 0; j < 2; j++) {
        c[i][j] = -1; 
      }//for
      continue;
    }//if

    //roughly locates the peaks in the spectrum
    TSpectrum *spec = new TSpectrum(2 * 2);
    Int_t num_found = spec->Search(hist[i], 2, "", 0.5);

    cout << "Found " << num_found << " peaks in histogram." << endl;

    //if too many or too few peaks have been found,
    //something has gone wrong and we move on to the next core
    if (num_found != 2) {
      for (int j = 0; j < 2; j++) {
        c[i][j] = -1;
      }//for
	    continue;
    }//if

    //we fit a gaussian distribution around each peak
    //using the TSpectrum's data as a starting point
    TF1 *fit[2];
    Double_t* x_pos_array = spec->GetPositionX();

    for (int j = 0; j < num_found; j++) {
      Double_t x_pos = x_pos_array[j];
      Int_t bin = hist[i]->GetXaxis()->FindBin(x_pos);
      Double_t y_pos = hist[i]->GetBinContent(bin);

      fit[j] = new TF1(Form("Fit %i-%i", i, j), "[0]*(exp(-((x-[1])^2/(2*[2]^2))))*(1+ROOT::Math::erf([5]*((x-[1]))/([2]*pow(2,0.5)))) + [3] + [4]*(0.5*(1-(ROOT::Math::erf(((x-[1])/([2]*2^(0.5)))))))", x_pos - 50, x_pos + 50);
      fit[j]->SetParameters(y_pos, x_pos, 1, 15, 1, -1);
      fit[j]->SetParLimits(0, 10, 1e6); //area
      fit[j]->SetParLimits(1, x_pos - 10, x_pos + 10); //centroid
      fit[j]->SetParLimits(2, 0.2, 15); //sigma
      fit[j]->SetParLimits(4, 0.1, 100); //magnitude of step in background noise
      fit[j]->SetParLimits(5, -10, -0.1); //background noise constant

      //fitting the equation and storing the calculated centroids
      hist[i]->Fit(fit[j], "RQ+");
      c[i][j] = fit[j]->GetParameter(1);
      c_er[i][j] = fit[j]->GetParError(1);
      cout << fit[j]->GetParameter(2) << endl;
      cout << "CENTROID PROCESSED: Graph " << i << " Guess: " << x_pos << " Actual: " << c[i][j] << endl;
    }//for

    //to make sure the centroids match up to the correct energies
    //we sort the centroids in ascending order
    //making sure the centroids_er keep the correspondance
    insertion_sort(c[i], c_er[i], 2);
  }//for



  //we graph the centroids vs their corresponding energies
  //and fit a linear equation on the points
  for (int i = 0; i < num_cores; i++) {
    //if the histogram is empty, skip this core
    if (c[i][1] == -1) {
      centroids[i] = -1;
      centroid_errors[i] = -1;
      continue;
    }//if

    centroids[i] = c[i][1];
    centroid_errors[i] = c_er[i][1];

  } //for
}//get_centroids

//using gamma ray spectra and actual energy peaks, calculates a linear equation to calibrate the detector's outputs
void get_moments(const char analysis_filepath[], char cal_filepath[], const Double_t centroids[num_cores], Double_t means[num_cores], Double_t variances[num_cores], Double_t skewnesses[num_cores]) {

  //opening analysis root file
  TFile *input_file = new TFile(analysis_filepath, "READ");
  if (!input_file->IsOpen()) {
    cout << "Cannot open input file!" << endl;
    return 0;
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
    tig = true;
  } else {
    if (analysis->FindBranch("TGriffin")) {
      analysis->SetBranchAddress("TGriffin", & griffin);
      tig = false;
    } else {
      cout << "ERROR: No TTigress or TGriffin Branch Found." << endl;
      exit(-1);
    }
  }

  //initialize values
  uint64_t numCounts[num_cores];
  memset(numCounts,0,sizeof(numCounts));
  for(int i = 0; i < num_cores; i++){
    means[i] = 0.0;
    variances[i] = 0.0;
    skewnesses[i] = 0.0;
  }

  //get means from hit data
  for (int i = 0; i < num_entries - 1; i++) {
    tree->GetEntry(i);
    if(tig){
      for (int j = 0; j < tigress->GetMultiplicity(); j++) {
        TTigressHit *tigress_hit = tigress->GetTigressHit(j);
        if(tigress_hit->GetCharge() > 0.0 && tigress_hit->GetCharge() < 2000.0){
          if(fabs(tigress_hit->GetCharge() - centroids[tigress_hit->GetArrayNumber()]) < (centroids[tigress_hit->GetArrayNumber()]/50.0)){
            numCounts[tigress_hit->GetArrayNumber()]++;
            means[tigress_hit->GetArrayNumber()] += tigress_hit->GetCharge();
          }
        }
      }//for
    }else{
      for (int j = 0; j < griffin->GetMultiplicity(); j++) {
        TGriffinHit *griffin_hit = griffin->GetGriffinHit(j);
        if(griffin_hit->GetCharge() > 0.0 && griffin_hit->GetCharge() < 2000.0){
          if(fabs(griffin_hit->GetCharge() - centroids[griffin_hit->GetArrayNumber()-1]) < (centroids[griffin_hit->GetArrayNumber()-1]/50.0)){
            numCounts[griffin_hit->GetArrayNumber()-1]++;
            means[griffin_hit->GetArrayNumber()-1] += griffin_hit->GetCharge();
          }
        }
      }//for
    }
    if (i % 10000 == 0) {
      cout << setiosflags(ios::fixed) << "Mean: entry " << i << " of " << num_entries << ", " << 100 * (i) / num_entries << "% complete" << "\r" << flush;
    }//if
  }//for
  for(int i = 0; i < num_cores; i++){
    means[i] /= 1.0*numCounts[i];
  }

  //get variances from hit data
  for (int i = 0; i < num_entries - 1; i++) {
    tree->GetEntry(i);
    if(tig){
      for (int j = 0; j < tigress->GetMultiplicity(); j++) {
        TTigressHit *tigress_hit = tigress->GetTigressHit(j);
        if(tigress_hit->GetCharge() > 0.0 && tigress_hit->GetCharge() < 2000.0){
          if(fabs(tigress_hit->GetCharge() - centroids[tigress_hit->GetArrayNumber()]) < (centroids[tigress_hit->GetArrayNumber()]/50.0)){
            variances[tigress_hit->GetArrayNumber()] += pow(tigress_hit->GetCharge() - means[tigress_hit->GetArrayNumber()],2);
          }
        }
      }//for
    }else{
      for (int j = 0; j < griffin->GetMultiplicity(); j++) {
        TGriffinHit *griffin_hit = griffin->GetGriffinHit(j);
        if(griffin_hit->GetCharge() > 0.0 && griffin_hit->GetCharge() < 2000.0){
          if(fabs(griffin_hit->GetCharge() - centroids[griffin_hit->GetArrayNumber()-1]) < (centroids[griffin_hit->GetArrayNumber()-1]/50.0)){
            variances[griffin_hit->GetArrayNumber()-1] += pow(griffin_hit->GetCharge() - means[griffin_hit->GetArrayNumber()-1],2);
          }
        }
      }//for
    }
    if (i % 10000 == 0) {
      cout << setiosflags(ios::fixed) << "Variance: entry " << i << " of " << num_entries << ", " << 100 * (i) / num_entries << "% complete" << "\r" << flush;
    }//if
  }//for
  for(int i = 0; i < num_cores; i++){
    variances[i] /= 1.0*numCounts[i];
  }

  //get skewnesses from hit data
  for (int i = 0; i < num_entries - 1; i++) {
    tree->GetEntry(i);
    if(tig){
      for (int j = 0; j < tigress->GetMultiplicity(); j++) {
        TTigressHit *tigress_hit = tigress->GetTigressHit(j);
        if(tigress_hit->GetCharge() > 0.0 && tigress_hit->GetCharge() < 2000.0){
          if(fabs(tigress_hit->GetCharge() - centroids[tigress_hit->GetArrayNumber()]) < (centroids[tigress_hit->GetArrayNumber()]/50.0)){
            skewnesses[tigress_hit->GetArrayNumber()] += pow((tigress_hit->GetCharge() - means[tigress_hit->GetArrayNumber()])/sqrt(variances[tigress_hit->GetArrayNumber()]),3);
          }
        }
      }//for
    }else{
      for (int j = 0; j < griffin->GetMultiplicity(); j++) {
        TGriffinHit *griffin_hit = griffin->GetGriffinHit(j);
        if(griffin_hit->GetCharge() > 0.0 && griffin_hit->GetCharge() < 2000.0){
          if(fabs(griffin_hit->GetCharge() - centroids[griffin_hit->GetArrayNumber()-1]) < (centroids[griffin_hit->GetArrayNumber()-1]/50.0)){
            skewnesses[griffin_hit->GetArrayNumber()-1] += pow((griffin_hit->GetCharge() - means[griffin_hit->GetArrayNumber()-1])/sqrt(variances[griffin_hit->GetArrayNumber()-1]),3);
          }
        }
      }//for
    }
    if (i % 10000 == 0) {
      cout << setiosflags(ios::fixed) << "Skewness: entry " << i << " of " << num_entries << ", " << 100 * (i) / num_entries << "% complete" << "\r" << flush;
    }//if
  }//for
  for(int i = 0; i < num_cores; i++){
    skewnesses[i] /= 1.0*numCounts[i];
  }

}//get_moments

//Main method to be executed by GRSISort
void co60_moments(const char *analysis_file, const char *cal_file) {

  Double_t co60_ener[2] = {1173.240, 1332.508};
  Double_t co60_ener_e[2] = {0.003, 0.004};
  char afile[64]; strncpy(afile,analysis_file,64);
  char empty_gains_cal[64]; strncpy(empty_gains_cal,cal_file,64);

  TH1F *lin_hist[num_cores]; //histograms for each core

  //linear calibration equation parameters for each core
  Double_t centroids[num_cores];
  Double_t centroid_errors[num_cores];
  Double_t means[num_cores];
  Double_t variances[num_cores];
  Double_t skewnesses[num_cores];
  Double_t mean_sum = 0;
  Double_t variance_sum = 0;
  Double_t skewness_sum = 0;

  //loading in histograms from analysis root file
  load_histograms(afile, empty_gains_cal, lin_hist, 1, num_cores);

  //linear fit to get the centroids of the 1332 keV peak
  get_centroids(lin_hist, centroids, centroid_errors);

  //compute moments from original data in the region around the 1332 keV peak centroid
  get_moments(afile, empty_gains_cal, centroids, means, variances, skewnesses);
  cout << endl;
  
  //write out linear equations parameters
  ofstream coeff_fw("moments.txt");

  for (int i = 0; i < num_cores; i++) {
    coeff_fw << means[i] << " ";
    coeff_fw << variances[i] << " ";
    coeff_fw << skewnesses[i] << endl;
    mean_sum += means[i];
    variance_sum += variances[i];
    skewness_sum += skewnesses[i];
    cout << "Position " << i << "\t mean: " << means[i] << "\t variance: " << variances[i] << "\t skewness: " << skewnesses[i] << "\t Resolution: " << sqrt(variances[i])*2.355*co60_ener[1]/means[i] << " keV" << endl;
  }//for
  cout << "-----------" << endl;
  cout << "Average values: " << mean_sum/num_cores << "\t " << variance_sum/num_cores << "\t " << skewness_sum/num_cores << endl;

  coeff_fw.close();

}//linear_energy_RF
