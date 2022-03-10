#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <iostream>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <dirent.h>
#include <algorithm>
//#include <variables.h>
//Source_names
string source_names[] = {
"eu152", "ba133", "co56", "co60", "bi207"
};

//number of peaks used for each source's calibration
Int_t num_peaks[] = {
11, 7, 15, 2, 5
};

//number of sum pairs
Int_t num_sum[] = {
2, 2, 1, 1, 1
};

double source_halfl[] = {4933.705, 3851.115, 77.236, 1925.28, 11515.75};
double source_halfl_e[] = {3.285, 4.015, 0.026, 0.14, 14.6};

// Use summing corrections if true
bool summingcorrections = false;
float fC_sing = 1.0;
float fC_add = 1.0;
// Compare Source IDs
int search_array(vector<string> array, string search, int len) {
  for (int i = 0; i < len; i++) {
    if (search == array.at(i)) {
      return i;
    }//if
  }//for
  return -1;
}//search_array

// Finds Sum in peak Energy from gate + fit
int search_array(Double_t array[], Double_t gate, Double_t fit, int len) {
  for (int i = 0; i < len; i++) {
    if (abs((gate + fit) - array[i]) < 2.0 ) {
      return i;
    }//if
  }//for
  return -1;
}//search_array

// Functions to calculate Areas + errors of Peaks
double GetArea(TH1F* h, TF1 *f, TF1 *bg){

  int firstbin = h->FindBin(f->GetParameter(1) - 4 * f->GetParameter(2));
  int lastbin = h->FindBin(f->GetParameter(1) + 4 * f->GetParameter(2));

  double firstedge = h->GetBinLowEdge(firstbin);
  double lastedge = h->GetBinLowEdge(lastbin) + h->GetBinWidth(lastbin);

  double bg_sub_integral = f->Integral(firstedge,lastedge)/h->GetBinWidth(1) - bg->Integral(firstedge,lastedge)/h->GetBinWidth(1);

  return bg_sub_integral;
}

double GetAreaError(TH1F* h, TF1 *f, TF1 *bg){

  int firstbin = h->FindBin(f->GetParameter(1) - 4 * f->GetParameter(2));
  int lastbin = h->FindBin(f->GetParameter(1) + 4 * f->GetParameter(2));

  double firstedge = h->GetBinLowEdge(firstbin);
  double lastedge = h->GetBinLowEdge(lastbin) + h->GetBinWidth(lastbin);

  double bg_sub_err = TMath::Sqrt( f->Integral(firstedge,lastedge)/h->GetBinWidth(1) + bg->Integral(firstedge,lastedge)/h->GetBinWidth(1));

  return bg_sub_err;
}
double GetArea(TH1D* h, TF1 *f, TF1 *bg){

  int firstbin = h->FindBin(f->GetParameter(1) - 4 * f->GetParameter(2));
  int lastbin = h->FindBin(f->GetParameter(1) + 4 * f->GetParameter(2));

  double firstedge = h->GetBinLowEdge(firstbin);
  double lastedge = h->GetBinLowEdge(lastbin) + h->GetBinWidth(lastbin);

  double bg_sub_integral = f->Integral(firstedge,lastedge)/h->GetBinWidth(1) - bg->Integral(firstedge,lastedge)/h->GetBinWidth(1);
  double bg_sub_err = TMath::Sqrt( f->Integral(firstedge,lastedge)/h->GetBinWidth(1) + bg->Integral(firstedge,lastedge)/h->GetBinWidth(1));

  return bg_sub_integral;
}

double GetAreaError(TH1D* h, TF1 *f, TF1 *bg){

  int firstbin = h->FindBin(f->GetParameter(1) - 4 * f->GetParameter(2));
  int lastbin = h->FindBin(f->GetParameter(1) + 4 * f->GetParameter(2));

  double firstedge = h->GetBinLowEdge(firstbin);
  double lastedge = h->GetBinLowEdge(lastbin) + h->GetBinWidth(lastbin);

  double bg_sub_integral = f->Integral(firstedge,lastedge)/h->GetBinWidth(1) - bg->Integral(firstedge,lastedge)/h->GetBinWidth(1);
  double bg_sub_err = TMath::Sqrt( f->Integral(firstedge,lastedge)/h->GetBinWidth(1) + bg->Integral(firstedge,lastedge)/h->GetBinWidth(1));

  return bg_sub_err;
}

// Loads histograms from file created by histogram_make.C and makes projection of coincidence matrix fot Sum out corrections
void load_histograms(const char histogram_filepath[], TH1F *singles[], TH1F *addback[], TH2F *singles_sum[], TH2F *addback_sum[], TH1D *singles_sum_pro[], TH1D *addback_sum_pro[], Int_t source_count) {
  TFile *infile = new TFile(histogram_filepath);
  singles[0] = (TH1F *)infile->Get(Form("singles_%i",source_count));
  addback[0] = (TH1F *)infile->Get(Form("addback_%i",source_count));
  singles_sum[0] = (TH2F *)infile->Get(Form("singles180_%i",source_count));
  addback_sum[0] = (TH2F *)infile->Get(Form("addback180_%i",source_count));
  singles_sum[0]->Scale(1.0/fC_sing);
  addback_sum[0]->Scale(1.0/fC_add);

  char hname[20];
  sprintf(hname,"singles_sum_pro_%i",source_count);
  singles_sum_pro[0] = new TH1D(hname,Form("singles_sum_pro_%i",source_count),singles_sum[0]->GetNbinsX(),singles_sum[0]->GetXaxis()->GetXmin(),singles_sum[0]->GetXaxis()->GetXmax());
  singles_sum[0]->ProjectionY(hname,0,-1);

  char aname[20];
  sprintf(aname,"addback_sum_pro_%i",source_count);
  addback_sum_pro[0] = new TH1D(aname,"",addback_sum[0]->GetNbinsX(),addback_sum[0]->GetXaxis()->GetXmin(),addback_sum[0]->GetXaxis()->GetXmax());
  addback_sum[0]->ProjectionY(aname,0,-1);
}

// Fits the spectra and gets peak areas
void fit_peaks(TH1F *hist[], Double_t energy[], Double_t energy_er[], Int_t num_peaks_used, Double_t area[], Double_t area_error[], Int_t fit_min_range[], Int_t fit_max_range[]) {
  for (int j = 0; j < num_peaks_used; j++) {
    if(area[j] != 0) continue; // Skips if already fit as a doublet
    if(signbit(energy[j]) == 0) { // If negative energy given fit peak as a doublet
      Int_t bin_guess = hist[0]->GetXaxis()->FindBin(energy[j]);
      Double_t y_guess = hist[0]->GetBinContent(bin_guess);
      if(fit_min_range[j] == 0)fit_min_range[j] = 10;
      if(fit_max_range[j] == 0)fit_max_range[j] = 10;
      TF1 *fit = new TF1(Form("fit %i", j),
	"[0]*((((exp(-((x-[1])^2/(2*[2]^2)))))*[4])+(exp((x-[1])/([3]))*(ROOT::Math::erfc(((x-[1])/([2]*2^(0.5)))+([2]/([3]*2^(0.5)))))*(1-[4]))) + [5] + [6]*(0.5*(1-(ROOT::Math::erf(((x-[1])/([2]*2^(0.5)))))))",
         energy[j] - fit_min_range[j], energy[j] + fit_max_range[j]);
      fit->SetParameters(y_guess, energy[j], 2, 5, 0.5,100,-5);
      fit->SetParLimits(0, y_guess*0.1, y_guess*5); //area
      fit->SetParLimits(1, energy[j] - 1, energy[j] + 1); //centroid
      fit->SetParLimits(2, 0.2, 5); // sigma of gaussian distribution
      fit->SetParLimits(3, 0.1, 100); // Controls low energy tail
      fit->SetParLimits(4, 0.01, 1); // Controls asymmetry 

      //fitting equation and saving centroids
      TFitResultPtr r = hist[0]->Fit(fit,"QR+");
      Int_t fitStatus = r;
      // Refits if fit fails
      if(r!=0) {
        cout << "\033[1;31mFit Failed on " << energy[j] << " Refitting \033[0m" << endl;
        fit->SetParLimits(0, fit->GetParameter(0)*0.5, fit->GetParameter(0)*5.0); //area
	fit->SetParameter(0, fit->GetParameter(0)* (1.5 - (rand() % 10)/10.0));
        fit->SetParameter(3, fit->GetParameter(3)* (1.5 - (rand() % 10)/10.0));
        fit->SetParameter(4, fit->GetParameter(4)* (1.5 - (rand() % 10)/10.0));
        TFitResultPtr r = hist[0]->Fit(fit,"QR+");
        if(r!=0) cout << "\033[1;31mFit Failed on " << energy[j] << "\033[0m" << endl;
      }

      TF1 *bg = new TF1(Form("bg_fit %i", j)," [0] + [1]*(0.5*(1-(ROOT::Math::erf(((x-[2])/([3]*2^(0.5)))))))", energy[j] - fit_min_range[j], energy[j] + fit_max_range[j]);
      bg->SetParameters(fit->GetParameter(5),fit->GetParameter(6),fit->GetParameter(1),fit->GetParameter(2));
      Double_t temp = GetArea(hist[0],fit,bg);
      area[j] = temp;
      temp = GetAreaError(hist[0],fit,bg);
      area_error[j] = temp;
    } else { // Double peak fit
      Int_t bin_guess = hist[0]->GetXaxis()->FindBin(abs(energy[j]));
      Int_t bin_guess2 = hist[0]->GetXaxis()->FindBin(abs(energy[j+1]));
      Double_t y_guess = hist[0]->GetBinContent(bin_guess);
      Double_t y_guess2 = hist[0]->GetBinContent(bin_guess2);
      if(fit_min_range[j] == 0)fit_min_range[j] = 10;
      if(fit_max_range[j] == 0)fit_max_range[j+1] = 10;
      TF1 *fit = new TF1(Form("Double fit %i %i", j, j + 1),
        "[0]*((((exp(-((x-[1])^2/(2*[2]^2)))))*[4])+(exp((x-[1])/([3]))*(ROOT::Math::erfc(((x-[1])/([2]*2^(0.5)))+([2]/([3]*2^(0.5)))))*(1-[4]))) + [5] + [6]*[0]*(0.5*(1-(ROOT::Math::erf(((x-[1])/([2]*2^(0.5))))))) + [7]*((((exp(-((x-[8])^2/(2*[2]^2)))))*[4])+(exp((x-[1])/([3]))*(ROOT::Math::erfc(((x-[1])/([2]*2^(0.5)))+([2]/([3]*2^(0.5)))))*(1-[4]))) + [7]*[6]*(0.5*(1-(ROOT::Math::erf(((x-[8])/([2]*2^(0.5)))))))",
         abs(energy[j]) - fit_min_range[j], abs(energy[j+1]) + fit_max_range[j+1]);
      fit->SetParameters(y_guess, abs(energy[j]), 2, 5, 0.5, 100, -0.05, y_guess, abs(energy[j+1]));
      fit->SetParLimits(0, y_guess*0.01, y_guess*5); //area
      fit->SetParLimits(1, abs(energy[j]) - 1, abs(energy[j]) + 1); //centroid
      fit->SetParLimits(2, 0.2, 5); //sigma of gaussian distribution
      fit->SetParLimits(3, 0.1, 100); //magnitude of step in background noise
      fit->SetParLimits(4, 0.01, 1); //magnitude of step in background noise
      fit->SetParLimits(7, y_guess2*0.01, y_guess2*5); //area
      fit->SetParLimits(8, abs(energy[j+1]) - 1, abs(energy[j+1]) + 1); //centroid

      TFitResultPtr r = hist[0]->Fit(fit,"QR+");
      Int_t fitStatus = r;
      // Refits if fit failed
      if(r!=0) {
        cout << "\033[1;31mFit Failed on Doublet " << abs(energy[j]) << " " <<  abs(energy[j+1]) << " Refitting \033[0m" << endl;
        TFitResultPtr r = hist[0]->Fit(fit,"QR+");
        if(r!=0) cout << "\033[1;31mFit Failed on " << energy[j] << "\033[0m" << endl;
      }

      TF1 *bg = new TF1(Form("bg_fit %i %i", j, j+1),"[0] + [1]*(0.5*(1-(ROOT::Math::erf(((x-[2])/([3]*2^(0.5))))))) + [4]*(0.5*(1-(ROOT::Math::erf(((x-[5])/([3]*2^(0.5)))))))", abs(energy[j]) - fit_min_range[j], abs(energy[j+1]));
      bg->SetParameters(fit->GetParameter(5), fit->GetParameter(6)*fit->GetParameter(0), fit->GetParameter(1), fit->GetParameter(2), fit->GetParameter(7)*fit->GetParameter(6), fit->GetParameter(8));

      TF1 *p1 = new TF1(Form("Double peak %i", j),"[0]*((((exp(-((x-[1])^2/(2*[2]^2)))))*[4])+(exp((x-[1])/([3]))*(ROOT::Math::erfc(((x-[1])/([2]*2^(0.5)))+([2]/([3]*2^(0.5)))))*(1-[4]))) + [5] + [6]*(0.5*(1-(ROOT::Math::erf(((x-[1])/([2]*2^(0.5))))))) + [7]*(0.5*(1-(ROOT::Math::erf(((x-[8])/([2]*2^(0.5)))))))",
	abs(energy[j]) - fit_min_range[j], abs(energy[j+1]) + fit_max_range[j+1]);
      p1->SetParameters(fit->GetParameter(0), fit->GetParameter(1), fit->GetParameter(2), fit->GetParameter(3), fit->GetParameter(4), fit->GetParameter(5), fit->GetParameter(6)*fit->GetParameter(0), fit->GetParameter(7)*fit->GetParameter(6), fit->GetParameter(8));


      TF1 *p2 = new TF1(Form("Double peak %i", j+1),"[0]*((((exp(-((x-[1])^2/(2*[2]^2)))))*[4])+(exp((x-[1])/([3]))*(ROOT::Math::erfc(((x-[1])/([2]*2^(0.5)))+([2]/([3]*2^(0.5)))))*(1-[4]))) + [5] + [6]*(0.5*(1-(ROOT::Math::erf(((x-[1])/([2]*2^(0.5))))))) + [7]*(0.5*(1-(ROOT::Math::erf(((x-[8])/([2]*2^(0.5)))))))",
         abs(energy[j]) - fit_min_range[j], abs(energy[j+1]) + fit_max_range[j+1]);
      p2->SetParameters(fit->GetParameter(7), fit->GetParameter(8), fit->GetParameter(2), fit->GetParameter(3), fit->GetParameter(4), fit->GetParameter(5), fit->GetParameter(7)*fit->GetParameter(6), fit->GetParameter(0)*fit->GetParameter(6), fit->GetParameter(1));

      Double_t temp = GetArea(hist[0],p1,bg);
      area[j] = temp;
      temp = GetAreaError(hist[0],p1,bg);
      area_error[j] = temp;

      temp = GetArea(hist[0],p2,bg);
      area[j+1] = temp;
      temp = GetAreaError(hist[0],p2,bg);
      area_error[j+1] = temp;
    }
   }//for
}//find_centroids

// Does the sumout corrections
void sumout_cor(TH1D *hist[], Double_t energy[], Int_t num_peaks_used, Double_t area[], Double_t area_error[], Int_t fit_min_range[], Int_t fit_max_range[]) {
  Double_t temp_area[num_peaks_used];
  Double_t temp_errors[num_peaks_used];
  for(int i = 0; i < num_peaks_used; i++) temp_area[i] = 0;
  for (int j = 0; j < num_peaks_used; j++) {
    if(temp_area[j] > 1) continue; // Skips if fit as doublet
    if(signbit(energy[j]) == 0) {
      Int_t bin_guess = hist[0]->GetXaxis()->FindBin(energy[j]);
      Double_t y_guess = hist[0]->GetBinContent(bin_guess);
      if(fit_min_range[j] == 0)fit_min_range[j] = 10;
      if(fit_max_range[j] == 0)fit_max_range[j] = 10;
      TF1 *fit = new TF1(Form("fit %i", j),
	"[0]*((((exp(-((x-[1])^2/(2*[2]^2)))))*[4])+(exp((x-[1])/([3]))*(ROOT::Math::erfc(((x-[1])/([2]*2^(0.5)))+([2]/([3]*2^(0.5)))))*(1-[4]))) + [5] + [6]*(0.5*(1-(ROOT::Math::erf(((x-[1])/([2]*2^(0.5)))))))",
         energy[j] - fit_min_range[j], energy[j] + fit_max_range[j]);
      fit->SetParameters(y_guess, energy[j], 2, 5, 0.5,100,-5);
      fit->SetParLimits(0, y_guess*0.1, y_guess*5); //area
      fit->SetParLimits(1, energy[j] - 1, energy[j] + 1); //centroid
      fit->SetParLimits(2, 0.2, 5); //sigma of gaussian distribution
      fit->SetParLimits(3, 0.1, 100);
      fit->SetParLimits(4, 0.01, 1);

      //fitting equation and saving centroids
      hist[0]->Fit(fit, "QR+");
      TFitResultPtr r = hist[0]->Fit(fit,"QR+");
      Int_t fitStatus = r;
      if(r!=0) {
        cout << "\033[1;31mFit Failed on Sumout " << energy[j] << " Refitting \033[0m" << endl;
        fit->SetParLimits(0, fit->GetParameter(0)*0.5, fit->GetParameter(0)*5.0); //area
	fit->SetParameter(0, fit->GetParameter(0)* (1.5 - (rand() % 10)/10.0));
        fit->SetParameter(3, fit->GetParameter(3)* (1.5 - (rand() % 10)/10.0));
        TFitResultPtr r = hist[0]->Fit(fit,"QR+");
        if(r!=0) cout << "\033[1;31mFit Failed on " << energy[j] << "\033[0m" << endl;
      }

      TF1 *bg = new TF1(Form("bg_fit %i", j)," [0] + [1]*(0.5*(1-(ROOT::Math::erf(((x-[2])/([3]*2^(0.5)))))))", energy[j] - fit_min_range[j], energy[j] + fit_max_range[j]);
      bg->SetParameters(fit->GetParameter(5),fit->GetParameter(6),fit->GetParameter(1),fit->GetParameter(2));

      temp_area[j] = GetArea(hist[0],fit,bg);
      Double_t temp = area[j];
      area[j] = temp + temp_area[j];

      temp_errors[j] = GetAreaError(hist[0],fit,bg);
      temp = area_error[j];
      area_error[j] = TMath::Sqrt(pow(temp,2)+pow(temp_errors[j],2));

    } else { // Double peak fit
      Int_t bin_guess = hist[0]->GetXaxis()->FindBin(abs(energy[j]));
      Int_t bin_guess2 = hist[0]->GetXaxis()->FindBin(abs(energy[j+1]));
      Double_t y_guess = hist[0]->GetBinContent(bin_guess);
      Double_t y_guess2 = hist[0]->GetBinContent(bin_guess2);
      if(fit_min_range[j] == 0)fit_min_range[j] = 10;
      if(fit_max_range[j] == 0)fit_max_range[j+1] = 10;
      TF1 *fit = new TF1(Form("Double fit %i %i", j, j + 1),
        "[0]*((((exp(-((x-[1])^2/(2*[2]^2)))))*[4])+(exp((x-[1])/([3]))*(ROOT::Math::erfc(((x-[1])/([2]*2^(0.5)))+([2]/([3]*2^(0.5)))))*(1-[4]))) + [5] + [6]*[0]*(0.5*(1-(ROOT::Math::erf(((x-[1])/([2]*2^(0.5))))))) + [7]*((((exp(-((x-[8])^2/(2*[2]^2)))))*[4])+(exp((x-[1])/([3]))*(ROOT::Math::erfc(((x-[1])/([2]*2^(0.5)))+([2]/([3]*2^(0.5)))))*(1-[4]))) + [7]*[6]*(0.5*(1-(ROOT::Math::erf(((x-[8])/([2]*2^(0.5)))))))",
         abs(energy[j]) - fit_min_range[j], abs(energy[j+1]) + fit_max_range[j+1]);
      fit->SetParameters(y_guess, abs(energy[j]), 2, 5, 0.5, 100, -0.05, y_guess, abs(energy[j+1]));
      fit->SetParLimits(0, y_guess*0.01, y_guess*5); //area
      fit->SetParLimits(1, abs(energy[j]) - 1, abs(energy[j]) + 1); //centroid
      fit->SetParLimits(2, 0.2, 5); //sigma of gaussian distribution
      fit->SetParLimits(3, 0.1, 100);
      fit->SetParLimits(4, 0.01, 1);
      fit->SetParLimits(7, y_guess2*0.01, y_guess2*5); //area
      fit->SetParLimits(8, abs(energy[j+1]) - 1, abs(energy[j+1]) + 1); //centroid

      TFitResultPtr r = hist[0]->Fit(fit,"QR+");
      Int_t fitStatus = r;
      if(r!=0) {
        cout << "\033[1;31mFit Failed on " << energy[j] << " Refitting \033[0m" << endl;
        TFitResultPtr r = hist[0]->Fit(fit,"QR+");
        if(r!=0) cout << "\033[1;31mFit Failed on " << energy[j] << "\033[0m" << endl;
      }

      TF1 *bg = new TF1(Form("bg_fit %i %i", j, j+1),"[0] + [1]*(0.5*(1-(ROOT::Math::erf(((x-[2])/([3]*2^(0.5))))))) + [4]*(0.5*(1-(ROOT::Math::erf(((x-[5])/([3]*2^(0.5)))))))", abs(energy[j]) - fit_min_range[j], abs(energy[j+1]));
      bg->SetParameters(fit->GetParameter(5), fit->GetParameter(6)*fit->GetParameter(0), fit->GetParameter(1), fit->GetParameter(2), fit->GetParameter(7)*fit->GetParameter(6), fit->GetParameter(8));

      TF1 *p1 = new TF1(Form("Double peak %i", j),"[0]*((((exp(-((x-[1])^2/(2*[2]^2)))))*[4])+(exp((x-[1])/([3]))*(ROOT::Math::erfc(((x-[1])/([2]*2^(0.5)))+([2]/([3]*2^(0.5)))))*(1-[4]))) + [5] + [6]*(0.5*(1-(ROOT::Math::erf(((x-[1])/([2]*2^(0.5))))))) + [7]*(0.5*(1-(ROOT::Math::erf(((x-[8])/([2]*2^(0.5)))))))",
	abs(energy[j]) - fit_min_range[j], abs(energy[j+1]) + fit_max_range[j+1]);
      p1->SetParameters(fit->GetParameter(0), fit->GetParameter(1), fit->GetParameter(2), fit->GetParameter(3), fit->GetParameter(4), fit->GetParameter(5), fit->GetParameter(6)*fit->GetParameter(0), fit->GetParameter(7)*fit->GetParameter(6), fit->GetParameter(8));


      TF1 *p2 = new TF1(Form("Double peak %i", j+1),"[0]*((((exp(-((x-[1])^2/(2*[2]^2)))))*[4])+(exp((x-[1])/([3]))*(ROOT::Math::erfc(((x-[1])/([2]*2^(0.5)))+([2]/([3]*2^(0.5)))))*(1-[4]))) + [5] + [6]*(0.5*(1-(ROOT::Math::erf(((x-[1])/([2]*2^(0.5))))))) + [7]*(0.5*(1-(ROOT::Math::erf(((x-[8])/([2]*2^(0.5)))))))",
         abs(energy[j]) - fit_min_range[j], abs(energy[j+1]) + fit_max_range[j+1]);
      p2->SetParameters(fit->GetParameter(7), fit->GetParameter(8), fit->GetParameter(2), fit->GetParameter(3), fit->GetParameter(4), fit->GetParameter(5), fit->GetParameter(7)*fit->GetParameter(6), fit->GetParameter(0)*fit->GetParameter(6), fit->GetParameter(1));

      temp_area[j] = GetArea(hist[0],p1,bg);
      Double_t temp = area[j];
      area[j] = temp + temp_area[j];
      temp_errors[j] = GetAreaError(hist[0],p1,bg);
      temp = area_error[j];
      area_error[j] = TMath::Sqrt(pow(temp,2)+pow(temp_errors[j],2));

      temp_area[j+1] = GetArea(hist[0],p2,bg);
      temp = area[j+1];
      area[j+1] = temp + temp_area[j+1];
      temp_errors[j+1] = GetAreaError(hist[0],p2,bg);
      temp = area_error[j+1];
      area_error[j+1] = TMath::Sqrt(pow(temp,2)+pow(temp_errors[j+1],2));
    }
  }//for
}
// Does the sumin corrections
void sumin_cor(TH1D *hist_pro[], TH1D *hist_gate[], Double_t energy[], Int_t num_peaks_used, Double_t area[], Double_t area_error[], Int_t source_id, Int_t fit_min_range[], Int_t fit_max_range[], Double_t sum_gate[], Double_t sum_fit[], Int_t num_sum_peaks) {

    for(int i = 0; i < num_sum_peaks; i++ ){
      if(sum_fit[i] == 0 )continue;
      Int_t bin_guess = hist_gate[i]->GetXaxis()->FindBin(sum_fit[i]);
      Double_t y_guess = hist_gate[i]->GetBinContent(bin_guess);
      TF1 *fit_sub = new TF1("fit_sub",
        "[0]*((((exp(-((x-[1])^2/(2*[2]^2)))))*[4])+(exp((x-[1])/([3]))*(ROOT::Math::erfc(((x-[1])/([2]*2^(0.5)))+([2]/([3]*2^(0.5)))))*(1-[4]))) + [5]",
      sum_fit[i] - 8, sum_fit[i] + 8);
      fit_sub->SetParameters(y_guess, sum_fit[i], 2, 5, 0.5,100,-5);
      fit_sub->SetParLimits(0, y_guess*0.4, y_guess*2); //area
      fit_sub->SetParLimits(1, sum_fit[i] - 1, sum_fit[i] + 1); //centroid
      fit_sub->SetParLimits(2, 0.2, 5); //sigma of gaussian distribution
      fit_sub->SetParLimits(3, 0.1, 100);
      fit_sub->SetParLimits(4, 0.01, 1);
      //fitting equation and saving centroids
      TFitResultPtr r = hist_gate[i]->Fit(fit_sub,"QLR+");
      Int_t fitStatus = r;
      if(r!=0) {
        cout << "\033[1;31mFit Failed on Sumin " << sum_fit[i] << " Refitting \033[0m" << endl;
        TFitResultPtr r = hist_gate[i]->Fit(fit_sub,"QLR+");
        if(r!=0) cout << "\033[1;31mFit Failed on " << sum_fit[i] << "\033[0m" << endl;
      }

      TF1 *bg_sub = new TF1("bg_sub"," [0]", sum_fit[i] - fit_min_range[i], sum_fit[i] + fit_max_range[i]);
      bg_sub->SetParameter(0,fit_sub->GetParameter(5));
      cout << fit_sub->GetParameter(1) << "\t" << GetArea(hist_gate[i],fit_sub,bg_sub) << "\t" << GetAreaError(hist_gate[i],fit_sub,bg_sub) << endl;
      int j = search_array(energy, sum_gate[i], sum_fit[i], num_peaks_used);
      if(j == -1) cout << "No entry found for Combination Gate " << sum_gate[i] << "+ Peak Fit " << sum_fit[i] << " sum total = " << sum_gate[i] + sum_fit[i] << endl;
      area[j] -= GetArea(hist_gate[i],fit_sub,bg_sub);
      double temp = GetAreaError(hist_gate[i],fit_sub,bg_sub);
      double temp_error = area_error[j];
      area_error[j] = TMath::Sqrt(pow(temp,2)+pow(temp_error,2));
  }
}
// Converts source data to unix time
int date_to_unix(int year, int month, int day, int hour, int min, int sec){

  time_t rawtime;
  struct tm * timeinfo;
  time ( &rawtime ); //or: rawtime = time(0);
  timeinfo = localtime ( &rawtime );
  timeinfo->tm_year   = year - 1900;
  timeinfo->tm_mon    = month - 1;    //months since January - [0,11]
  timeinfo->tm_mday   = day;          //day of the month - [1,31]
  timeinfo->tm_hour   = hour;         //hours since midnight - [0,23]
  timeinfo->tm_min    = min;          //minutes after the hour - [0,59]
  timeinfo->tm_sec    = sec;          //seconds after the minute - [0,59]

  return mktime ( timeinfo );
}
// Calculated the source activity at measurement time
double RadActivity(int source_date, int run_info, double halflife, double radact) {

  // Calculating Time Difference (days)
  Double_t time_dif = (run_info - source_date)/(3600*24);

  // Calculating current radioactivity (Bq)
  return radact * pow(2, -1 * time_dif / halflife);
}
// Calculated the error in source activity 
double RadActivityError(int source_date, int run_info, double halflife, double halflife_er, double radact) {

  // Calculating Time Difference (days)
  Double_t time_dif = (run_info - source_date)/(3600*24);

  // Calculating current radioactivity (Bq)
  Double_t radact_cur = radact * pow(2, -1 * time_dif / halflife);

  return log(2) * radact_cur * halflife_er * time_dif / (halflife * halflife);
}

// Converts Peak areas to efficience to absolute efficiencies
void abseff(Int_t num_peaks_used, Double_t area[], Double_t area_error[], Double_t eff[], Double_t eff_error[], Double_t intensity[], Double_t intensity_error[], Int_t runtime, Double_t source_activity, Double_t source_activity_error) {
  for(int i = 0; i < num_peaks_used; i++) {
    eff[i] = area[i] / (intensity[i] * runtime * source_activity);
    eff_error[i] = sqrt(pow(area_error[i] / (intensity[i] * runtime * source_activity), 2) + pow(intensity_error[i] * area[i] / (pow(intensity[i], 2) * runtime * source_activity), 2) + pow(source_activity_error * area[i] / (intensity[i] * runtime * pow(source_activity, 2)), 2));
  }
return;
}//abseff

// Applies gates to make the spectra for the sum in corrections
void make_sum_hist(TH2F *hist[], TH1D *gate[], TH2F *addhist[], TH1D *addgate[], Double_t sum_gate[], Int_t source_count, Int_t num_sum_peaks){
  for(int i = 0; i < num_sum_peaks; i++) {
    char hname[20];
    sprintf(hname,"sum_gate_%i_%i",source_count,i);
    gate[i] = new TH1D(hname,"",hist[0]->GetNbinsX(),hist[0]->GetXaxis()->GetXmin(),hist[0]->GetXaxis()->GetXmax());
    hist[0]->ProjectionY(hname, hist[0]->GetXaxis()->FindBin(sum_gate[i] - 8), hist[0]->GetXaxis()->FindBin(sum_gate[i] + 8));

    char aname[20];
    sprintf(aname,"addback_sum_gate_%i_%i",source_count,i);
    addgate[i] = new TH1D(aname,"",addhist[0]->GetNbinsX(),addhist[0]->GetXaxis()->GetXmin(),addhist[0]->GetXaxis()->GetXmax());
    addhist[0]->ProjectionY(aname, addhist[0]->GetXaxis()->FindBin(sum_gate[i] - 8), addhist[0]->GetXaxis()->FindBin(sum_gate[i] + 8));
    TH1D *temp_left = new TH1D("temp_left","",addhist[0]->GetNbinsX(),addhist[0]->GetXaxis()->GetXmin(),addhist[0]->GetXaxis()->GetXmax());
    TH1D *temp_right = new TH1D("temp_right","",addhist[0]->GetNbinsX(),addhist[0]->GetXaxis()->GetXmin(),addhist[0]->GetXaxis()->GetXmax());
    addhist[0]->ProjectionY("temp_left",addhist[0]->GetXaxis()->FindBin(sum_gate[i] - 16),addhist[0]->GetXaxis()->FindBin(sum_gate[i] - 8));
    addhist[0]->ProjectionY("temp_right",addhist[0]->GetXaxis()->FindBin(sum_gate[i] + 8), addhist[0]->GetXaxis()->FindBin(sum_gate[i] + 16));
    addgate[i]->Add(temp_left,-0.5);
    addgate[i]->Add(temp_right,-0.5);
    temp_left->Delete();
    temp_right->Delete();
  }
}
void fit_efficiency() {

  /*******Initialization of Source Energy Peak Data********/
  Double_t** source_energy = new Double_t*[5];
  source_energy[0] = new Double_t[num_peaks[0]]{121.7817, 244.6974, 344.2785, 367.7891, 443.9606, 656.489, 778.9045, 964.057, 1112.076, 1212.948, 1408.018};
  source_energy[1] = new Double_t[num_peaks[1]]{80.9979, 160.612, 223.2368, 276.3989, 302.851, 356.013, 383.848}; //ba133
  source_energy[2] = new Double_t[num_peaks[2]]{846.770, 1037.843, 1175.101, 1238.2883, 1360.212, 1771.3567, 2015.215, 2034.791, 2598.500, 3009.645, 3202.029, -3253.416, -3273.079, 3451.232, 3548.05}; //co56
  source_energy[3] = new Double_t[num_peaks[3]]{1173.240, 1332.508}; //co60
  source_energy[4] = new Double_t[num_peaks[4]]{-72.805, -74.969, 569.698, 1063.656, 1770.228}; //bi207

  Double_t** source_energy_er = new Double_t*[5];
  source_energy_er[0] = new Double_t[num_peaks[0]]{0.0003, 0.0008, 0.0012, 0.0016, 0.0020, 0.005, 0.0024, 0.005, 0.003, 0.011, 0.003};
  source_energy_er[1] = new Double_t[num_peaks[1]]{0.0011, 0.0016, 0.0013, 0.0012, 0.0016, 0.0017, 0.0012}; //ba133
  source_energy_er[2] = new Double_t[num_peaks[2]]{0.002, 0.004, 0.004, 0.003, 0.004, 0.004, 0.005, 0.005, 0.004, 0.004, 0.009, 0.004, 0.004, 0.004}; //co56
  source_energy_er[3] = new Double_t[num_peaks[3]]{0.003, 0.004}; //co60
  source_energy_er[4] = new Double_t[num_peaks[4]]{0.00001, 0.00001, 0.002, 0.003, 0.009}; //bi207

  // For Sum in correction define the gate and peak combination to fit
  Double_t** source_sum_gate = new Double_t*[4];
  source_sum_gate[0] = new Double_t[num_sum[0]]{867.4, 1408.1}; //eu 152
  source_sum_gate[1] = new Double_t[num_sum[1]]{302.9, 276.4}; //ba133
  source_sum_gate[2] = new Double_t[num_sum[2]]{}; //co56
  source_sum_gate[3] = new Double_t[num_sum[3]]{}; //co60
  source_sum_gate[4] = new Double_t[num_sum[4]]{}; //bi207

  Double_t** source_sum_fit = new Double_t*[5];
  source_sum_fit[0] = new Double_t[num_sum[0]]{244.6974, 121.7817}; //eu 152
  source_sum_fit[1] = new Double_t[num_sum[1]]{80.9979,79.6}; //ba133
  source_sum_fit[2] = new Double_t[num_sum[2]]{}; //co56
  source_sum_fit[3] = new Double_t[num_sum[3]]{}; //co60
  source_sum_fit[4] = new Double_t[num_sum[4]]{}; //bi207

  // Absolute intensities of sources
  Double_t** source_abs = new Double_t*[5];
  source_abs[0] = new Double_t[num_peaks[0]]{0.2853, 0.0755, 0.2659, 0.00859, 0.03125, 0.001441, 0.1293,0.1465,0.13859,0.01415,0.2087};
  source_abs[1] = new Double_t[num_peaks[1]]{0.3555, 0.00638, 0.00453, 0.0716, 0.1834, 0.6205, 0.0894}; //ba133
  source_abs[2] = new Double_t[num_peaks[2]]{0.999399, 0.1405, 0.02252, 0.6646, 0.04283, 0.1541, 0.03016, 0.0777, 0.1697, 0.01036, 0.03209, 0.07923, 0.018759, 0.00949, 0.001955};  //co56
  source_abs[3] = new Double_t[num_peaks[3]]{0.9985, 0.999826};	//co60
  source_abs[4] = new Double_t[num_peaks[4]]{0.214, 0.357, 0.9775, 0.745, 0.0687}; //bi207

  Double_t** source_abs_e = new Double_t*[5];
  source_abs_e[0] = new Double_t[num_peaks[0]]{0.0016, 0.0004, 0.002, 0.0006, 0.00017, 0.000020, 0.0008,0.0013,0.00014,0.00008,0.00009};
  source_abs_e[1] = new Double_t[num_peaks[1]]{0.0035, 0.00005, 0.00003, 0.0005, 0.0013, 0., 0.0006}; //ba133
  source_abs_e[2] = new Double_t[num_peaks[2]]{0.00846261, 0.0004, 0.00004, 0.0012, 0.00012, 0.0006, 0.00012, 0.0003, 0.0004, 0.00013, 0.00012, 0.00021,  0.000020, 0.00005, 0.000015}; //co56
  source_abs_e[3] = new Double_t[num_peaks[3]]{0.0003, 0.000006}; 	//co60
  source_abs_e[4] = new Double_t[num_peaks[4]]{0.005, 0.007, 0.0003, 0.003, 0.0003};	//bi207

  // Custom fit ranges, if 0 uses a default of 10
  Int_t** fit_min_range = new Int_t*[5];
  fit_min_range[0] = new Int_t[num_peaks[0]]{0, 6, 8, 0, 0, 0, 6, 0, 8, 10}; //eu152
  fit_min_range[1] = new Int_t[num_peaks[1]]{10,5,0,0,0, 0, 0}; //ba133
  fit_min_range[2] = new Int_t[num_peaks[2]]{0, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 20, 10, 20, 20}; //co56
  fit_min_range[3] = new Int_t[num_peaks[3]]{15, 15}; //co60
  fit_min_range[4] = new Int_t[num_peaks[4]]{8, 8, 0, 15, 18}; //bi207

  Int_t** fit_max_range = new Int_t*[5];
  fit_max_range[0] = new Int_t[num_peaks[0]]{0, 6, 5, 0, 0, 0, 0, 5, 6, 0, 10}; //eu152
  fit_max_range[1] = new Int_t[num_peaks[1]]{4,5, 0, 0, 0, 0,5}; //ba133
  fit_max_range[2] = new Int_t[num_peaks[2]]{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 20, 0, 0}; //co56
  fit_max_range[3] = new Int_t[num_peaks[3]]{0, 0}; //co60
  fit_max_range[4] = new Int_t[num_peaks[4]]{6, 6}; //bi207

  // Containers for things
  Double_t** peak_areas = new Double_t*[5];
  peak_areas[0] = new Double_t[num_peaks[0]]{}; //eu152
  peak_areas[1] = new Double_t[num_peaks[1]]{}; //ba133
  peak_areas[2] = new Double_t[num_peaks[2]]{}; //co56
  peak_areas[3] = new Double_t[num_peaks[3]]{}; //co60
  peak_areas[4] = new Double_t[num_peaks[4]]{}; //bi207

  Double_t** peak_errors = new Double_t*[5];
  peak_errors[0] = new Double_t[num_peaks[0]]{}; //eu152
  peak_errors[1] = new Double_t[num_peaks[1]]{}; //ba133
  peak_errors[2] = new Double_t[num_peaks[2]]{}; //co56
  peak_errors[3] = new Double_t[num_peaks[3]]{}; //co60
  peak_errors[4] = new Double_t[num_peaks[4]]{}; //bi207

  Double_t** eff_cor_areas = new Double_t*[5];
  eff_cor_areas[0] = new Double_t[num_peaks[0]]{}; //eu152
  eff_cor_areas[1] = new Double_t[num_peaks[1]]{}; //ba133
  eff_cor_areas[2] = new Double_t[num_peaks[2]]{}; //co56
  eff_cor_areas[3] = new Double_t[num_peaks[3]]{}; //co60
  eff_cor_areas[4] = new Double_t[num_peaks[4]]{}; //bi207

  Double_t** eff_cor_errors = new Double_t*[5];
  eff_cor_errors[0] = new Double_t[num_peaks[0]]{}; //eu152
  eff_cor_errors[1] = new Double_t[num_peaks[1]]{}; //ba133
  eff_cor_errors[2] = new Double_t[num_peaks[2]]{}; //co56
  eff_cor_errors[3] = new Double_t[num_peaks[3]]{}; //co60
  eff_cor_errors[4] = new Double_t[num_peaks[4]]{}; //bi207

  Double_t** addback_areas = new Double_t*[5];
  addback_areas[0] = new Double_t[num_peaks[0]]{}; //eu152
  addback_areas[1] = new Double_t[num_peaks[1]]{}; //ba133
  addback_areas[2] = new Double_t[num_peaks[2]]{}; //co56
  addback_areas[3] = new Double_t[num_peaks[3]]{}; //co60
  addback_areas[4] = new Double_t[num_peaks[4]]{}; //bi207

  Double_t** addback_errors = new Double_t*[5];
  addback_errors[0] = new Double_t[num_peaks[0]]{}; //eu152
  addback_errors[1] = new Double_t[num_peaks[1]]{}; //ba133
  addback_errors[2] = new Double_t[num_peaks[2]]{}; //co56
  addback_errors[3] = new Double_t[num_peaks[3]]{}; //co60
  addback_errors[4] = new Double_t[num_peaks[4]]{}; //bi207

  Double_t** eff_cor_add_areas = new Double_t*[5];
  eff_cor_add_areas[0] = new Double_t[num_peaks[0]]{}; //eu152
  eff_cor_add_areas[1] = new Double_t[num_peaks[1]]{}; //ba133
  eff_cor_add_areas[2] = new Double_t[num_peaks[2]]{}; //co56
  eff_cor_add_areas[3] = new Double_t[num_peaks[3]]{}; //co60
  eff_cor_add_areas[4] = new Double_t[num_peaks[4]]{}; //bi207

  Double_t** eff_cor_add_errors = new Double_t*[5];
  eff_cor_add_errors[0] = new Double_t[num_peaks[0]]{}; //eu152
  eff_cor_add_errors[1] = new Double_t[num_peaks[1]]{}; //ba133
  eff_cor_add_errors[2] = new Double_t[num_peaks[2]]{}; //co56
  eff_cor_add_errors[3] = new Double_t[num_peaks[3]]{}; //co60
  eff_cor_add_errors[4] = new Double_t[num_peaks[4]]{}; //bi207

  /*************************************************************/

  // Loads histogram information
  string sourceID[10];
  long runstart[10];
  long runstop[10];
  int iii = 0;
  string histfile = "histinfo.dat";
  ifstream fp;
  fp.open(histfile.c_str());
  while (fp.good()) {
    fp >> sourceID[iii] >> runstart[iii] >> runstop[iii];
    iii++;
  }
  int num_sources = iii-1;
  if (num_sources == 0) {
    cout << "No sources inputted. Exiting program..." << endl;
    return;
  }//
  // Reads source information from sourceinfo.dat
  cout << "Loading Source Information" << endl;
  vector<string> known_sourceNames;
  vector<string> known_sourceIDs;
  vector<int> source_ref_activities;
  vector<int> source_ref_year;
  vector<int> source_ref_month;
  vector<int> source_ref_day;
  vector<int> source_ref_hour;
  vector<int> source_ref_minute;
  vector<int> source_ref_second;

  vector<double> energies;
  vector<double> energy_errors;
  vector<double> abs_eff;
  vector<double> abs_eff_errors;
  vector<double> addback_eff;
  vector<double> addback_eff_errors;
  vector<double> abs_eff_fit_areas;
  vector<double> abs_eff_fit_errors;
  vector<double> addback_eff_fit_areas;
  vector<double> addback_eff_fit_errors;

  iii = 0;
  string sourcefile = "sourceinfo.dat";
  ifstream fpsource;
  fpsource.open(sourcefile.c_str());
  string tmpString, tmpString2;
  int tmpSID[7];
  while (fpsource.good()) {
    fpsource >> tmpString >> tmpString2 >> tmpSID[0] >> tmpSID[1] >> tmpSID[2] >> tmpSID[3] >> tmpSID[4] >> tmpSID[5] >> tmpSID[6];
    known_sourceNames.push_back(tmpString);
    known_sourceIDs.push_back(tmpString2);
    source_ref_activities.push_back(tmpSID[0]);
    source_ref_year.push_back(tmpSID[1]);
    source_ref_month.push_back(tmpSID[2]);
    source_ref_day.push_back(tmpSID[3]);
    source_ref_hour.push_back(tmpSID[4]);
    source_ref_minute.push_back(tmpSID[5]);
    source_ref_second.push_back(tmpSID[6]);
  }
  cout << "Finished Loading source Information" << endl;

  int IDSearch;
  int sourcetype;
  TList *list = new TList();
  TH1F *singles_hist[num_sources][1];
  TH1F *addback_hist[num_sources][1];
  TH2F *singles_sum[num_sources][1];
  TH2F *addback_sum[num_sources][1];

  TH1D *singles_sum_pro[num_sources][1];
  TH1D *addback_sum_pro[num_sources][1];

  TH1D *singles_gated[num_sources][10];
  TH1D *addback_gated[num_sources][10];
  for(int i = 0; i < num_sources; i++ ){
    // Compares sourceID in histinfo.dat with known sources to identify source type for fitting
    IDSearch = search_array(known_sourceIDs, sourceID[i], known_sourceIDs.size());
    if(IDSearch == -1 ) {
      cout << "Source R" << sourceID[i] << " not found, aborting program" << endl;
      return;
    }
    if(strcmp(known_sourceNames.at(IDSearch).c_str(),"co60") == 0){
      cout << "60Co source R" << sourceID[i] << endl;
      sourcetype = 3;
    } else if(strcmp(known_sourceNames.at(IDSearch).c_str(),"co56") == 0){
      cout << "56Co source R" << sourceID[i] << endl;
      sourcetype = 2;
    } else if(strcmp(known_sourceNames.at(IDSearch).c_str(),"eu152") == 0){
      cout << "152Eu source R" << sourceID[i] << endl;
      sourcetype = 0;
    } else if(strcmp(known_sourceNames.at(IDSearch).c_str(),"ba133") == 0){
      cout << "133Ba source R" << sourceID[i] << endl;
      sourcetype = 1;
    } else if(strcmp(known_sourceNames.at(IDSearch).c_str(),"bi207") == 0){
      cout << "207Bi source R" << sourceID[i] << endl;
      sourcetype = 4;
    }
    // Loads the efficiency histograms
    load_histograms("EffHist.root", singles_hist[i], addback_hist[i], singles_sum[i], addback_sum[i], singles_sum_pro[i], addback_sum_pro[i],i);
    list->Add(singles_hist[i][0]);
    list->Add(addback_hist[i][0]);
    list->Add(singles_sum[i][0]);
    list->Add(addback_sum[i][0]);
    list->Add(singles_sum_pro[i][0]);
    list->Add(addback_sum_pro[i][0]);
    // Makes the histograms for the sum in corrections 
    make_sum_hist(singles_sum[i], singles_gated[i], addback_sum[i], addback_gated[i], source_sum_gate[sourcetype], i, num_sum[sourcetype]);
    for(int k = 0; k < num_sum[sourcetype]; k++) {
      list->Add(singles_gated[i][k]);
      list->Add(addback_gated[i][k]);
    }
    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(500000);
    // Clears arrays if using multiple of the same source
    for(int l = 0; l < num_peaks[sourcetype]; l++) {
      peak_areas[sourcetype][l] = 0;
      addback_areas[sourcetype][l] = 0;
    }
    // Fits peaks and gets areas
    fit_peaks(singles_hist[i], source_energy[sourcetype], source_energy_er[sourcetype], num_peaks[sourcetype], peak_areas[sourcetype],peak_errors[sourcetype],fit_min_range[sourcetype],fit_max_range[sourcetype]);
    fit_peaks(addback_hist[i], source_energy[sourcetype], source_energy_er[sourcetype], num_peaks[sourcetype], addback_areas[sourcetype],addback_errors[sourcetype],fit_min_range[sourcetype],fit_max_range[sourcetype]);
    if(summingcorrections) { // Summing corrections are optional
      sumout_cor(singles_sum_pro[i], source_energy[sourcetype], num_peaks[sourcetype], peak_areas[sourcetype], peak_errors[sourcetype], fit_min_range[sourcetype],fit_max_range[sourcetype]);
      sumout_cor(addback_sum_pro[i], source_energy[sourcetype], num_peaks[sourcetype], addback_areas[sourcetype], addback_errors[sourcetype], fit_min_range[sourcetype],fit_max_range[sourcetype]);

      sumin_cor(singles_sum_pro[i], singles_gated[i], source_energy[sourcetype], num_peaks[sourcetype], peak_areas[sourcetype], peak_errors[sourcetype], sourcetype, fit_min_range[sourcetype],fit_max_range[sourcetype],source_sum_gate[sourcetype],source_sum_fit[sourcetype],num_sum[sourcetype]);
      sumin_cor(addback_sum_pro[i], addback_gated[i], source_energy[sourcetype], num_peaks[sourcetype], addback_areas[sourcetype], addback_errors[sourcetype], sourcetype, fit_min_range[sourcetype],fit_max_range[sourcetype], source_sum_gate[sourcetype], source_sum_fit[sourcetype],num_sum[sourcetype]);
    }
    // Calculates the source strength and absolute efficiency 
    int source_ref_date = date_to_unix(source_ref_year.at(IDSearch), source_ref_month.at(IDSearch), source_ref_day.at(IDSearch), source_ref_hour.at(IDSearch), source_ref_minute.at(IDSearch), source_ref_second.at(IDSearch));
    double activity = RadActivity(source_ref_date,(runstart[i]+runstop[i])/2,source_halfl[sourcetype],source_ref_activities[IDSearch]);
    double activity_error = RadActivityError(source_ref_date,(runstart[i]+runstop[i])/2,source_halfl[sourcetype],source_halfl_e[sourcetype],source_ref_activities[IDSearch]);
    abseff(num_peaks[sourcetype], peak_areas[sourcetype], peak_errors[sourcetype], eff_cor_areas[sourcetype], eff_cor_errors[sourcetype], source_abs[sourcetype], source_abs_e[sourcetype], (runstop[i] - runstart[i]), activity, activity_error);
    abseff(num_peaks[sourcetype], addback_areas[sourcetype], addback_errors[sourcetype], eff_cor_add_areas[sourcetype], eff_cor_add_errors[sourcetype], source_abs[sourcetype], source_abs_e[sourcetype], (runstop[i] - runstart[i]), activity, activity_error);
    // Adds all values to vecttors for readout later
    for(int j=0;j<num_peaks[sourcetype];j++){
      energies.push_back(abs(source_energy[sourcetype][j]));
      energy_errors.push_back(source_energy_er[sourcetype][j]);
      abs_eff.push_back(eff_cor_areas[sourcetype][j]);
      abs_eff_errors.push_back(eff_cor_errors[sourcetype][j]);
      addback_eff.push_back(eff_cor_add_areas[sourcetype][j]);
      addback_eff_errors.push_back(eff_cor_add_errors[sourcetype][j]);

      abs_eff_fit_areas.push_back(peak_areas[sourcetype][j]);
      abs_eff_fit_errors.push_back(peak_errors[sourcetype][j]);
      addback_eff_fit_areas.push_back(addback_areas[sourcetype][j]);
      addback_eff_fit_errors.push_back(addback_errors[sourcetype][j]);
    }
  }

  ofstream efftxt("eff_points.dat");
  ofstream addtxt("addBack_eff_points.dat");

  // Draws absolute efficiency curve and outputs the efficiency values to text files to be used in curve_fit.cxx
  TGraphErrors *eff_curve = new TGraphErrors();
  TGraphErrors *addback_curve = new TGraphErrors();
  int counter = 0;
  for(int i = 0; i < energies.size(); i++) {
    eff_curve->SetPoint(counter,energies.at(i),abs_eff.at(i));
    eff_curve->SetPointError(counter,energy_errors.at(i),abs_eff_errors.at(i));
    addback_curve->SetPoint(counter,energies.at(i),addback_eff.at(i));
    addback_curve->SetPointError(counter,energy_errors.at(i),addback_eff_errors.at(i));
    efftxt << energies.at(i) << "\t" << abs_eff.at(i) << "\t" << energy_errors.at(i) << "\t" << abs_eff_errors.at(i) << "\t" << abs_eff_fit_areas.at(i) << "\t" << abs_eff_fit_errors.at(i) << endl;
    addtxt << energies.at(i) << "\t" << addback_eff.at(i) << "\t" << energy_errors.at(i) << "\t" << addback_eff_errors.at(i) << "\t" << addback_eff_fit_areas.at(i) << "\t" << addback_eff_fit_errors.at(i) << endl;
    counter++;
  }
  eff_curve->SetName("Efficiency_Curve_Singles");
  eff_curve->Draw("a*");
  addback_curve->SetName("Addback_Curve_Singles");
  addback_curve->Draw("a*");
  list->Add(eff_curve);
  list->Add(addback_curve);

  TFile *out = new TFile("Fits.root", "RECREATE");
  out->cd();
  list->Write();
  out->Close();//closes files
}

