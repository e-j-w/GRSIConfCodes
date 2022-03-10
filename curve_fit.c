//g++ curve_fit.c -std=c++0x `root-config --cflags --libs` -lTreePlayer -lMathMore -lSpectrum -lMinuit -lPyROOT -o EffCurve

#include <iostream>
#include <iomanip>
#include "TH1.h"
#include "TF1.h"
#include "TTree.h"
#include "TChain.h"
#include "TH2.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TFile.h"
#include "TMath.h"
#include "TFitResult.h"
#include <cmath>
#include "Math/PdfFuncMathCore.h"
#include "Math/MinimizerOptions.h"
//#include "/opt/root_v6.14.06/include/Math/MinimizerOptions.h"
using namespace std;

Double_t RadwareEfficiency(Double_t *dim, Double_t *par) {

	if (dim[0] < 1) {
		return 0;
	}//if

	Double_t x = dim[0];

	Double_t A = par[0];
	Double_t B = par[1];
	Double_t C = par[2];
	Double_t D = par[3];
	Double_t E = par[4];
	Double_t F = par[5];
	Double_t G = par[6];
	Double_t H = par[7];

	return H * exp(pow((pow(A + B * log(x / 100.) + C * log(x / 100.) * log(x / 100.), -G) + pow(D + E * log(x / 1000.) + F * log(x / 1000.) * log(x / 1000.), -G)), -1/G));
}
Double_t CarlottaEfficiency(Double_t *dim, Double_t *par) {

        if (dim[0] < 1) {
                return 0;
        }//if

        Double_t x = dim[0];

        Double_t A = par[0];
        Double_t B = par[1];
        Double_t C = par[2];
        Double_t D = par[3];
        Double_t E = par[4];
        Double_t F = par[5];

        return exp(A + B*x*0.05 + C*log(x*0.05)/(x*0.05)+ D/(x*0.05) + E/pow((x*0.05),2) + F/pow((x*0.05),3));
}

Double_t PolEfficiency(Double_t *dim, Double_t *par) {
  Double_t sum = 0;
  Double_t x = dim[0];

  for(int i = 0; i < 8; i++) {
    sum += par[i]*pow(log(x),i);
  }
  return log(sum);
}

Double_t crystal_ball(Double_t *dim, Double_t *par) {

  Double_t x = -dim[0];
  Double_t alpha = par[0];
  Double_t n = par[1];
  Double_t sigma = par[2];
  Double_t H = par[3];
  Double_t mu = par[4];

  return H*ROOT::Math::crystalball_function(x, alpha, n, sigma, mu);
}

Double_t (*func_ptr[4]) (Double_t *dim, Double_t *par) = {RadwareEfficiency, CarlottaEfficiency, PolEfficiency, crystal_ball};

//void curve_fit(float energy = 0) {
void curve_fit(float energy[] = {}) {
  TList *list = new TList();
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(500000);
  int Npar[4] = {8, 6, 8, 5};
  TF1 *eff_functions[2][4];
  string name[2] = {"Singles", "Addback"};
  TGraphErrors *effcurve[2];
  effcurve[0] = new TGraphErrors("eff_points.dat");
  effcurve[1] = new TGraphErrors("addBack_eff_points.dat");

  double Xmin = 0;
  double Xmax = 0;
  double tmpX, tmpY;
  for(int i = 0; i < effcurve[0]->GetN(); i++) {
    if(i == 0) {
      effcurve[0]->GetPoint(i, tmpX, tmpY);
      Xmin = tmpX;
      Xmax = tmpX;
    } else {
      effcurve[0]->GetPoint(i, tmpX, tmpY);
      if(Xmin > tmpX) Xmin = tmpX;
      if(Xmax < tmpX) Xmax = tmpX;
    }
  }

  double par[4][8] = {
    {7.04,0.70,0,5.273,-0.863,0.01,11,1.25e-4},
    {-3.84,-0.0046,16,-23,43,-28},
    {18.4215,19.9859,10,5.81387,-0.520934,-0.0162382,2.01455,0.000141736},
    {0.6, 5, 10, 0.2, 150}
  };

  double par_lower[4][8] = {
    {5.0, 0.0, 0.0, 5.0, -5, -1 , -1 , 1.0e-6},
    {-10.0, -1, -5, -50, -100 , -100},
    {10.0, 10.0 , 0, 0, -2, -0.1, -1 , 0.00004},
    {-10.0, 0.1, 1, 0.05, -250}
  };

  double par_upper[4][8] = {
    {25.0, 20.0, 0.0, 20.0, 5, 5 , 10, 1.0e-3},
    {10.0, 1.0, 25, 50.0, 100, 100.0},
    {30, 30, 10, 10, 0, 0 , 5, 0.0005},
    {10, 10, 500, 1.0, 250}
  };

  double err[1];
  double x[1] = {energy[0]};
  vector<double> eff[2];
  vector<double> ef_er[2];

  for(int k = 0; k < 2; k++) {
    char gname[64];
    sprintf(gname,"%s_Curve", name[k].c_str());
    effcurve[k]->SetName(gname);
    list->Add(effcurve[k]);
    for(int i = 0; i < 4; i++) {
      char fname[64];
      sprintf(fname,"%s_func_%i", name[k].c_str(), i);
      eff_functions[k][i] = new TF1(fname, *func_ptr[i], Xmin-5, Xmax+5, Npar[i]);
      for(int j = 0; j < Npar[i]; j++) {
        eff_functions[k][i]->SetParameter(j, par[i][j]);
        eff_functions[k][i]->SetParLimits(j, par_lower[i][j], par_upper[i][j]);
      };
      TFitResultPtr r = effcurve[k]->Fit(eff_functions[k][i],"QSR");
      if(energy[0] > 1) {
        r->GetConfidenceIntervals(1, 1, 1, x, err, 0.683, false);
	eff[k].push_back(eff_functions[k][i]->Eval(x[0]));
        ef_er[k].push_back(err[0]);
      }

      list->Add(eff_functions[k][i]);
    }
  }

  if(energy[0] > 1) {
    cout << energy[0] << "\n";
    for(int i = 0; i < 2; i++) {
      cout << name[i] << "\n";
      for(int j = 0; j < 4; j++) {
        cout << eff[i].at(j)  << " \u00B1 " << ef_er[i].at(j) << endl;
      }
    }
  }
  TFile *outfile = new TFile("curves.root","RECREATE");
  outfile->cd();
  list->Write();
  outfile->Close();
}

int main(int argc, char **argv){

  float aa[argc];
  if(argc == 1) {
    aa[0] = 0.0;
    curve_fit(aa);
  } else if (argc == 2) {
    for(int i = 1; i < argc; i++) aa[i-1] = atof(argv[i]);
    curve_fit(aa);
//    curve_fit(atof(argv[1]));
  }
  return 0;
}

