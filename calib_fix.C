Int_t num_cores = 64;

//fits a equation on the centroids vs energies for each core of the TIGRESS detector
void fit_equation(TList *list, Double_t energy[], Double_t energy_er[], Double_t** centroids, Double_t** centroids_er, Int_t num_peaks_used,
					Double_t thirds[], Double_t gains[], Double_t offsets[], char use_quad) {
	
	
	string equation;

	if (use_quad == 'Y' || use_quad == 'y') {
		equation = "[0] + [1]*x + [2]*x*x";
	} else {
		equation = "[0] + [1]*x";
	}//else

	ofstream quad_fit("quad_fit_points.txt");
		for (int i = 0; i < num_peaks_used; i++) {
		quad_fit << energy[i] << " " << energy_er[i] << endl;
	}//for

	for (int i = 0; i < num_cores; i++) {

		Double_t energy_copy[num_peaks_used];

		Double_t sample_good_energy = -1;
		Double_t sample_good_centroid = -1;
		for (int k = 0; k < num_peaks_used; k++) {
			if (centroids[i][k] != -1) {
				sample_good_energy = energy[k];
				sample_good_centroid = centroids[i][k];
			}//if
		}//for

		for (int j = 0; j < num_peaks_used; j++) {
			energy_copy[j] = energy[j];

			quad_fit << centroids[i][j] << " " << centroids_er[i][j] << endl;

			if (centroids[i][j] == -1) {
				energy_copy[j] = sample_good_energy;
				centroids[i][j] = sample_good_centroid;
			}//if
		}//for

		//if this centroid is null, then skip this core	
		if (centroids[i][0] == -1 || centroids_er[i][0] == -1) {
			thirds[i] = 0;
			gains[i] = 0;
			offsets[i] = -1;
			continue;
		}//if

		//creating graph of centroids vs energy
		TGraphErrors *gr = new TGraphErrors(num_peaks_used, centroids[i], energy_copy, centroids_er[i], energy_er);
		gr->SetName(Form("Equation Fit - %i", i));
		gr->Draw("AP");

		list->Add(gr);

		//fitting quadratic on data points
		
		TF1 *coeffFit = new TF1("coeffFit", &equation[0]);
		gr->Fit(coeffFit, "Q+");

		if (use_quad == 'Y' || use_quad == 'y') {		
			thirds[i] = coeffFit->GetParameter(2);
		} else {
			thirds[i] = 0;
		}//else


		//saving curve parameters		
		gains[i] = coeffFit->GetParameter(1);
		offsets[i] = coeffFit->GetParameter(0);

		double x[num_peaks_used];
		double resid[num_peaks_used];
		
		for (int j = 0; j < num_peaks_used; j++) {
	      x[j] = j;
		  resid[j] = energy[j] - coeffFit->Eval(centroids[i][j]);  
		}//for
		
		TGraph *resid_g = new TGraph(num_peaks_used, x, resid);
		resid_g->SetName(Form("Fit Residuals - %i", i));
	  
	  	list->Add(resid_g);	

	}//for

	cout << __LINE__ << endl;

	quad_fit.close();
}//fit_quadratic

Int_t search_closest(Int_t num_peaks, Double_t centroids[], Double_t val) {

	if (val < centroids[0]) {
		return 0;
	}//if

	for (int i = 0; i < num_peaks - 1; i++) {		
		if (val == centroids[i]) {
			return i;
		}//if	

		if (val > centroids[i] && val < centroids[i + 1]
			&& TMath::Abs(centroids[i] - val) < TMath::Abs(centroids[i + 1] - val)) {
			return i;
		} else if (val > centroids[i] && val < centroids[i + 1]
			&& TMath::Abs(centroids[i] - val) > TMath::Abs(centroids[i + 1] - val)) {
			return i + 1;
		}//elseif
	}//for

	return num_peaks - 1;

}//search_closest

void calib_fix() {

	char rewrite_root_file;
	cout << "Rewrite quad_cal.root? (Y/N) ";
	cin >> rewrite_root_file;

	char use_quad;
	cout << "Use quadratic fit equation? (Y/N) ";
	cin >> use_quad;

	Int_t num_peaks;	

	cout << "How many peaks were used in the calibration process?" << endl;
	cin >> num_peaks;

	ifstream data_points ("quad_fit_points.txt");
	
	Double_t **centroids = new Double_t*[num_cores];
	Double_t **centroids_er = new Double_t*[num_cores];

	for (int i = 0; i < num_cores; i++) {
		centroids[i] = new Double_t[num_peaks];
		centroids_er[i] = new Double_t[num_peaks];
	}//for

	Double_t energy[num_peaks];
	Double_t energy_er[num_peaks];

	for (int i = 0; i < num_peaks; i++) {
		data_points >> energy[i] >> energy_er[i];
	}//for

	for (int i = 0; i < num_cores; i++) {
		for (int j = 0; j < num_peaks; j++) {
			data_points >> centroids[i][j] >> centroids_er[i][j];
		}//for
	}//for
	

	cout << "Enter the graph number, the centroid of the old peak, and its edited centroid." << endl;
	while (true) {
		Int_t graph_num;
		Double_t old_centr;
		Double_t new_centr;

		cout << "Graph Number: ";
		cin >> graph_num;
		cout << "Old Centroid Approx: ";
		cin >> old_centr;
		cout << "New Centroid Loc: ";
		cin >> new_centr;
	
		if (graph_num == -1) {
			break;
		}//if

		cout << "Confirm Change " << centroids[graph_num][search_closest(num_peaks, centroids[graph_num], old_centr)] << " to " << new_centr << "? (Y/N) ";
		string answer;
		cin >> answer;
		cout << "--------------" << endl;
	
		if (answer == "y" || answer == "Y") {	
			centroids[graph_num][search_closest(num_peaks, centroids[graph_num], old_centr)] = new_centr;
			centroids_er[graph_num][search_closest(num_peaks, centroids[graph_num], old_centr)] = 0.5;	
		}//if

	}//while
	
	Double_t quad_thirds[num_cores];
	Double_t quad_gains[num_cores];
	Double_t quad_offsets[num_cores];
	
	TList *list = new TList;

	fit_equation(list, energy, energy_er, centroids, centroids_er, num_peaks, quad_thirds, quad_gains, quad_offsets, use_quad);
	
	string outfilepath = "quad_fit.root";
	
	if (rewrite_root_file == 'Y' || rewrite_root_file == 'y') {

		outfilepath = "quad_cal.root";

		TFile *infile = new TFile("quad_cal.root", "READ");
		
		for (int i = 0; i < num_cores; i++) {
			for (int j = 1; ; j++) {

				TObject *obj = infile->Get(Form("Histogram - %i;%i", i, j));

				if (obj != 0) {
					((TH1F *) obj)->SetDirectory(0);
					list->Add((TH1F *)obj);
				} else {
					break;
				}//else

			}//for

		}//for

		infile->Close();
	}//if


	TFile *out = new TFile(&outfilepath[0], "RECREATE");
	out->cd();
	list->Write();
	out->Close();
	
	ofstream quad_out("quad_energy_coeff.txt");

	quad_out << "float non_lin[" << num_cores << "] = {" << quad_thirds[0];

	for (int i = 1; i < num_cores; i++) {
		quad_out << ", " << quad_thirds[i];
	}//for

	quad_out << "};" << endl;
	
	quad_out << "float gain[" << num_cores << "] = {" << quad_gains[0];

	for (int i = 1; i < num_cores; i++) {
		quad_out << ", " << quad_gains[i];
	}//for

	quad_out << "};" << endl;

	quad_out << "float offset[" << num_cores << "] = {" << quad_offsets[0];

	for (int i = 1; i < num_cores; i++) {
		quad_out << ", " << quad_offsets[i];
	}//for

	quad_out << "};" << endl;
	
	quad_out.close();

}//calib_fix


